

###########
library(gdata)
library(datasets)
library(nlme)
library(TSA)
library(maxLik)
library(boot)
library(turner)
library(pracma)
library(dplyr)
library(purrr)
library(ggplot2)


##################################################
###make sure the data frame just contain the following variables and it is followed by this order:
## exposure, mediator,covariates, time, event.###########


##create the big function###
##bootstrap data: The original data that will be used to estimate the natural indirect and direct effects. It also
##                will be used in the bootstrap to have the confidence interval.

##num: how many times to run in the bootstrap to create the confidence interval.
##seed: for the bootstrap. Default is 1.
##cov: The baseline covariates you want to specify when it calculates the natural indirect and direct effects,
##     The default is a vector of 1.


one_mediator_bootstrap=function(bootstrap_data,num,cov=c(rep(1,ncol(bootstrap_data)-4)),seed=1){

  set.seed(seed)


  ###create the function to

  ##1:estimate the parameters of Cox model
  ##2 calculate baseline hazard and cumulative baseline hazard for each event time points
  ###in bootstrap######
  cox_est=function( data, indices){

    dat1=data[indices,]
    dat1=dat1[,rev(names(dat1))]
    dat1=dat1[order(dat1[,2],dat1[,1]),]
    dat1$d=0
    dat1$zlag_t=zlag(dat1[,2])
    dat1$zlag_t[is.na(dat1$zlag_t)] = 0

    for (i in 1:nrow(dat1)){

      if (dat1[i,1]==1){

        if (dat1[i,2]!=dat1$zlag_t[i]){
          dat1$d[i]=1
        }
        else {dat1$d[i]=dat1$d[i-1]+1}
      }
      else {dat1$d[i]=0}

    }


    dat1=dat1[order(dat1[,2], -dat1[,1]),]

    dat1$d_first=ifelse((dat1[,1]==1 & !duplicated(dat1[,2])),1,0)

    dat1=dat1[order(dat1[,2], dat1$d),]

    dat1$d_last=ifelse((dat1[,1]==1 & !duplicated(dat1[,2],fromLast=TRUE)),1,0)

    ##########estimate################################################
    dat1=dat1[order(-dat1[,2], dat1$d),]


    dat4=cbind((dat1[,rev(names(dat1[1:(ncol(dat1)-4)]))]),dat1[(ncol(dat1)-3):ncol(dat1)])

    #convert the dataframe to matrix######


    l_cox=function(par,dat4){

      dat4=data.matrix(dat4,rownames.force = NA)

      ll=sum((dat4[,ncol(dat4)-4])*((dat4[,1])*par[1]
                                    +(dat4[,2])*par[2]
                                    +(dat4[,1]*dat4[,2])*par[3]
                                    +dat4[,3:(ncol(dat4)-6)]%*%par[4:(ncol(dat4)-5)])
             -(dat4[,"d_last"])*(log((cumsum(exp((dat4[,1])*par[1]
                                                 + (dat4[,2])*par[2]
                                                 + (dat4[,1]*dat4[,2])*par[3]
                                                 +dat4[,3:(ncol(dat4)-6)]%*%par[4:(ncol(dat4)-5)])))^(dat4[,"d"])))

             -log(par[ncol(dat4)+1])-0.5*(((dat4[,2]-
                                              (par[ncol(dat4)-4]+(dat4[,1]*par[ncol(dat4)-3])
                                               +dat4[,3:(ncol(dat4)-6)]%*%par[(ncol(dat4)-2):(ncol(dat4))]))/par[ncol(dat4)+1])^2))


      return(ll)

    }

    h44= maxLik(l_cox, dat4 = dat4,start = c(0,0,0,rep(0,(ncol(dat4)-8)),0,0,rep(0,(ncol(dat4)-8)),1))

    h44_est=h44$estimate
    h44_est= as.vector(h44_est)

    ####!!!!!!I prefer using  h44_est as matrix in the following calculation!!!!!!!!!!!!!!!######

    t_last=ifelse(( !duplicated(dat4[,ncol(dat4)-5],fromLast=TRUE)),1,0)

    dat4=cbind(dat4,t_last)
    dat4=data.matrix(dat4,rownames.force = NA)

    hazard= dat4[,"t_last"]*( dat4[,"d"]/(cumsum(exp((dat4[,1])* h44_est[1]
                                                     +(dat4[,2])* h44_est[2]
                                                     +(dat4[,1]*dat4[,2])* h44_est[3]
                                                     +dat4[,3:(ncol(dat4)-7)]%*% h44_est[4:(ncol(dat4)-6)]))))


    dat4=cbind(dat4,hazard)

    dat4=dat4[order(dat4[,(ncol(dat4)-7)]),]

    cum_hazard=cumsum(dat4[,"hazard"])

    dat4=cbind(dat4[,-c(ncol(dat4)-6,ncol(dat4)-4,ncol(dat4)-3,ncol(dat4)-2,ncol(dat4)-1)]
               ,cum_hazard=cumsum(dat4[,"hazard"]),
               t(replicate(nrow(dat4), h44_est)))


    return(dat4)



  }
  ####end of create the function for bootstrap######
  #####run the bootstrap and organize the results. Prepare for the integration####
  results = boot(bootstrap_data, statistic=cox_est,R=num)
  results_t=results$t
  results_t0=results$t0

  results_t0=results_t0[results_t0[,ncol(bootstrap_data)+1] != 0,]

  t_results_t=t(results_t)
  row_sets=c(rep(nrow(bootstrap_data),(3*ncol(bootstrap_data))))
  blocks=matrix_to_blocks(t_results_t, row_sets)
  un_t=unlist(blocks)

  restructure_t=matrix(un_t,nrow=(3*ncol(bootstrap_data)),byrow=T)

  restructure_t=t(restructure_t)
  restructure_t=restructure_t[restructure_t[,ncol(bootstrap_data)+1] != 0,]
  ###end of preparation ###
  ##############VanderWeele estimate###############
  indirect_vander_t0=exp((results_t0[,(ncol(bootstrap_data)+4)]* results_t0[,(2*ncol(bootstrap_data)+3)]
                          + results_t0[,(ncol(bootstrap_data)+5)]* results_t0[,(2*ncol(bootstrap_data)+3)]*1)*(1-0))

  direct_vander_t0=exp((results_t0[,(ncol(bootstrap_data)+3)]
                        + results_t0[,(ncol(bootstrap_data)+5)]*(results_t0[,(2*ncol(bootstrap_data)+2)]
                                                                 + results_t0[,(2*ncol(bootstrap_data)+3)]*0
                                                                 + rowSums(results_t0[,(2*ncol(bootstrap_data)+4):(3*ncol(bootstrap_data)-1)]%*%cov)
                                                                 + results_t0[,(ncol(bootstrap_data)+4)]*(results_t0[,(3*ncol(bootstrap_data))])^2))*(1-0)
                       +0.5*(results_t0[,(ncol(bootstrap_data)+5)])^2
                       *( results_t0[,(3*ncol(bootstrap_data))])^2*(1-0) )





  indirect_vander_t=exp((restructure_t[,(ncol(bootstrap_data)+4)]* restructure_t[,(2*ncol(bootstrap_data)+3)]
                         + restructure_t[,(ncol(bootstrap_data)+5)]* restructure_t[,(2*ncol(bootstrap_data)+3)]*1)*(1-0))

  direct_vander_t=exp((restructure_t[,(ncol(bootstrap_data)+3)]
                       + restructure_t[,(ncol(bootstrap_data)+5)]*(restructure_t[,(2*ncol(bootstrap_data)+2)]
                                                                   + restructure_t[,(2*ncol(bootstrap_data)+3)]*0
                                                                   + rowSums(restructure_t[,(2*ncol(bootstrap_data)+4): (3*ncol(bootstrap_data)-1)]%*%cov)
                                                                   + restructure_t[,(ncol(bootstrap_data)+4)]*(restructure_t[,(3*ncol(bootstrap_data))])^2))*(1-0)
                      +0.5*(restructure_t[,(ncol(bootstrap_data)+5)])^2
                      *( restructure_t[,(3*ncol(bootstrap_data))])^2*(1-0) )

  #############



  ###create the function for integration######

  fun=function(dd){

    fun_a1_ma0_n=function(m)

      dd[(ncol(bootstrap_data)+1)]*exp(dd[(ncol(bootstrap_data)+3)]*1
                                       +dd[(ncol(bootstrap_data)+4)]*m
                                       +dd[(ncol(bootstrap_data)+5)]*1*m
                                       +sum(dd[(ncol(bootstrap_data)+6):(2*ncol(bootstrap_data)+1)]%*%cov)
                                       - dd[(ncol(bootstrap_data)+2)]*exp(dd[(ncol(bootstrap_data)+3)]*1
                                                                          +dd[(ncol(bootstrap_data)+4)]*m
                                                                          +dd[(ncol(bootstrap_data)+5)]*1*m
                                                                          +sum(dd[(ncol(bootstrap_data)+6):(2*ncol(bootstrap_data)+1)]%*%cov))
                                       -(m-dd[(2*ncol(bootstrap_data)+2)]
                                         -dd[(2*ncol(bootstrap_data)+3)]*0
                                         -(sum(dd[(2*ncol(bootstrap_data)+4):(3*ncol(bootstrap_data)-1)]%*%cov)))^2/(2*dd[(3*ncol(bootstrap_data))]^2))

    fun_a1_ma0_d=function(m)

      exp(- dd[(ncol(bootstrap_data)+2)]*exp(dd[(ncol(bootstrap_data)+3)]*1
                                             +dd[(ncol(bootstrap_data)+4)]*m
                                             +dd[(ncol(bootstrap_data)+5)]*1*m
                                             +sum(dd[(ncol(bootstrap_data)+6):(2*ncol(bootstrap_data)+1)]%*%cov))
          -(m-dd[(2*ncol(bootstrap_data)+2)]
            -dd[(2*ncol(bootstrap_data)+3)]*0
            -(sum(dd[(2*ncol(bootstrap_data)+4):(3*ncol(bootstrap_data)-1)]%*%cov)))^2/(2*dd[(3*ncol(bootstrap_data))]^2))

    ##############################################################################

    fun_a0_ma0_n=function(m)

      dd[(ncol(bootstrap_data)+1)]*exp(dd[(ncol(bootstrap_data)+3)]*0
                                       +dd[(ncol(bootstrap_data)+4)]*m
                                       +dd[(ncol(bootstrap_data)+5)]*0*m
                                       +sum(dd[(ncol(bootstrap_data)+6):(2*ncol(bootstrap_data)+1)]%*%cov)
                                       - dd[(ncol(bootstrap_data)+2)]*exp(dd[(ncol(bootstrap_data)+3)]*0
                                                                          +dd[(ncol(bootstrap_data)+4)]*m
                                                                          +dd[(ncol(bootstrap_data)+5)]*0*m
                                                                          +sum(dd[(ncol(bootstrap_data)+6):(2*ncol(bootstrap_data)+1)]%*%cov))
                                       -(m-dd[(2*ncol(bootstrap_data)+2)]
                                         -dd[(2*ncol(bootstrap_data)+3)]*0
                                         -(sum(dd[(2*ncol(bootstrap_data)+4):(3*ncol(bootstrap_data)-1)]%*%cov)))^2/(2*dd[(3*ncol(bootstrap_data))]^2))

    fun_a0_ma0_d=function(m)

      exp(- dd[(ncol(bootstrap_data)+2)]*exp(dd[(ncol(bootstrap_data)+3)]*0
                                             +dd[(ncol(bootstrap_data)+4)]*m
                                             +dd[(ncol(bootstrap_data)+5)]*0*m
                                             +sum(dd[(ncol(bootstrap_data)+6):(2*ncol(bootstrap_data)+1)]%*%cov))
          -(m-dd[(2*ncol(bootstrap_data)+2)]
            -dd[(2*ncol(bootstrap_data)+3)]*0
            -(sum(dd[(2*ncol(bootstrap_data)+4):(3*ncol(bootstrap_data)-1)]%*%cov)))^2/(2*dd[(3*ncol(bootstrap_data))]^2))
    ##############################################################################

    fun_a1_ma1_n=function(m)

      dd[(ncol(bootstrap_data)+1)]*exp(dd[(ncol(bootstrap_data)+3)]*1
                                       +dd[(ncol(bootstrap_data)+4)]*m
                                       +dd[(ncol(bootstrap_data)+5)]*1*m
                                       +sum(dd[(ncol(bootstrap_data)+6):(2*ncol(bootstrap_data)+1)]%*%cov)
                                       - dd[(ncol(bootstrap_data)+2)]*exp(dd[(ncol(bootstrap_data)+3)]*1
                                                                          +dd[(ncol(bootstrap_data)+4)]*m
                                                                          +dd[(ncol(bootstrap_data)+5)]*1*m
                                                                          +sum(dd[(ncol(bootstrap_data)+6):(2*ncol(bootstrap_data)+1)]%*%cov))
                                       -(m-dd[(2*ncol(bootstrap_data)+2)]
                                         -dd[(2*ncol(bootstrap_data)+3)]*1
                                         -(sum(dd[(2*ncol(bootstrap_data)+4):(3*ncol(bootstrap_data)-1)]%*%cov)))^2/(2*dd[(3*ncol(bootstrap_data))]^2))

    fun_a1_ma1_d=function(m)

      exp(- dd[(ncol(bootstrap_data)+2)]*exp(dd[(ncol(bootstrap_data)+3)]*1
                                             +dd[(ncol(bootstrap_data)+4)]*m
                                             +dd[(ncol(bootstrap_data)+5)]*1*m
                                             +sum(dd[(ncol(bootstrap_data)+6):(2*ncol(bootstrap_data)+1)]%*%cov))
          -(m-dd[(2*ncol(bootstrap_data)+2)]
            -dd[(2*ncol(bootstrap_data)+3)]*1
            -(sum(dd[(2*ncol(bootstrap_data)+4):(3*ncol(bootstrap_data)-1)]%*%cov)))^2/(2*dd[(3*ncol(bootstrap_data))]^2))

    ####################################################################



    a1_ma0_n=integrate(fun_a1_ma0_n, -Inf, Inf)
    a1_ma0_d=integrate(fun_a1_ma0_d, -Inf, Inf)

    a0_ma0_n=integrate(fun_a0_ma0_n, -Inf, Inf)
    a0_ma0_d=integrate(fun_a0_ma0_d, -Inf, Inf)

    a1_ma1_n=integrate(fun_a1_ma1_n, -Inf, Inf)
    a1_ma1_d=integrate(fun_a1_ma1_d, -Inf, Inf)


    a1_ma0_n_value=a1_ma0_n$value
    a1_ma0_d_value=a1_ma0_d$value

    a0_ma0_n_value=a0_ma0_n$value
    a0_ma0_d_value=a0_ma0_d$value

    a1_ma1_n_value=a1_ma1_n$value
    a1_ma1_d_value=a1_ma1_d$value

    a1_ma0=(a1_ma0_n_value/a1_ma0_d_value)
    a0_ma0=(a0_ma0_n_value/a0_ma0_d_value)
    a1_ma1=(a1_ma1_n_value/a1_ma1_d_value)

    indirect=a1_ma1/a1_ma0
    direct=a1_ma0/a0_ma0

    effects=cbind(indirect,direct)

    return(effects)
  }

  effects_for_t0=apply(results_t0,1, fun)
  effects_for_t0=t(effects_for_t0)
  rownames(effects_for_t0)<-NULL

  effects_for_t=apply(restructure_t,1, fun)
  effects_for_t=t(effects_for_t)
  rownames(effects_for_t)<-NULL

  ###effects_for_t0###is the results for original dataset####
  effects_for_t0=cbind(results_t0[,(ncol(bootstrap_data)-1)],effects_for_t0,indirect_vander_t0,direct_vander_t0)
  colnames(effects_for_t0) <- c("event_time", "indirect_t0", "direct_t0", "indirect_vander_t0", "direct_vander_t0")

   ###effects_for_t is the results for bootstrap dataset######
  effects_for_t=cbind(restructure_t[,(ncol(bootstrap_data)-1)],effects_for_t,indirect_vander_t,direct_vander_t)
  colnames(effects_for_t) <- c("event_time", "indirect_t", "direct_t", "indirect_vander_t", "direct_vander_t")


  effects_for_t0=data.frame(effects_for_t0)
  effects_for_t=data.frame(effects_for_t)


  #####plot the results from bootstrap. It contains the mean of estimates of natural indirect and direct effects,and confidence intervals
  ###compared with VanderWeele's method###

  effects_for_t_cl=effects_for_t %>%
    group_by(event_time) %>%
    summarise(indirect_t_2.5 = quantile(indirect_t, probs = 0.025),
              indirect_t_97.5 = quantile(indirect_t, probs = 0.975),
              direct_t_2.5 = quantile(direct_t, probs = 0.025),
              direct_t_97.5 = quantile(direct_t, probs = 0.975))


  indirect_vander_t_cl= quantile(effects_for_t$indirect_vander_t, probs = c(0.025,0.975))
  indirect_vander_t_cl=data.frame(indirect_vander_t_cl)
  indirect_vander_t_cl=t(indirect_vander_t_cl)
  colnames(indirect_vander_t_cl)=c("indirect_vander_t_2.5","indirect_vander_t_97.5")
  rownames(indirect_vander_t_cl)=NULL

  direct_vander_t_cl= quantile(effects_for_t$direct_vander_t, probs = c(0.025,0.975))
  direct_vander_t_cl=data.frame(direct_vander_t_cl)
  direct_vander_t_cl=t(direct_vander_t_cl)
  colnames(direct_vander_t_cl)=c("direct_vander_t_2.5","direct_vander_t_97.5")
  rownames(direct_vander_t_cl)=NULL

  effects_for_t_mean=effects_for_t %>%
    group_by(event_time) %>%
    dplyr::summarize(indirect_t_mean = mean(indirect_t, na.rm=TRUE),
                     direct_t_mean = mean(direct_t, na.rm=TRUE) )

  indirect_vander_t_m=mean(effects_for_t$indirect_vander_t)
  indirect_vander_t_m=data.frame(indirect_vander_t_m)
  colnames(indirect_vander_t_m)=c("indirect_vander_t_mean")

  direct_vander_t_m=mean(effects_for_t$direct_vander_t)
  direct_vander_t_m=data.frame(direct_vander_t_m)
  colnames(direct_vander_t_m)=c("direct_vander_t_mean")


  effects_for_t_mean_cl=cbind(effects_for_t_mean,
                              indirect_vander_t_m,direct_vander_t_m,
                              effects_for_t_cl[,2:5],
                              indirect_vander_t_cl,
                              direct_vander_t_cl)

  effects_for_t_mean_cl=data.frame(effects_for_t_mean_cl)

  indirect_plot=ggplot(effects_for_t_mean_cl, aes(x=event_time)) +
    geom_line(aes(y = indirect_t_mean), color = "blue",size=0.5) +
    geom_line(aes(y = indirect_t_2.5), color = "blue",linetype="twodash")+
    geom_line(aes(y = indirect_t_97.5), color = "blue",linetype="twodash")+
    geom_line(aes(y = indirect_vander_t_mean), color="cyan",size=0.5)+
    geom_line(aes(y = indirect_vander_t_2.5), color="cyan",linetype="twodash")+
    geom_line(aes(y = indirect_vander_t_97.5), color="cyan",linetype="twodash")+
    ggtitle("Plot for indirect on hazard ratio") +
    xlab(" time") + ylab("hazard ratio")

  direct_plot=ggplot(effects_for_t_mean_cl, aes(x=event_time)) +
    geom_line(aes(y = direct_t_mean), color = "red",size=0.5) +
    geom_line(aes(y = direct_t_2.5), color = "red",linetype="twodash")+
    geom_line(aes(y = direct_t_97.5), color = "red",linetype="twodash")+
    geom_line(aes(y = direct_vander_t_mean), color="hotpink",size=0.5)+
    geom_line(aes(y = direct_vander_t_2.5), color="hotpink",linetype="twodash")+
    geom_line(aes(y = direct_vander_t_97.5), color="hotpink",linetype="twodash")+
    ggtitle("Plot for direct on hazard ratio") +
    xlab(" time") + ylab("hazard ratio")



  results_list <- list(effects_for_t0,indirect_plot,direct_plot)

  return(results_list)



}


####The end of the big function####
###the result is a list contain:

#1. The estimate of natural indirect and direct effects from original dataset
#from 1 to the last event time.
#2. The plot of natural indirect with mean and confidence interval from our approach
#compared with VanderWeele' estimate
#3.The plot of natural direct with mean and confidence interval from our approach
#compared with VanderWeele' estimate


##run the big function with example data to get the results##
bootstrap_results=one_mediator_bootstrap(example_data,500)



