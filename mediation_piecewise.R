

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
###make sure the data frame is followed by this order:
## exposure, mediator,covariates, time, event.###########



one_mediator_bootstrap_p=function(bootstrap_data,n_piece,num,seed=1,cov=c(rep(1,ncol(bootstrap_data)-4))){

  set.seed(seed)

  piecewise_est=function( data, indices){

    dat1=data[indices,]

    dat1=dat1[,rev(names(dat1))]
    dat1_event=subset(dat1, dat1[,1] == 1)
    dat1_noevent=subset(dat1, dat1[,1] ==0)
    percentile_event1= quantile(dat1_event[,2], probs=seq(0,1,by=1/n_piece),type = 3)

    percentile_all=quantile(dat1[,2],probs=1,type = 3)


    percentile_event1[1]=0
    percentile_event1[length(percentile_event1)]=percentile_all[1]

    percentile_event1=data.frame(percentile_event1)
    percentile_event1=t(percentile_event1)
    rownames(percentile_event1)=NULL
    dat1=dat1[order(dat1[,2],dat1[,1]),]

    group_name=c(paste0("group_",seq(1,n_piece,by=1)))

    group=c(rep(0,(length(group_name)*nrow(dat1))))

    dim(group)=c(nrow(dat1),length(group_name))
    group=data.frame(group)
    colnames(group)=group_name

    percentile_event=percentile_event1[rep(seq_len(nrow(percentile_event1)), nrow(dat1)),]


    group[,1] = ifelse((percentile_event[,1]<=dat1[,2] & dat1[,2]<=percentile_event[,2]), 1, 0)

    for (i  in 2:(ncol(percentile_event)-2)){

      group[,i] = ifelse(( percentile_event[,i]<dat1[,2] & dat1[,2]<=percentile_event[,i+1]), 1, 0)

    }

    group[,n_piece]=ifelse((percentile_event[,(ncol(percentile_event)-1)]<dat1[,2]),1,0)

    dat1=dat1[,rev(names(dat1))]
    #######################################
    percentile_event2=percentile_event1[-1]
    percentile_event3=percentile_event1[-length(percentile_event1)]
    percentile_event4=matrix(rep(percentile_event3,nrow(dat1)), nrow=nrow(dat1),byrow=TRUE)

    time_group=percentile_event3[-1]-percentile_event3[-length(percentile_event3)]
    time_in_group=dat1[,ncol(dat1)-1]-percentile_event4
    time_in_group=time_in_group*group
    colnames(time_in_group)=NULL
    #########################

    cum_group=c(rep(0,(n_piece-1)*nrow(dat1)))

    dim(cum_group)=c(nrow(dat1),(n_piece-1))

    cum_percentile_event=percentile_event[,-c(1,(n_piece+1))]

    time_group=t(time_group)
    time_group1=time_group[rep(seq_len(nrow(time_group)),nrow(dat1)),]


    for (i  in 1:(ncol(cum_percentile_event))){
      cum_group[,i] = ifelse(( dat1[,(ncol(dat1)-1)]>cum_percentile_event[,i]), 1, 0)
    }

    cum_group= cum_group*(time_group1)
    dat4=cbind(dat1,cum_group,time_in_group,group)
    dat4=data.matrix(dat4)
    ################################################

    l_piecewise=function(par, dat4){

      mediator_mean=par[ncol(dat1)]+(dat4[,1]*par[ncol(dat1)+1])+(dat4[,(3:(ncol(dat1)-2))])%*%(par[(ncol(dat1)+2):(2*ncol(dat1)-3)])


      ll_mediator=-log(par[2*ncol(dat1)-2])-0.5*(((dat4[,2]-mediator_mean)/par[2*ncol(dat1)-2])^2)
      cox=exp(dat4[,1]*par[1]+dat4[,2]*par[2]+dat4[,1]*dat4[,2]*par[3]+(dat4[,(3:(ncol(dat1)-2))])%*%(par[4:(ncol(dat1)-1)]))
      hazard=dat4[,((2*n_piece+ncol(dat1)):(3*n_piece-1+ncol(dat1)))]%*%par[(2*ncol(dat1)-1):(2*ncol(dat1)-2+ncol(group))]*cox
      cum_hazard=(dat4[,((ncol(dat1)+1):(ncol(dat1)+n_piece-1))]%*%par[(2*ncol(dat1)-1):(2*ncol(dat1)-3+n_piece)]
                  +(dat4[,((ncol(dat1)+n_piece):(ncol(dat1)-1+2*n_piece))])%*%par[(2*ncol(dat1)-1):(2*ncol(dat1)-2+n_piece)])*cox

      log_survival=-cum_hazard
      log_pdf=log(hazard)+log_survival


      ll_t=dat4[,ncol(dat1)]*log_pdf+(1-dat4[,ncol(dat1)])*log_survival
      ll=sum(ll_t+ll_mediator)

      return(-ll)

    }


    #  res1 <- optim(c(0,0,0,rep(0,(ncol(bootstrap_data)-4)),0,0,rep(0,(ncol(bootstrap_data)-4)),1,rep(1,n_piece)), l_piecewise, dat4=dat4, method = "BFGS", control=list(reltol=1e-8,maxit=300))
    #
    # AIC=2*(res1$value)+2*length(res1$par)
    # para=data.matrix(res1$par)

    # res1 <- optim(c(0,0,0,rep(0,(ncol(dat)-4)),0,0,rep(0,(ncol(dat)-4)),1,rep(1,n_piece)), l_piecewise, dat4=dat4, method = "BFGS", control=list(reltol=1e-8,maxit=300))
    res1=nlm(l_piecewise,c(0,0,0,rep(0,(ncol(bootstrap_data)-4)),0,0,rep(0,(ncol(dat1)-4)),1,rep(1,n_piece)),dat4=dat4,hessian = FALSE)

    #  AIC=2*(res1$value)+2*length(res1$par)
    AIC=2*res1$minimum+2*length(res1$estimate)

    #para=data.matrix(res1$par)
    para=data.matrix(res1$estimate)
    ##################################################################################


    # baseline_hazard=dat4[,((2*n_piece+ncol(dat1)):(3*n_piece-1+ncol(dat1)))]%*%para[(2*ncol(dat1)-1):(2*ncol(dat1)-2+ncol(group))]
    # base_cum_hazard=dat4[,((ncol(dat1)+1):(ncol(dat1)+n_piece-1))]%*%para[(2*ncol(dat1)-1):(2*ncol(dat1)-3+n_piece)]
    # +(dat4[,((ncol(dat1)+n_piece):(ncol(dat1)-1+2*n_piece))])%*%para[(2*ncol(dat1)-1):(2*ncol(dat1)-2+n_piece)]
    #
    # para=t(para)
    # para=para[rep(seq_len(nrow(para)), nrow(dat4)),]
    # AIC=rep(AIC,nrow(dat4))
    #
    #
    # dat2=cbind(dat4[,((ncol(dat1)-1):ncol(dat1))],para,baseline_hazard,base_cum_hazard,AIC)
    # dat2=data.matrix(dat2)
    #
    #  return(dat2)
    ####################################################################

    timeset=seq(1,quantile(bootstrap_data[,(ncol(bootstrap_data)-1)],probs=1,type = 3), by=1)
    percentile_event4=percentile_event4[1:length(timeset),]
    time_in_group=timeset-percentile_event4

    group=c(rep(0,(length(timeset)*n_piece)))
    dim(group)=c(length(timeset),n_piece)
    percentile_event=percentile_event[(1:length(timeset)),]

    group[,1] = ifelse((percentile_event[,1]<=timeset & timeset<=percentile_event[,2]), 1, 0)

    for (i  in 2:(ncol(percentile_event)-2)){

      group[,i] = ifelse(( percentile_event[,i]<timeset & timeset<=percentile_event[,i+1]), 1, 0)

    }
    group[,n_piece]=ifelse((percentile_event[,(ncol(percentile_event)-1)]<timeset),1,0)


    time_in_group=time_in_group*group


    cum_group=c(rep(0,(length(percentile_event1)-2)*length(timeset)))
    dim(cum_group)=c(length(timeset),(length(percentile_event1)-2))
    cum_percentile_event=percentile_event[,-c(1,ncol(percentile_event))]
    time_group2=time_group[rep(seq_len(nrow(time_group)),length(timeset)),]


    for (i  in 1:(ncol(cum_percentile_event))){
      cum_group[,i] = ifelse(( timeset>cum_percentile_event[,i]), 1, 0)
    }
    baseline_hazard=group%*%para[(2*ncol(bootstrap_data)-1):(2*ncol(bootstrap_data)-2+n_piece)]

    base_cum_hazard=(cum_group*(time_group2))%*%para[(2*ncol(bootstrap_data)-1):(2*ncol(bootstrap_data)-3+n_piece)]+(time_in_group%*%para[(2*ncol(bootstrap_data)-1):(2*ncol(bootstrap_data)-2+ncol(group))])

    para=t(para)
    para=para[rep(seq_len(nrow(para)), length(timeset)),]

    AIC=rep(AIC,length(timeset))

    dat2=cbind(timeset,para, baseline_hazard,base_cum_hazard,AIC)
    dat2=data.matrix(dat2)

    return(dat2)


  }


  results = boot(data=bootstrap_data, statistic=piecewise_est,R=num)
  results_t=results$t
  results_t0=results$t0



  t_results_t=t(results_t)
  row_sets=c(rep(nrow(results_t0),(2*ncol(bootstrap_data)+2+n_piece)))
  blocks=matrix_to_blocks(t_results_t, row_sets)
  un_t=unlist(blocks)

  restructure_t=matrix(un_t,nrow=(2*ncol(bootstrap_data)+2+n_piece),byrow=T)

  restructure_t=t(restructure_t)

  ##############VanderWeele estimate###############

  indirect_vander_t0=exp((results_t0[,3]* results_t0[,(ncol(bootstrap_data)+2)]
                          + results_t0[,4]* results_t0[,(ncol(bootstrap_data)+2)]*1)*(1-0))

  direct_vander_t0=exp((results_t0[,2]
                        + results_t0[,4]*(results_t0[,(ncol(bootstrap_data)+1)]
                                          + results_t0[,(ncol(bootstrap_data)+2)]*0
                                          +rowSums(results_t0[,((ncol(bootstrap_data)+3):(2*ncol(bootstrap_data)-2))]%*%cov)
                                          + results_t0[,3]*(results_t0[,(2*ncol(bootstrap_data)-1)])^2))*(1-0)
                       +0.5*(results_t0[,4])^2
                       *( results_t0[,(2*ncol(bootstrap_data)-1)])^2*(1-0) )



  indirect_vander_t=exp((restructure_t[,3]* restructure_t[,(ncol(bootstrap_data)+2)]
                         + restructure_t[,4]* restructure_t[,(ncol(bootstrap_data)+2)]*1)*(1-0))

  direct_vander_t=exp((restructure_t[,2]
                       + restructure_t[,4]*(restructure_t[,(ncol(bootstrap_data)+1)]
                                            + restructure_t[,(ncol(bootstrap_data)+2)]*0
                                            +rowSums(restructure_t[,((ncol(bootstrap_data)+3):(2*ncol(bootstrap_data)-2))]%*%cov)
                                            + restructure_t[,3]*(restructure_t[,(2*ncol(bootstrap_data)-1)])^2))*(1-0)
                      +0.5*(restructure_t[,4])^2
                      *( restructure_t[,(2*ncol(bootstrap_data)-1)])^2*(1-0) )

  #############



  fun=function(dd){

    fun_a1_ma0_n=function(m)

      dd[(2*ncol(bootstrap_data)+n_piece)]*exp(dd[2]*1
                                               +dd[3]*m
                                               +dd[4]*1*m
                                               +sum(dd[5:ncol(bootstrap_data)]%*%cov)
                                               - dd[(2*ncol(bootstrap_data)+n_piece+1)]*exp(dd[2]*1
                                                                                            +dd[3]*m
                                                                                            +dd[4]*1*m
                                                                                            +sum(dd[5:ncol(bootstrap_data)]%*%cov))
                                               -(m-dd[(ncol(bootstrap_data)+1)]
                                                 -dd[(ncol(bootstrap_data)+2)]*0
                                                 -(sum(dd[(ncol(bootstrap_data)+3):(2*ncol(bootstrap_data)-2)]%*%cov)))^2/(2*dd[(2*ncol(bootstrap_data)-1)]^2))

    fun_a1_ma0_d=function(m)

      exp(-dd[(2*ncol(bootstrap_data)+n_piece+1)]*exp(dd[2]*1
                                                      +dd[3]*m
                                                      +dd[4]*1*m
                                                      +sum(dd[5:ncol(bootstrap_data)]%*%cov))
          -(m-dd[(ncol(bootstrap_data)+1)]
            -dd[(ncol(bootstrap_data)+2)]*0
            -(sum(dd[(ncol(bootstrap_data)+3):(2*ncol(bootstrap_data)-2)]%*%cov)))^2/(2*dd[(2*ncol(bootstrap_data)-1)]^2))

    ##############################################################################

    fun_a0_ma0_n=function(m)

      dd[(2*ncol(bootstrap_data)+n_piece)]*exp(dd[2]*0
                                               +dd[3]*m
                                               +dd[4]*0*m
                                               +sum(dd[5:ncol(bootstrap_data)]%*%cov)
                                               - dd[(2*ncol(bootstrap_data)+n_piece+1)]*exp(dd[2]*0
                                                                                            +dd[3]*m
                                                                                            +dd[4]*0*m
                                                                                            +sum(dd[5:ncol(bootstrap_data)]%*%cov))
                                               -(m-dd[(ncol(bootstrap_data)+1)]
                                                 -dd[(ncol(bootstrap_data)+2)]*0
                                                 -(sum(dd[(ncol(bootstrap_data)+3):(2*ncol(bootstrap_data)-2)]%*%cov)))^2/(2*dd[(2*ncol(bootstrap_data)-1)]^2))

    fun_a0_ma0_d=function(m)

      exp(-dd[(2*ncol(bootstrap_data)+n_piece+1)]*exp(dd[2]*0
                                                      +dd[3]*m
                                                      +dd[4]*0*m
                                                      +sum(dd[5:ncol(bootstrap_data)]%*%cov))
          -(m-dd[(ncol(bootstrap_data)+1)]
            -dd[(ncol(bootstrap_data)+2)]*0
            -(sum(dd[(ncol(bootstrap_data)+3):(2*ncol(bootstrap_data)-2)]%*%cov)))^2/(2*dd[(2*ncol(bootstrap_data)-1)]^2))

    ##############################################################################

    fun_a1_ma1_n=function(m)


      dd[(2*ncol(bootstrap_data)+n_piece)]*exp(dd[2]*1
                                               +dd[3]*m
                                               +dd[4]*1*m
                                               +sum(dd[5:ncol(bootstrap_data)]%*%cov)
                                               - dd[(2*ncol(bootstrap_data)+n_piece+1)]*exp(dd[2]*1
                                                                                            +dd[3]*m
                                                                                            +dd[4]*1*m
                                                                                            +sum(dd[5:ncol(bootstrap_data)]%*%cov))
                                               -(m-dd[(ncol(bootstrap_data)+1)]
                                                 -dd[(ncol(bootstrap_data)+2)]*1
                                                 -(sum(dd[(ncol(bootstrap_data)+3):(2*ncol(bootstrap_data)-2)]%*%cov)))^2/(2*dd[(2*ncol(bootstrap_data)-1)]^2))

    fun_a1_ma1_d=function(m)

      exp(-dd[(2*ncol(bootstrap_data)+n_piece+1)]*exp(dd[2]*1
                                                      +dd[3]*m
                                                      +dd[4]*1*m
                                                      +sum(dd[5:ncol(bootstrap_data)]%*%cov))
          -(m-dd[(ncol(bootstrap_data)+1)]
            -dd[(ncol(bootstrap_data)+2)]*1
            -(sum(dd[(ncol(bootstrap_data)+3):(2*ncol(bootstrap_data)-2)]%*%cov)))^2/(2*dd[(2*ncol(bootstrap_data)-1)]^2))

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


  #effects_for_all=cbind(results_t0[ncol(bootstrap_data)-1],effects_for_all)
  effects_for_t0=cbind(results_t0[,1],effects_for_t0,indirect_vander_t0,direct_vander_t0)
  colnames(effects_for_t0) <- c("event_time", "indirect_t0", "direct_t0", "indirect_vander_t0", "direct_vander_t0")

  effects_for_t=cbind(restructure_t[,1],effects_for_t,indirect_vander_t,direct_vander_t)
  colnames(effects_for_t) <- c("event_time", "indirect_t", "direct_t", "indirect_vander_t", "direct_vander_t")


  effects_for_t0=data.frame(effects_for_t0)
  effects_for_t=data.frame(effects_for_t)

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

  AIC_t0=results_t0[1,ncol(results_t0)]

  results_list <- list(effects_for_t0,indirect_plot,direct_plot,AIC_t0)

  return(results_list)



}


bootstrap_results=one_mediator_bootstrap_p(mediation_hui1,n_piece=6,num=500,seed=1,cov=c(rep(1,ncol(mediation_hui1)-4)))

end_time <- Sys.time()

bbbb=end_time - start_time
bbbb

# results_t
#
# results_t0
