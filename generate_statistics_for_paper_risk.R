##########################################
# prepare dataset
##########################################

setwd(dir_ltla)

##########################################
# extract exposure risk
##########################################

table_outcome_riskScores <- read_csv(file=paste0("table_riskScores.",extension_csv)) %>%
  filter(maxRiskScore>100 & date>="2021-04-01" & date<"2022-03-01") 

##########################################
# general statistics
##########################################

summarystats_individual <- apply_by_ltla("table_individual",
                                         function(x){x %>% 
                                             filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive_tests=tests_positive(positive,number_exposures)) %>%
                                             summarise(rows=n(),nExp=sum(number_exposures),timeDur=sum(number_exposures*meanDuration),positive=sum(positive*number_exposures),positive_tests=sum(positive_tests))
                                         },
                                         function(y){
                                           (y %>% summarise(total_rows=sum(rows),total_nExp=sum(nExp),total_Time=sum(timeDur)/3600,positive=sum(positive),positive_tests=sum(positive_tests)))
                                         },
                                         gz=save_as_gz)
SaveTable(summarystats_individual)

summarystats_all <- apply_by_ltla("filtered_events_data_all",
                                  function(x){x %>% 
                                      filter(getmonthhq(received_date) %in% monthshq) %>% mutate(positive=0+(type=="exposureWindowPositiveTest")) %>%
                                      summarise(rows=sum(1-positive),timeDur=sum(duration_scans*(1-positive)),positive=sum(positive))
                                  },
                                  function(y){
                                    (y %>% summarise(total_rows=sum(rows),total_Time=sum(timeDur)/3600,positive=sum(positive)))
                                  },
                                  gz=save_as_gz)
SaveTable(summarystats_all)

summarystats_AAE <- AAE %>%
  filter(getmonthhq(date) %in% monthshq) %>% 
  summarise(mean_uptake=sum(users)/sum(population,na.rm=TRUE),tot_positive=sum(test_positive),tot_notifications=sum(notifications),tot_positive_after_EN=sum(positive_after_EN))
SaveTable(summarystats_AAE)

monthly_summarystats_individual <- apply_by_ltla("table_individual",
                                                 function(x){x %>% 
                                                     mutate(month=getmonthhq(notification_date)) %>% 
                                                     filter(month %in% monthshq) %>% 
                                                     mutate(positive_tests=tests_positive(positive,number_exposures)) %>%
                                                     group_by(month) %>%
                                                     summarise(rows=n(),nExp=sum(number_exposures),timeDur=sum(number_exposures*meanDuration),positive=sum(positive*number_exposures),positive_tests=sum(positive_tests))
                                                 },
                                                 function(y){
                                                   (y %>% group_by(month) %>% summarise(total_rows=sum(rows),total_nExp=sum(nExp),total_Time=sum(timeDur)/3600,positive=sum(positive),positive_tests=sum(positive_tests)))
                                                 },
                                                 gz=save_as_gz) %>%
  mutate(month=factor(month,ordered=TRUE,levels=monthshq)) %>%
  arrange(month)
SaveTable(monthly_summarystats_individual)

monthly_summarystats_all <- apply_by_ltla("filtered_events_data_all",
                                          function(x){x  %>%
                                              mutate(month=getmonthhq(received_date)) %>% 
                                              filter(month %in% monthshq) %>% 
                                              mutate(positive=0+(type=="exposureWindowPositiveTest")) %>%
                                              group_by(month) %>%
                                              summarise(rows=sum(1-positive),timeDur=sum(duration_scans*(1-positive)),positive=sum(positive))
                                          },
                                          function(y){
                                            (y %>% group_by(month) %>% summarise(total_rows=sum(rows),total_Time=sum(timeDur)/3600,positive=sum(positive)))
                                          },
                                          gz=save_as_gz) %>%
  mutate(month=factor(month,ordered=TRUE,levels=monthshq)) %>%
  arrange(month)
SaveTable(monthly_summarystats_all)

monthly_summarystats_AAE <- AAE %>%
  mutate(month=getmonthhq(date)) %>% 
  filter(month %in% monthshq) %>% 
  group_by(month) %>%
  summarise(mean_uptake=sum(users)/sum(population,na.rm=TRUE),tot_positive=sum(test_positive),tot_notifications=sum(notifications),tot_positive_after_EN=sum(positive_after_EN)) %>%
  mutate(month=factor(month,ordered=TRUE,levels=monthshq)) %>%
  arrange(month)
SaveTable(monthly_summarystats_AAE)

# correlations between quantities 
mycorq<-c("maxProximityScore","avgProximityScore","maxRiskScore","avgRiskScore","duration","cumRiskScore","bg_rate_cases_app","positive")
table_correlations_pearson <- as.data.frame(round(cor(table_outcome_riskScores[,mycorq],method = "pearson"),2))
table_correlations_spearman <- as.data.frame(round(cor(table_outcome_riskScores[,mycorq],method = "spearman"),2))
SaveTable(table_correlations_pearson)
SaveTable(table_correlations_spearman)

##########################################
# define functions 2 infer exposure risk
##########################################

printflush <- function(x){
  print(x); flush.console()
}

# inference of transmission risk
MLinference<-function(ME,MB,pos){
  nE<-dim(ME)[2]
  nB<-dim(MB)[2]
  neglogL<-function(r){
    -(
      -sum(1000*MB[!pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))-sum(ME[!pos,]%*%abs(r[1:nE]))+
        sum(log(1-exp(-(as.vector(1000*MB[pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))+as.vector(ME[pos,]%*%abs(r[1:nE]))))))
      -sum((r[1:(nE-1)]-r[2:nE])/abs(r[2:nE])*(r[2:nE]<r[1:(nE-1)]))*dim(ME)[1]*1000000
    )
  }
  mypars<-c(rep(0.001,nE),rep(0.001,nB))
  res_opt<-optim(par = mypars, fn = neglogL)$par
  for(i in 1:25){ print(i)
    old_res<-res_opt
    mypars<-c(rep(0.001,nE),rep(0.001,nB))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  printflush(res_opt)
  for(i in 1:25){ print(i)
    old_res<-res_opt
    mypars<-c(sort(rlnorm(nE,meanlog=log(0.001),sdlog=2)),rlnorm(nB,meanlog=log(0.001),sdlog=2))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  old_res<-res_opt
  res_opt<-optim(par = res_opt, fn = neglogL)$par
  printflush(res_opt)
  while(max(abs((res_opt-old_res)/res_opt))>0.01){
    old_res<-res_opt
    res_opt<-optim(par = res_opt, fn = neglogL)$par
    printflush(res_opt)
  }
  final_logL<-(-neglogL(res_opt))
  res_opt<-cbind(res_opt,res_opt,res_opt)
  res_opt[-c(1:nE),]<-res_opt[-c(1:nE),,drop=F]*1000
  res_opt[nE+(1:nB),]<-abs(res_opt[nE+(1:nB),,drop=F])
  res_opt[1:nE,]<-1-exp(-abs(res_opt[1:nE,,drop=F]))
  return(c(res_opt[,1],res_opt[,2],res_opt[,3],final_logL,length(pos)))
}

# inference of transmission risk with CIs
MLinference_profiled_nonmonotonic<-function(ME,MB,pos){
  nE<-dim(ME)[2]
  nB<-dim(MB)[2]
  neglogL<-function(r){
    -(
      -sum(1000*MB[!pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))-sum(ME[!pos,]%*%abs(r[1:nE]))+
        sum(log(1-exp(-(as.vector(1000*MB[pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))+as.vector(ME[pos,]%*%abs(r[1:nE]))))))
    )
  }
  mypars<-c(rep(0.001,nE),rep(0.001,nB))
  res_opt<-optim(par = mypars, fn = neglogL)$par
  for(i in 1:25){ print(i)
    old_res<-res_opt
    mypars<-c(rep(0.001,nE),rep(0.001,nB))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  printflush(res_opt)
  for(i in 1:25){ print(i)
    old_res<-res_opt
    mypars<-c(sort(rlnorm(nE,meanlog=log(0.001),sdlog=2)),rlnorm(nB,meanlog=log(0.001),sdlog=2))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  nvars<-length(res_opt)
  printflush(res_opt)
  old_res<-res_opt
  res_opt<-optim(par = res_opt, fn = neglogL)$par
  printflush(res_opt)
  while(max(abs((res_opt-old_res)/res_opt))>0.01){
    old_res<-res_opt
    res_opt<-optim(par = res_opt, fn = neglogL)$par
    printflush(res_opt)
  }
  final_logL<-(-neglogL(res_opt))
  res_opt<-cbind(res_opt,res_opt,res_opt)
  for(var in 1:nvars){
    printflush(var)
    temp_res_opt<-res_opt
    mask<-rep(1,nE+nB);mask[var]<-0
    neglogLprof<-function(x){neglogL(x*mask+res_opt[,2]*(1-mask))}
    factor<-1+0.1
    res_opt[var,2]<-res_opt[var,2]/factor
    if(neglogLprof(optim(par = temp_res_opt, fn = neglogLprof)$par)
       >1.92-final_logL){ res_opt[var,2]<-res_opt[var,2]*factor; factor<-1.01 } else { factor<-1.05 }
    while(
      neglogLprof(temp_res_opt<-optim(par = temp_res_opt, fn = neglogLprof)$par)
      <=1.92-final_logL){
      res_opt[var,2]<-res_opt[var,2]/factor
    }
    temp_res_opt<-res_opt
    mask<-rep(1,nE+nB);mask[var]<-0
    neglogLprof<-function(x){neglogL(x*mask+res_opt[,3]*(1-mask))}
    factor<-1+0.1
    res_opt[var,2]<-res_opt[var,2]/factor
    if(neglogLprof(optim(par = temp_res_opt, fn = neglogLprof)$par)
       >1.92-final_logL){ res_opt[var,2]<-res_opt[var,2]*factor; factor<-1.01 } else { factor<-1.05 }
    while(
      neglogLprof(temp_res_opt<-optim(par = temp_res_opt, fn = neglogLprof)$par)
      <=1.92-final_logL){
      res_opt[var,3]<-res_opt[var,3]*factor
    }
  }
  res_opt[-c(1:nE),]<-res_opt[-c(1:nE),,drop=F]*1000
  res_opt[nE+(1:nB),]<-abs(res_opt[nE+(1:nB),,drop=F])
  res_opt[1:nE,]<-1-exp(-abs(res_opt[1:nE,,drop=F]))
  return(c(res_opt[,1],res_opt[,2],res_opt[,3],final_logL,length(pos)))
}

# inference of transmission risk with CIs
MLinference_nonmonotonic<-function(ME,MB,pos){
  nE<-dim(ME)[2]
  nB<-dim(MB)[2]
  neglogL<-function(r){
    -(
      -sum(1000*MB[!pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))-sum(ME[!pos,]%*%abs(r[1:nE]))+
        sum(log(1-exp(-(as.vector(1000*MB[pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))+as.vector(ME[pos,]%*%abs(r[1:nE]))))))
    )
  }
  mypars<-c(rep(0.001,nE),rep(0.001,nB))
  res_opt<-optim(par = mypars, fn = neglogL)$par
  for(i in 1:25){ print(i)
    old_res<-res_opt
    mypars<-c(rep(0.001,nE),rep(0.001,nB))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  printflush(res_opt)
  for(i in 1:25){ print(i); print(Sys.time())
    old_res<-res_opt
    mypars<-c(sort(rlnorm(nE,meanlog=log(0.001),sdlog=2)),rlnorm(nB,meanlog=log(0.001),sdlog=2))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  nvars<-length(res_opt)
  printflush(res_opt)
  old_res<-res_opt
  res_opt<-optim(par = res_opt, fn = neglogL)$par
  printflush(res_opt)
  while(max(abs((res_opt-old_res)/res_opt))>0.01){
    old_res<-res_opt
    res_opt<-optim(par = res_opt, fn = neglogL)$par
    printflush(res_opt)
  }
  final_logL<-(-neglogL(res_opt))
  res_opt<-cbind(res_opt,res_opt,res_opt)
  res_opt[-c(1:nE),]<-res_opt[-c(1:nE),,drop=F]*1000
  res_opt[nE+(1:nB),]<-abs(res_opt[nE+(1:nB),,drop=F])
  res_opt[1:nE,]<-1-exp(-abs(res_opt[1:nE,,drop=F]))
  return(c(res_opt[,1],res_opt[,2],res_opt[,3],final_logL,length(pos)))
}

# inference of transmission risk with CIs
MLinference_profiled<-function(ME,MB,pos){
  nE<-dim(ME)[2]
  nB<-dim(MB)[2]
  neglogL<-function(r){
    -(
      -sum(1000*MB[!pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))-sum(ME[!pos,]%*%abs(r[1:nE]))+
        sum(log(1-exp(-(as.vector(1000*MB[pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))+as.vector(ME[pos,]%*%abs(r[1:nE]))))))
      -sum((r[1:(nE-1)]-r[2:nE])/abs(r[2:nE])*(r[2:nE]<r[1:(nE-1)]))*dim(ME)[1]*1000000
    )
  }
  res_opt<-optim(par = c(rep(0.001,nE),rep(0.001,nB)), fn = neglogL)$par
  nvars<-length(res_opt)
  for(i in 1:25){ print(i)
    old_res<-res_opt
    mypars<-c(rep(0.001,nE),rep(0.001,nB))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  printflush(res_opt)
  old_res<-res_opt
  res_opt<-optim(par = res_opt, fn = neglogL)$par
  printflush(res_opt)
  while(max(abs((res_opt-old_res)/res_opt))>0.01){
    old_res<-res_opt
    res_opt<-optim(par = res_opt, fn = neglogL)$par
    printflush(res_opt)
  }
  final_logL<-(-neglogL(res_opt))
  res_opt<-cbind(res_opt,res_opt,res_opt)
  for(var in 1:nvars){
    printflush(var)
    temp_res_opt<-res_opt
    mask<-rep(1,nE+nB);mask[var]<-0
    neglogLprof<-function(x){neglogL(x*mask+res_opt[,2]*(1-mask))}
    factor<-1+0.1
    res_opt[var,2]<-res_opt[var,2]/factor
    if(neglogLprof(optim(par = temp_res_opt, fn = neglogLprof)$par)
       >1.92-final_logL){ res_opt[var,2]<-res_opt[var,2]*factor; factor<-1.01 } else { factor<-1.05 }
    while(
      neglogLprof(temp_res_opt<-optim(par = temp_res_opt, fn = neglogLprof)$par)
      <=1.92-final_logL){
      res_opt[var,2]<-res_opt[var,2]/factor
    }
    temp_res_opt<-res_opt
    mask<-rep(1,nE+nB);mask[var]<-0
    neglogLprof<-function(x){neglogL(x*mask+res_opt[,3]*(1-mask))}
    factor<-1+0.1
    res_opt[var,2]<-res_opt[var,2]/factor
    if(neglogLprof(optim(par = temp_res_opt, fn = neglogLprof)$par)
       >1.92-final_logL){ res_opt[var,2]<-res_opt[var,2]*factor; factor<-1.01 } else { factor<-1.05 }
    while(
      neglogLprof(temp_res_opt<-optim(par = temp_res_opt, fn = neglogLprof)$par)
      <=1.92-final_logL){
      res_opt[var,3]<-res_opt[var,3]*factor
    }
  }
  res_opt[-c(1:nE),]<-res_opt[-c(1:nE),,drop=F]*1000
  res_opt[nE+(1:nB),]<-abs(res_opt[nE+(1:nB),,drop=F])
  res_opt[1:nE,]<-1-exp(-abs(res_opt[1:nE,,drop=F]))
  return(c(res_opt[,1],res_opt[,2],res_opt[,3],final_logL,length(pos)))
}

# inference of transmission risk with ascertainment
MLinference_nonmonotonic_ascertainment<-function(ME,MB,pos){
  nE<-dim(ME)[2]
  nB<-dim(MB)[2]
  neglogL<-function(r){
    temp_a<-1/(1+exp(-1000*r[nE+nB+1]))
    -(
      +sum(log(1-temp_a+temp_a*exp(-1000*MB[!pos,,drop=FALSE]%*%abs(r[nE+(1:nB)])-ME[!pos,]%*%abs(r[1:nE]))))+
        sum(log(temp_a*(1-exp(-(as.vector(1000*MB[pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))+as.vector(ME[pos,]%*%abs(r[1:nE])))))))
    )
  }
  mypars<-c(rep(0.001,nE),rep(0.001,nB),0.001)
  res_opt<-optim(par = mypars, fn = neglogL)$par
  for(i in 1:50){ print(i)
    old_res<-res_opt
    mypars<-c(rep(0.001,nE),rep(0.001,nB),0.001)
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  printflush(res_opt)
  for(i in 1:50){ print(i)
    old_res<-res_opt
    mypars<-c(sort(rlnorm(nE,meanlog=log(0.001),sdlog=2)),rlnorm(nB,meanlog=log(0.001),sdlog=2),rlnorm(1,meanlog=log(0.001),sdlog=2))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  nvars<-length(res_opt)-1
  printflush(res_opt)
  old_res<-res_opt
  res_opt<-optim(par = res_opt, fn = neglogL)$par
  printflush(res_opt)
  while(max(abs((res_opt-old_res)/res_opt))>0.01){
    old_res<-res_opt
    res_opt<-optim(par = res_opt, fn = neglogL)$par
    printflush(res_opt)
  }
  final_logL<-(-neglogL(res_opt))
  res_opt<-cbind(res_opt,res_opt,res_opt)
  res_opt[-c(1:nE),]<-res_opt[-c(1:nE),,drop=F]*1000
  res_opt[nE+(1:nB),]<-abs(res_opt[nE+(1:nB),,drop=F])
  res_opt[1:nE,]<-1-exp(-abs(res_opt[1:nE,,drop=F]))
  res_opt[nE+nB+1,]<-1/(1+exp(-(res_opt[nE+nB+1,,drop=F])))
  return(c(res_opt[,1],final_logL,length(pos)))
}

# inference of transmission risk with gamma heterogeneities
MLinference_nonmonotonic_gammahet<-function(ME,MB,pos){
  nE<-dim(ME)[2]
  nB<-dim(MB)[2]
  neglogL<-function(r){
    temp_p0<-1/(1+exp(-1000*r[nE+nB+1]))
    temp_k<-exp(1000*r[nE+nB+2])
    temp_theta<-exp(1000*r[nE+nB+3])
    avgexp<-function(x){temp_p0+(1-temp_p0)/(1-temp_theta*x)^temp_k}
    -(
      +sum(log(exp(-1000*MB[!pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))*avgexp(-ME[!pos,]%*%abs(r[1:nE]))))+
        sum(log((1-exp(-as.vector(1000*MB[pos,,drop=FALSE]%*%abs(r[nE+(1:nB)])))*avgexp(-as.vector(ME[pos,]%*%abs(r[1:nE]))))))
    )
  }
  mypars<-c(rep(0.001,nE),rep(0.001,nB),rep(0.001,3))
  res_opt<-optim(par = mypars, fn = neglogL)$par
  for(i in 1:50){ print(i)
    old_res<-res_opt
    mypars<-c(rep(0.001,nE),rep(0.001,nB),rep(0.001,3))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  printflush(res_opt)
  for(i in 1:50){ print(i)
    old_res<-res_opt
    mypars<-c(sort(rlnorm(nE,meanlog=log(0.001),sdlog=2)),rlnorm(nB,meanlog=log(0.001),sdlog=2),rlnorm(3,meanlog=log(0.001),sdlog=2))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  nvars<-length(res_opt)-1
  printflush(res_opt)
  old_res<-res_opt
  res_opt<-optim(par = res_opt, fn = neglogL)$par
  printflush(res_opt)
  while(max(abs((res_opt-old_res)/res_opt))>0.01){
    old_res<-res_opt
    res_opt<-optim(par = res_opt, fn = neglogL)$par
    printflush(res_opt)
  }
  final_logL<-(-neglogL(res_opt))
  res_opt<-cbind(res_opt,res_opt,res_opt)
  res_opt[-c(1:nE),]<-res_opt[-c(1:nE),,drop=F]*1000
  res_opt[nE+(1:nB),]<-abs(res_opt[nE+(1:nB),,drop=F])
  res_opt[1:nE,]<-1-exp(-abs(res_opt[1:nE,,drop=F]))
  res_opt[nE+nB+1,]<-1/(1+exp(-(res_opt[nE+nB+1,,drop=F])))
  res_opt[nE+nB+2,]<-exp((res_opt[nE+nB+2,,drop=F]))
  res_opt[nE+nB+3,]<-exp((res_opt[nE+nB+3,,drop=F]))
  return(c(res_opt[,1],final_logL,length(pos)))
}

# inference of transmission risk with gamma heterogeneities and ascertainment
MLinference_nonmonotonic_gammahet_ascertainment<-function(ME,MB,pos){
  nE<-dim(ME)[2]
  nB<-dim(MB)[2]
  neglogL<-function(r){
    temp_a<-1/(1+exp(-1000*r[nE+nB+1]))
    temp_p0<-1/(1+exp(-1000*r[nE+nB+2]))
    temp_k<-exp(1000*r[nE+nB+3])
    temp_theta<-exp(1000*r[nE+nB+4])
    avgexp<-function(x){temp_p0+(1-temp_p0)/(1-temp_theta*x)^temp_k}
    -(
      +sum(log(1-temp_a+temp_a*exp(-1000*MB[!pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))*avgexp(-ME[!pos,]%*%abs(r[1:nE]))))+
        sum(log(temp_a*(1-exp(-as.vector(1000*MB[pos,,drop=FALSE]%*%abs(r[nE+(1:nB)])))*avgexp(-as.vector(ME[pos,]%*%abs(r[1:nE]))))))
    )
  }
  mypars<-c(rep(0.001,nE),rep(0.001,nB),rep(0.001,4))
  res_opt<-optim(par = mypars, fn = neglogL)$par
  for(i in 1:100){ print(i)
    old_res<-res_opt
    mypars<-c(rep(0.001,nE),rep(0.001,nB),rep(0.001,4))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  printflush(res_opt)
  for(i in 1:100){ print(i)
    old_res<-res_opt
    mypars<-c(sort(rlnorm(nE,meanlog=log(0.001),sdlog=2)),rlnorm(nB,meanlog=log(0.001),sdlog=2),rlnorm(4,meanlog=log(0.001),sdlog=2))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  nvars<-length(res_opt)-1
  printflush(res_opt)
  old_res<-res_opt
  res_opt<-optim(par = res_opt, fn = neglogL)$par
  printflush(res_opt)
  while(max(abs((res_opt-old_res)/res_opt))>0.01){
    old_res<-res_opt
    res_opt<-optim(par = res_opt, fn = neglogL)$par
    printflush(res_opt)
  }
  final_logL<-(-neglogL(res_opt))
  res_opt<-cbind(res_opt,res_opt,res_opt)
  res_opt[-c(1:nE),]<-res_opt[-c(1:nE),,drop=F]*1000
  res_opt[nE+(1:nB),]<-abs(res_opt[nE+(1:nB),,drop=F])
  res_opt[1:nE,]<-1-exp(-abs(res_opt[1:nE,,drop=F]))
  res_opt[nE+nB+1,]<-1/(1+exp(-(res_opt[nE+nB+1,,drop=F])))
  res_opt[nE+nB+2,]<-1/(1+exp(-(res_opt[nE+nB+2,,drop=F])))
  res_opt[nE+nB+3,]<-exp((res_opt[nE+nB+3,,drop=F]))
  res_opt[nE+nB+4,]<-exp((res_opt[nE+nB+4,,drop=F]))
  return(c(res_opt[,1],final_logL,length(pos)))
}

# inference of transmission risk with uniform heterogeneities
MLinference_nonmonotonic_unifhet<-function(ME,MB,pos){
  nE<-dim(ME)[2]
  nB<-dim(MB)[2]
  neglogL<-function(r){
    temp_p0<-1/(1+exp(-1000*r[nE+nB+1]))
    temp_pf<-1/(1+exp(-1000*r[nE+nB+2]))
    avgexp<-function(x){temp_p0+(1-temp_p0)*temp_pf*exp(x)-(1-temp_p0)*(1-temp_pf)*(1-exp(x))/x}
    -(
      +sum(log(exp(-1000*MB[!pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))*avgexp(-ME[!pos,]%*%abs(r[1:nE]))))+
        sum(log((1-exp(-as.vector(1000*MB[pos,,drop=FALSE]%*%abs(r[nE+(1:nB)])))*avgexp(-as.vector(ME[pos,]%*%abs(r[1:nE]))))))
    )
  }
  mypars<-c(rep(0.001,nE),rep(0.001,nB),rep(0.001,2))
  res_opt<-optim(par = mypars, fn = neglogL)$par
  for(i in 1:50){ print(i)
    old_res<-res_opt
    mypars<-c(rep(0.001,nE),rep(0.001,nB),rep(0.001,2))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  printflush(res_opt)
  for(i in 1:50){ print(i)
    old_res<-res_opt
    mypars<-c(sort(rlnorm(nE,meanlog=log(0.001),sdlog=2)),rlnorm(nB,meanlog=log(0.001),sdlog=2),rlnorm(2,meanlog=log(0.001),sdlog=2))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  nvars<-length(res_opt)-1
  printflush(res_opt)
  old_res<-res_opt
  res_opt<-optim(par = res_opt, fn = neglogL)$par
  printflush(res_opt)
  while(max(abs((res_opt-old_res)/res_opt))>0.01){
    old_res<-res_opt
    res_opt<-optim(par = res_opt, fn = neglogL)$par
    printflush(res_opt)
  }
  final_logL<-(-neglogL(res_opt))
  res_opt<-cbind(res_opt,res_opt,res_opt)
  res_opt[-c(1:nE),]<-res_opt[-c(1:nE),,drop=F]*1000
  res_opt[nE+(1:nB),]<-abs(res_opt[nE+(1:nB),,drop=F])
  res_opt[1:nE,]<-1-exp(-abs(res_opt[1:nE,,drop=F]))
  res_opt[nE+nB+1,]<-1/(1+exp(-(res_opt[nE+nB+1,,drop=F])))
  res_opt[nE+nB+2,]<-1/(1+exp(-(res_opt[nE+nB+2,,drop=F])))
  return(c(res_opt[,1],final_logL,length(pos)))
}


# inference of transmission risk with uniform heterogeneities and ascertainment
MLinference_nonmonotonic_unifhet_ascertainment<-function(ME,MB,pos){
  nE<-dim(ME)[2]
  nB<-dim(MB)[2]
  neglogL<-function(r){
    temp_a<-1/(1+exp(-1000*r[nE+nB+1]))
    temp_p0<-1/(1+exp(-1000*r[nE+nB+2]))
    temp_pf<-1/(1+exp(-1000*r[nE+nB+3]))
    avgexp<-function(x){temp_p0+(1-temp_p0)*temp_pf*exp(x)-(1-temp_p0)*(1-temp_pf)*(1-exp(x))/x}
    -(
      +sum(log(1-temp_a+temp_a*exp(-1000*MB[!pos,,drop=FALSE]%*%abs(r[nE+(1:nB)]))*avgexp(-ME[!pos,]%*%abs(r[1:nE]))))+
        sum(log(temp_a*(1-exp(-as.vector(1000*MB[pos,,drop=FALSE]%*%abs(r[nE+(1:nB)])))*avgexp(-as.vector(ME[pos,]%*%abs(r[1:nE]))))))
    )
  }
  mypars<-c(rep(0.001,nE),rep(0.001,nB),rep(0.001,3))
  res_opt<-optim(par = mypars, fn = neglogL)$par
  for(i in 1:100){ print(i)
    old_res<-res_opt
    mypars<-c(rep(0.001,nE),rep(0.001,nB),rep(0.001,3))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  printflush(res_opt)
  for(i in 1:100){ print(i)
    old_res<-res_opt
    mypars<-c(sort(rlnorm(nE,meanlog=log(0.001),sdlog=2)),rlnorm(nB,meanlog=log(0.001),sdlog=2),rlnorm(3,meanlog=log(0.001),sdlog=2))
    res_opt<-optim(par = mypars, fn = neglogL)$par
    if((-neglogL(res_opt))<(-neglogL(old_res))){
      res_opt<-old_res
    }
  }
  nvars<-length(res_opt)-1
  printflush(res_opt)
  old_res<-res_opt
  res_opt<-optim(par = res_opt, fn = neglogL)$par
  printflush(res_opt)
  while(max(abs((res_opt-old_res)/res_opt))>0.01){
    old_res<-res_opt
    res_opt<-optim(par = res_opt, fn = neglogL)$par
    printflush(res_opt)
  }
  final_logL<-(-neglogL(res_opt))
  res_opt<-cbind(res_opt,res_opt,res_opt)
  res_opt[-c(1:nE),]<-res_opt[-c(1:nE),,drop=F]*1000
  res_opt[nE+(1:nB),]<-abs(res_opt[nE+(1:nB),,drop=F])
  res_opt[1:nE,]<-1-exp(-abs(res_opt[1:nE,,drop=F]))
  res_opt[nE+nB+1,]<-1/(1+exp(-(res_opt[nE+nB+1,,drop=F])))
  res_opt[nE+nB+2,]<-1/(1+exp(-(res_opt[nE+nB+2,,drop=F])))
  res_opt[nE+nB+3,]<-1/(1+exp(-(res_opt[nE+nB+3,,drop=F])))
  return(c(res_opt[,1],final_logL,length(pos)))
}


##########################################
# Maximum Likelihood of transmission risk
##########################################

# run ML inferences and extract probability of transmission per exposure

y<-table_outcome_riskScores %>%
  filter(!is.na(bg_rate_cases_app)) %>%
  filter(number_exposures>=2 & number_exposures<=6) %>%
  select(starts_with("riskScore_"),bg_rate_cases_app,positive)
ME <- data.matrix(y%>%select(starts_with("riskScore_")))
MB <- (-log(1-data.matrix(y%>%select(bg_rate_cases_app))))
pos <- pull(y,positive)==1

# non-monotonic inference
start_time<-Sys.time()
MLrisk_figX1<-MLinference_profiled_nonmonotonic(ME,MB,pos)
end_time <- Sys.time()
end_time - start_time
MLrisk_figX1 <- as.data.frame(MLrisk_figX1)
SaveTable(MLrisk_figX1)

# non-monotonic inference without CIs
start_time<-Sys.time()
MLrisk_nonmonotonic<-MLinference_nonmonotonic(ME,MB,pos)
end_time <- Sys.time()
end_time - start_time
MLrisk_nonmonotonic <- as.data.frame(MLrisk_nonmonotonic)
SaveTable(MLrisk_nonmonotonic)

ReadTable("MLrisk_figX1")
MLres<-t(data.matrix(MLrisk_figX1[1:27,1]))
MLres<-bind_rows(cbind(rep(c(bucket_labels[-1],"background"),times=3),rep(c("ML","lowerCI","upperCI"),each=9),as_tibble(t(MLres))))
colnames(MLres)<-c("riskScore","content","transmissionProbability")
MLres<-MLres %>% pivot_wider(names_from = "content", values_from = "transmissionProbability") %>% rename (transmissionProbability=ML)
table_MLrisk_figX1 <-MLres %>%
  mutate(riskScore=factor(riskScore,levels=c(bucket_labels[-1],"background")))
SaveTable(table_MLrisk_figX1)

# monotonic inference
start_time<-Sys.time()
MLrisk_figX1.1<-MLinference_profiled(ME,MB,pos)
end_time <- Sys.time()
end_time - start_time
MLrisk_figX1.1 <- as.data.frame(MLrisk_figX1.1)
SaveTable(MLrisk_figX1.1)

ReadTable("MLrisk_figX1.1")
MLres<-t(data.matrix(MLrisk_figX1.1[1:27,1]))
rho_ml<-MLres[1,9]
MLres<-bind_rows(cbind(rep(c(bucket_labels[-1],"background"),times=3),rep(c("ML","lowerCI","upperCI"),each=9),as_tibble(t(MLres))))
colnames(MLres)<-c("riskScore","content","transmissionProbability")
MLres<-MLres %>% pivot_wider(names_from = "content", values_from = "transmissionProbability") %>% rename (transmissionProbability=ML)
table_MLrisk_figX1.1 <-MLres %>%
  mutate(riskScore=factor(riskScore,levels=c(bucket_labels[-1],"background")))
SaveTable(table_MLrisk_figX1.1)

# include average risk scores in table.
mean_risk_per_bucket <- apply_by_ltla("filtered_unique_individuals_in_events",
                                      function(x){x %>% filter(getmonthhq(notification_date) %in% monthshq) %>% filter(peak_duration>=2 & peak_duration<=6) %>% select(risk_bucket,riskScore)},
                                      function(y){y %>% group_by(risk_bucket) %>% summarise(n=n(),meanRisk=mean(riskScore))},
                                      gz=save_as_gz)
table_relation_riskscore_transmissionrisk <- as_tibble(bind_cols(tibble(n=mean_risk_per_bucket %>% pull(n), x=mean_risk_per_bucket %>% pull(meanRisk)/90),read_csv("../FiguresTables/OutputTable_table_MLrisk_figX1.1.csv")[1:8,])) %>% 
  mutate(invwidth=1/((upperCI-transmissionProbability)^2+(transmissionProbability-lowerCI)^2))
SaveTable(table_relation_riskscore_transmissionrisk)

# reference transmission risk at 2 meters, 15 minutes
transmissionProbability_2m15min <- exp(predict(lm(I(log(transmissionProbability)) ~ x1+x2,
                                                  data=table_relation_riskscore_transmissionrisk %>% 
                                                    filter(x<3) %>% 
                                                    mutate(x1=x-1,x2=(x-1)^2),
                                                  weights=transmissionProbability^2*invwidth),
                                               newdata=data.frame(x1=0,x2=0)))

# risk profile by day of the week:
MLrisk_fig4a <- c()
for(weekday in 1:7){
  y<-table_outcome_riskScores %>%
    filter(!is.na(bg_rate_cases_app)) %>%
    filter(number_exposures>=2 & number_exposures<=6 & peak_day==weekday) %>%
    select(starts_with("riskScore_"),bg_rate_cases_app,positive)
  ME<-data.matrix(y%>%select(starts_with("riskScore_")))
  MB<-(-log(1-data.matrix(y%>%select(bg_rate_cases_app))))
  pos<-pull(y,positive)==1
  tempMLrisk_fig4a<-MLinference_profiled(ME,MB,pos)
  MLres<-matrix(as.data.frame(tempMLrisk_fig4a)[1:27,1],nrow=1)
  MLrisk_fig4a<-bind_rows(MLrisk_fig4a,cbind(rep(c(bucket_labels[-1],"background"),times=3),rep(c("ML","lowerCI","upperCI"),each=9),
                                             as_tibble(t(MLres)),rep(weekday,27)))
}
colnames(MLrisk_fig4a)<-c("riskScore","content","transmissionProbability","weekday")
SaveTable(MLrisk_fig4a)

MLres<-MLrisk_fig4a %>% pivot_wider(names_from = "content", values_from = "transmissionProbability") %>% rename (transmissionProbability=ML)
table_MLrisk_fig4a <- MLres %>%
  mutate(riskScore=factor(riskScore,levels=c(bucket_labels[-1],"background")))
SaveTable(table_MLrisk_fig4a)


# risk profile by day of the week, grouped:
MLrisk_fig4a2 <- c()
weekday_list<-list("Mon-Thu"=1:4,"Fri-Sun"=5:7)
for(weekday in c("Mon-Thu","Fri-Sun")){
  y<-table_outcome_riskScores %>%
    filter(!is.na(bg_rate_cases_app)) %>%
    filter(number_exposures>=2 & number_exposures<=6 & (peak_day %in% weekday_list[[weekday]])) %>%
    select(starts_with("riskScore_"),bg_rate_cases_app,positive)
  tempMLrisk_fig4a2<-MLinference_profiled(
    data.matrix(y%>%select(starts_with("riskScore_"))),
    data.matrix(y%>%select(bg_rate_cases_app)),
    pull(y,positive)==1
  )
  MLres<-matrix(as.data.frame(tempMLrisk_fig4a2)[1:27,1],nrow=1)
  MLrisk_fig4a2<-bind_rows(MLrisk_fig4a2,cbind(rep(c(bucket_labels[-1],"background"),times=3),rep(c("ML","lowerCI","upperCI"),each=9),
                                               as_tibble(t(MLres)),rep(weekday,27)))
}
colnames(MLrisk_fig4a2)<-c("riskScore","content","transmissionProbability","weekday")
SaveTable(MLrisk_fig4a2)

MLres<-MLrisk_fig4a2 %>% pivot_wider(names_from = "content", values_from = "transmissionProbability") %>% rename (transmissionProbability=ML)
table_MLrisk_fig4a2 <- MLres %>%
  mutate(riskScore=factor(riskScore,levels=c(bucket_labels[-1],"background")))
SaveTable(table_MLrisk_fig4a2)

# # risk profile by month:
# MLrisk_fig4am <- c()
# for(mymonth in monthshq){
#   y<-table_outcome_riskScores %>%
#     filter(!is.na(bg_rate_cases_app)) %>%
#     filter(number_exposures>=2 & number_exposures<=6 & month==mymonth) %>%
#     select(starts_with("riskScore_"),bg_rate_cases_app,positive)
#   tempMLrisk_fig4am<-MLinference(
#     data.matrix(y%>%select(starts_with("riskScore_"))),
#     data.matrix(y%>%select(bg_rate_cases_app)),
#     pull(y,positive)==1
#   )
#   MLres<-matrix(as.data.frame(tempMLrisk_fig4am)[1:27,1],nrow=1)
#   MLrisk_fig4am<-bind_rows(MLrisk_fig4am,cbind(rep(c(bucket_labels[-1],"background"),times=1),rep(c("ML"),each=9),
#                                                as_tibble(t(MLres)),rep(mymonth,9)))
# }
# colnames(MLrisk_fig4am)<-c("riskScore","content","transmissionProbability","month")
# SaveTable(MLrisk_fig4am)
# 
# MLres<-MLrisk_fig4am %>% select(-content)
# table_MLrisk_fig4am <- MLres %>%
#   mutate(riskScore=factor(riskScore,levels=c(bucket_labels[-1],"background")))
# SaveTable(table_MLrisk_fig4am)
# 
# # risk profile by month, grouped:
# MLrisk_fig4am2 <- c()
# month_list<-list("Jan-May"=monthshq[c(1,9:10)],"Jun-Aug"=monthshq[2:4],"Sep-Nov"=monthshq[5:7],"Dec"=monthshq[8])
# for(mymonth in c("Jan-May","Jun-Aug","Sep-Nov","Dec")){
#   y<-table_outcome_riskScores %>%
#     filter(!is.na(bg_rate_cases_app)) %>%
#     filter(number_exposures>=2 & number_exposures<=6 & (month %in% month_list[[mymonth]])) %>%
#     select(starts_with("riskScore_"),bg_rate_cases_app,positive)
#   tempMLrisk_fig4am2<-MLinference_profiled(
#     data.matrix(y%>%select(starts_with("riskScore_"))),
#     data.matrix(y%>%select(bg_rate_cases_app)),
#     pull(y,positive)==1
#   )
#   MLres<-matrix(as.data.frame(tempMLrisk_fig4am2)[1:27,1],nrow=1)
#   MLrisk_fig4am2<-bind_rows(MLrisk_fig4am2,cbind(rep(c(bucket_labels[-1],"background"),times=3),rep(c("ML","lowerCI","upperCI"),each=9),
#                                                  as_tibble(t(MLres)),rep(mymonth,27)))
# }
# colnames(MLrisk_fig4am2)<-c("riskScore","content","transmissionProbability","month")
# SaveTable(MLrisk_fig4am2)
# 
# MLres<-MLrisk_fig4am2 %>% pivot_wider(names_from = "content", values_from = "transmissionProbability") %>% rename (transmissionProbability=ML)
# table_MLrisk_fig4am2 <- MLres %>%
#   mutate(riskScore=factor(riskScore,levels=c(bucket_labels[-1],"background")))
# SaveTable(table_MLrisk_fig4am2)


# risk profile by rural/urban score:
MLrisk_fig_urban <- c()
for(location in list_locations){
  print(location); flush.console()
  y<-table_outcome_riskScores %>%
    filter(!is.na(bg_rate_cases_app)) %>%
    filter(number_exposures>=2 & number_exposures<=6 & ruc_from_ltla[ltla]==location) %>%
    select(starts_with("riskScore_"),bg_rate_cases_app,positive)
  ME<-data.matrix(y%>%select(starts_with("riskScore_")))
  MB<-(-log(1-data.matrix(y%>%select(bg_rate_cases_app))))
  pos<-pull(y,positive)==1
  tempMLrisk_fig_urban<-MLinference_profiled(ME,MB,pos) 
  MLres<-matrix(as.data.frame(tempMLrisk_fig_urban)[1:27,1],nrow=1)
  MLrisk_fig_urban<-bind_rows(MLrisk_fig_urban,cbind(rep(c(bucket_labels[-1],"background"),times=3),rep(c("ML","lowerCI","upperCI"),each=9),
                                                     as_tibble(t(MLres)),rep(location,27)))
}
colnames(MLrisk_fig_urban)<-c("riskScore","content","transmissionProbability","location")
SaveTable(MLrisk_fig_urban)

MLres<-MLrisk_fig_urban %>% pivot_wider(names_from = "content", values_from = "transmissionProbability") %>% rename (transmissionProbability=ML)
table_MLrisk_fig_urban <- MLres %>%
  mutate(riskScore=factor(riskScore,levels=c(bucket_labels[-1],"background")))
SaveTable(table_MLrisk_fig_urban)

# risk profile by rural/urban score:
MLrisk_fig_urban3 <- c()
for(location in c("rural","urban","conurbation")){
  print(location); flush.console()
  y<-table_outcome_riskScores %>%
    filter(!is.na(bg_rate_cases_app)) %>%
    filter(number_exposures>=2 & number_exposures<=6 & ruc_from_ltla_main[ltla]==location) %>%
    select(starts_with("riskScore_"),bg_rate_cases_app,positive)
  ME<-data.matrix(y%>%select(starts_with("riskScore_")))
  MB<-(-log(1-data.matrix(y%>%select(bg_rate_cases_app))))
  pos<-pull(y,positive)==1
  tempMLrisk_fig_urban3<-MLinference_profiled(ME,MB,pos) 
  MLres<-matrix(as.data.frame(tempMLrisk_fig_urban3)[1:27,1],nrow=1)
  MLrisk_fig_urban3<-bind_rows(MLrisk_fig_urban3,cbind(rep(c(bucket_labels[-1],"background"),times=3),rep(c("ML","lowerCI","upperCI"),each=9),
                                                     as_tibble(t(MLres)),rep(location,27)))
}
colnames(MLrisk_fig_urban3)<-c("riskScore","content","transmissionProbability","location")
SaveTable(MLrisk_fig_urban3)

MLres<-MLrisk_fig_urban3 %>% pivot_wider(names_from = "content", values_from = "transmissionProbability") %>% rename (transmissionProbability=ML)
table_MLrisk_fig_urban3 <- MLres %>%
  mutate(riskScore=factor(riskScore,levels=c(bucket_labels[-1],"background")))
SaveTable(table_MLrisk_fig_urban3)


# non-monotonic inference with ascertainment and different heterogeneities:
y<-table_outcome_riskScores %>%
  filter(!is.na(bg_rate_cases_app)) %>%
  filter(number_exposures>=2 & number_exposures<=6) %>%
  select(starts_with("riskScore_"),bg_rate_cases_app,positive)
ME <- data.matrix(y%>%select(starts_with("riskScore_")))
MB <- (-log(1-data.matrix(y%>%select(bg_rate_cases_app))))
pos <- pull(y,positive)==1
# ascertainment
start_time<-Sys.time()
MLrisk_asc<-MLinference_nonmonotonic_ascertainment(ME,MB,pos)
end_time <- Sys.time()
end_time - start_time
MLrisk_asc <- as.data.frame(MLrisk_asc)
SaveTable(MLrisk_asc)
# gamma heterogeneity
start_time<-Sys.time()
MLrisk_gammahet<-MLinference_nonmonotonic_gammahet(ME,MB,pos)
end_time <- Sys.time()
end_time - start_time
MLrisk_gammahet <- as.data.frame(MLrisk_gammahet)
SaveTable(MLrisk_gammahet)
# gamma heterogeneity + ascertainment
start_time<-Sys.time()
MLrisk_gammahet_asc<-MLinference_nonmonotonic_gammahet_ascertainment(ME,MB,pos)
end_time <- Sys.time()
end_time - start_time
MLrisk_gammahet_asc <- as.data.frame(MLrisk_gammahet_asc)
SaveTable(MLrisk_gammahet_asc)
# uniform heterogeneity
start_time<-Sys.time()
MLrisk_unifhet<-MLinference_nonmonotonic_unifhet(ME,MB,pos)
end_time <- Sys.time()
end_time - start_time
MLrisk_unifhet <- as.data.frame(MLrisk_unifhet)
SaveTable(MLrisk_unifhet)
# uniform heterogeneity + ascertainment
start_time<-Sys.time()
MLrisk_unifhet_asc<-MLinference_nonmonotonic_unifhet_ascertainment(ME,MB,pos)
end_time <- Sys.time()
end_time - start_time
MLrisk_unifhet_asc <- as.data.frame(MLrisk_unifhet_asc)
SaveTable(MLrisk_unifhet_asc)

ReadTable("table_relation_riskscore_transmissionrisk")
ReadTable("MLrisk_nonmonotonic")
ReadTable("MLrisk_asc")
ReadTable("MLrisk_unifhet")
ReadTable("MLrisk_unifhet_asc")
ReadTable("MLrisk_gammahet")
ReadTable("MLrisk_gammahet_asc")
table_ML_with_het_asc <- tibble(
  risk_score_bin=unlist(table_relation_riskscore_transmissionrisk[1:8,"riskScore"]), 
  average_risk_score=unlist(table_relation_riskscore_transmissionrisk[1:8,"x"]),       
  nonmonotonic_ML=(function(x){
    (1-exp(x))
  })(
    log(1-unlist(MLrisk_nonmonotonic[1:8,1]))
  ),
  with_ascertainment=(function(x){
    temp_a<-unlist(MLrisk_asc[10,1])
    temp_a*(1-exp(x))
  })(
    log(1-unlist(MLrisk_asc[1:8,1]))
  ),
  with_uniform_heterogeneities=(function(x){
    temp_p0<-unlist(MLrisk_unifhet[10,1])
    temp_pf<-unlist(MLrisk_unifhet[11,1])
    1-(temp_p0+(1-temp_p0)*temp_pf*exp(x)-(1-temp_p0)*(1-temp_pf)*(1-exp(x))/x)
  })(
    log(1-unlist(MLrisk_unifhet[1:8,1]))
  ),
  with_uniform_heterogeneities_and_ascertainment=(function(x){
    temp_a<-unlist(MLrisk_unifhet_asc[10,1])
    temp_p0<-unlist(MLrisk_unifhet_asc[11,1])
    temp_pf<-unlist(MLrisk_unifhet_asc[12,1])
    temp_a*(1-(temp_p0+(1-temp_p0)*temp_pf*exp(x)-(1-temp_p0)*(1-temp_pf)*(1-exp(x))/x))
  })(
    log(1-unlist(MLrisk_unifhet_asc[1:8,1]))
  ),
  with_gamma_heterogeneities=(function(x){
    temp_p0<-unlist(MLrisk_gammahet[10,1])
    temp_k<-unlist(MLrisk_gammahet[11,1])
    temp_theta<-unlist(MLrisk_gammahet[12,1])
    1-(temp_p0+(1-temp_p0)/(1-temp_theta*x)^temp_k)
  })(
    log(1-unlist(MLrisk_gammahet[1:8,1]))
  ),
  with_gamma_heterogeneities_and_ascertainment=(function(x){
    temp_a<-unlist(MLrisk_gammahet_asc[10,1])
    temp_p0<-unlist(MLrisk_gammahet_asc[11,1])
    temp_k<-unlist(MLrisk_gammahet_asc[12,1])
    temp_theta<-unlist(MLrisk_gammahet_asc[13,1])
    temp_a*(1-(temp_p0+(1-temp_p0)/(1-temp_theta*x)^temp_k))
  })(
    log(1-unlist(MLrisk_gammahet_asc[1:8,1]))
  )
)
SaveTable(table_ML_with_het_asc)

##########################################
# statistical analysis of exposure risk
##########################################

# reference probability of testing positive for 2 meters, 15 minutes exposure

table_maxrisk_1exposure_all_log <- apply_by_ltla("table_individual",
                                                 function(x){x %>% 
                                                     filter(number_exposures==1) %>%
                                                     filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)) %>%
                                                     select(notification_date,positive,maxRiskScore)},
                                                 function(y){y %>% 
                                                     mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmin(0.75,cut_interval(log(maxRiskScore),n=20,labels=FALSE)/20)*diff(range(log(maxRiskScore)))+min(log(maxRiskScore)))) %>%
                                                     group_by(score) %>%
                                                     summarise(TPAEN=mean(positive), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                                 },
                                                 gz=save_as_gz)
SaveTable(table_maxrisk_1exposure_all_log)

ReadTable("table_maxrisk_1exposure_all_log")
ReadTable("table_relation_riskscore_transmissionrisk")
TPAEN_2m15min <- table_maxrisk_1exposure_all_log %>%
  mutate(x=score/90) %>% #OLD VERSION: mutate(riskScore=cut(score,breaks = bucket_boundaries,labels = bucket_labels)) %>% left_join(table_relation_riskscore_transmissionrisk,.,by=c("riskScore")) %>%
  mutate(invwidth=1/((upperCI-TPAEN)^2+(TPAEN-lowerCI)^2)) %>%
  filter(x<3) %>% 
  mutate(x1=x-1,x2=(x-1)^2) %>% 
  lm(I(log(TPAEN)) ~ x1 + x2,data =.,weights=TPAEN^2*invwidth) %>% 
  predict(.,newdata=data.frame(x1=0,x2=0)) %>% 
  exp(.)

# run basic inferences and extract typical probability of testing positive

# max risk score
table_maxrisk_all_log <- apply_by_ltla("table_individual",
                                       function(x){x %>% 
                                           filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)) %>%
                                           select(notification_date,positive,maxRiskScore)},
                                       function(y){y %>% 
                                           mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmin(0.75,cut_interval(log(maxRiskScore),n=20,labels=FALSE)/20)*diff(range(log(maxRiskScore)))+min(log(maxRiskScore)))) %>%
                                           group_by(score) %>%
                                           summarise(n=n(),TPAEN=mean(positive),meanscore=mean(maxRiskScore),upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                       },gz=save_as_gz)
SaveTable(table_maxrisk_all_log)
# correcting for ML background
table_maxrisk_all_log_nobg <- apply_by_ltla("table_individual",
                                        function(x){x %>% 
                                            filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)-(1-(1-bg_rate_cases_app)^rho_ml)) %>%
                                            select(notification_date,positive,maxRiskScore)},
                                        function(y){y %>% 
                                            mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmin(0.75,cut_interval(log(maxRiskScore),n=20,labels=FALSE)/20)*diff(range(log(maxRiskScore)))+min(log(maxRiskScore)))) %>%
                                            group_by(score) %>%
                                            summarise(n=n(),TPAEN=mean(positive),meanscore=mean(maxRiskScore), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                        },gz=save_as_gz)
SaveTable(table_maxrisk_all_log_nobg)
# max risk score by month
table_maxrisk_log <- apply_by_ltla("table_individual",
                                   function(x){x %>% 
                                       filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)) %>%
                                       select(notification_date,positive,maxRiskScore)},
                                   function(y){y %>% 
                                       mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmin(0.75,cut_interval(log(maxRiskScore),n=20,labels=FALSE)/20)*diff(range(log(maxRiskScore)))+min(log(maxRiskScore)))) %>%
                                       group_by(score,month) %>%
                                       summarise(n=n(),TPAEN=mean(positive),meanscore=mean(maxRiskScore),upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                   },gz=save_as_gz)
SaveTable(table_maxrisk_log)
# correcting for ML background
table_maxrisk_log_nobg <- apply_by_ltla("table_individual",
                                        function(x){x %>% 
                                            filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)-(1-(1-bg_rate_cases_app)^rho_ml)) %>%
                                            select(notification_date,positive,maxRiskScore)},
                                        function(y){y %>% 
                                            mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmin(0.75,cut_interval(log(maxRiskScore),n=20,labels=FALSE)/20)*diff(range(log(maxRiskScore)))+min(log(maxRiskScore)))) %>%
                                            group_by(score,month) %>%
                                            summarise(n=n(),TPAEN=mean(positive),meanscore=mean(maxRiskScore), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                        },gz=save_as_gz)
SaveTable(table_maxrisk_log_nobg)
# correcting for full background
table_maxrisk_log_nofullbg <- apply_by_ltla("table_individual",
                                            function(x){x %>% 
                                                filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)-bg_rate_cases_app) %>%
                                                select(notification_date,positive,maxRiskScore)},
                                            function(y){y %>% 
                                                mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmin(0.75,cut_interval(log(maxRiskScore),n=20,labels=FALSE)/20)*diff(range(log(maxRiskScore)))+min(log(maxRiskScore)))) %>%
                                                group_by(score,month) %>%
                                                summarise(n=n(),TPAEN=mean(positive),meanscore=mean(maxRiskScore), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                            },gz=save_as_gz)
SaveTable(table_maxrisk_log_nofullbg)

# cumulative risk
table_cumrisk_all_log <- apply_by_ltla("table_individual",
                                       function(x){x %>% 
                                           filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)) %>%
                                           select(notification_date,positive,number_exposures,meanRiskScore)},
                                       function(y){y %>% 
                                           mutate(score=exp(pmin(0.75,cut_interval(log(meanRiskScore*number_exposures),n=20,labels=FALSE)/20)*diff(range(log(meanRiskScore*number_exposures)))+min(log(meanRiskScore*number_exposures)))) %>%
                                           group_by(score) %>%
                                           summarise(n=n(),TPAEN=mean(positive),meanscore=mean(meanRiskScore*number_exposures), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                       },gz=save_as_gz)
SaveTable(table_cumrisk_all_log)
# correcting for ML background
table_cumrisk_all_log_nobg <- apply_by_ltla("table_individual",
                                        function(x){x %>% 
                                            filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)-(1-(1-bg_rate_cases_app)^rho_ml)) %>%
                                            select(notification_date,positive,number_exposures,meanRiskScore)},
                                        function(y){y %>% 
                                            mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmin(0.75,cut_interval(log(meanRiskScore*number_exposures),n=20,labels=FALSE)/20)*diff(range(log(meanRiskScore*number_exposures)))+min(log(meanRiskScore*number_exposures)))) %>%
                                            group_by(score) %>%
                                            summarise(n=n(),TPAEN=mean(positive),meanscore=mean(meanRiskScore*number_exposures), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                        },gz=save_as_gz)
SaveTable(table_cumrisk_all_log_nobg)
# cumulative risk per month
table_cumrisk_log <- apply_by_ltla("table_individual",
                                   function(x){x %>% 
                                       filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)) %>%
                                       select(notification_date,positive,number_exposures,meanRiskScore)},
                                   function(y){y %>% 
                                       mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmin(0.75,cut_interval(log(meanRiskScore*number_exposures),n=20,labels=FALSE)/20)*diff(range(log(meanRiskScore*number_exposures)))+min(log(meanRiskScore*number_exposures)))) %>%
                                       group_by(score,month) %>%
                                       summarise(n=n(),TPAEN=mean(positive),meanscore=mean(meanRiskScore*number_exposures), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                   },gz=save_as_gz)
SaveTable(table_cumrisk_log)
# correcting for ML background
table_cumrisk_log_nobg <- apply_by_ltla("table_individual",
                                        function(x){x %>% 
                                            filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)-(1-(1-bg_rate_cases_app)^rho_ml)) %>%
                                            select(notification_date,positive,number_exposures,meanRiskScore)},
                                        function(y){y %>% 
                                            mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmin(0.75,cut_interval(log(meanRiskScore*number_exposures),n=20,labels=FALSE)/20)*diff(range(log(meanRiskScore*number_exposures)))+min(log(meanRiskScore*number_exposures)))) %>%
                                            group_by(score,month) %>%
                                            summarise(n=n(),TPAEN=mean(positive),meanscore=mean(meanRiskScore*number_exposures), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                        },gz=save_as_gz)
SaveTable(table_cumrisk_log_nobg)
# correcting for full background
table_cumrisk_log_nofullbg <- apply_by_ltla("table_individual",
                                            function(x){x %>% 
                                                filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)-bg_rate_cases_app) %>%
                                                select(notification_date,positive,number_exposures,meanRiskScore)},
                                            function(y){y %>% 
                                                mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmin(0.75,cut_interval(log(meanRiskScore*number_exposures),n=20,labels=FALSE)/20)*diff(range(log(meanRiskScore*number_exposures)))+min(log(meanRiskScore*number_exposures)))) %>%
                                                group_by(score,month) %>%
                                                summarise(n=n(),TPAEN=mean(positive),meanscore=mean(meanRiskScore*number_exposures), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                            },gz=save_as_gz)
SaveTable(table_cumrisk_log_nofullbg)


# duration
table_plot_TPAEN_duration_all_log<- apply_by_ltla("table_individual",
                                                  function(x){x %>%
                                                      filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)) %>%
                                                      select(notification_date,positive,meanDuration,number_exposures)},
                                                  function(y){y %>%
                                                      mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmax(0.25,pmin(0.9,cut_interval(log(number_exposures*meanDuration/3600),n=20,labels=FALSE)/20))*diff(range(log(number_exposures*meanDuration/3600)))+min(log(number_exposures*meanDuration/3600)))) %>%
                                                      group_by(score) %>%
                                                      summarise(n=n(),TPAEN=mean(positive),meanscore=mean(number_exposures*meanDuration/3600), count=n(), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                                  },gz=save_as_gz)
SaveTable(table_plot_TPAEN_duration_all_log)
# correcting for ML background
table_plot_TPAEN_duration_all_log_nobg<- apply_by_ltla("table_individual",
                                                   function(x){x %>%
                                                       filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)-(1-(1-bg_rate_cases_app)^rho_ml)) %>%
                                                       select(notification_date,positive,meanDuration,number_exposures)},
                                                   function(y){y %>%
                                                       mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmax(0.25,pmin(0.8,cut_interval(log(number_exposures*meanDuration/3600),n=20,labels=FALSE)/20))*diff(range(log(number_exposures*meanDuration/3600)))+min(log(number_exposures*meanDuration/3600)))) %>%
                                                       group_by(score) %>%
                                                       summarise(n=n(),TPAEN=mean(positive),meanscore=mean(number_exposures*meanDuration/3600), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                                   },gz=save_as_gz)
SaveTable(table_plot_TPAEN_duration_all_log_nobg)
# duration per month
table_plot_TPAEN_duration_log<- apply_by_ltla("table_individual",
                                              function(x){x %>%
                                                  filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)) %>%
                                                  select(notification_date,positive,meanDuration,number_exposures)},
                                              function(y){y %>%
                                                  mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmax(0.25,pmin(0.9,cut_interval(log(number_exposures*meanDuration/3600),n=20,labels=FALSE)/20))*diff(range(log(number_exposures*meanDuration/3600)))+min(log(number_exposures*meanDuration/3600)))) %>%
                                                  group_by(score,month) %>%
                                                  summarise(n=n(),TPAEN=mean(positive),meanscore=mean(number_exposures*meanDuration/3600), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                              },gz=save_as_gz)
SaveTable(table_plot_TPAEN_duration_log)
# correcting for ML background
table_plot_TPAEN_duration_log_nobg<- apply_by_ltla("table_individual",
                                                   function(x){x %>%
                                                       filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)-(1-(1-bg_rate_cases_app)^rho_ml)) %>%
                                                       select(notification_date,positive,meanDuration,number_exposures)},
                                                   function(y){y %>%
                                                       mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmax(0.25,pmin(0.8,cut_interval(log(number_exposures*meanDuration/3600),n=20,labels=FALSE)/20))*diff(range(log(number_exposures*meanDuration/3600)))+min(log(number_exposures*meanDuration/3600)))) %>%
                                                       group_by(score,month) %>%
                                                       summarise(n=n(),TPAEN=mean(positive),meanscore=mean(number_exposures*meanDuration/3600), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                                   },gz=save_as_gz)
SaveTable(table_plot_TPAEN_duration_log_nobg)
# correcting for full background
table_plot_TPAEN_duration_log_nofullbg<- apply_by_ltla("table_individual",
                                                       function(x){x %>%
                                                           filter(getmonthhq(notification_date) %in% monthshq) %>% mutate(positive=tests_positive(positive,number_exposures)-bg_rate_cases_app) %>%
                                                           select(notification_date,positive,meanDuration,number_exposures)},
                                                       function(y){y %>%
                                                           mutate(month=factor(getmonthhq(notification_date), levels = monthshq, ordered = T),score=exp(pmax(0.25,pmin(0.8,cut_interval(log(number_exposures*meanDuration/3600),n=20,labels=FALSE)/20))*diff(range(log(number_exposures*meanDuration/3600)))+min(log(number_exposures*meanDuration/3600)))) %>%
                                                           group_by(score,month) %>%
                                                           summarise(n=n(),TPAEN=mean(positive),meanscore=mean(number_exposures*meanDuration/3600), upperCI=qbinom(p = 0.975, size = n(), prob = pmax(1,sum(positive))/n())/n(), lowerCI=qbinom(p = 0.025, size = n(), prob = pmax(1,sum(positive))/n())/n())
                                                       },gz=save_as_gz)
SaveTable(table_plot_TPAEN_duration_log_nofullbg)


# joint tables
ReadTable("table_maxrisk_log")
ReadTable("table_maxrisk_log_nobg")
ReadTable("table_maxrisk_log_nofullbg")
ReadTable("table_cumrisk_log")
ReadTable("table_cumrisk_log_nobg")
ReadTable("table_cumrisk_log_nofullbg")
ReadTable("table_plot_TPAEN_duration_log")
ReadTable("table_plot_TPAEN_duration_log_nobg")
ReadTable("table_plot_TPAEN_duration_log_nofullbg")
ReadTable("table_maxrisk_all_log")
ReadTable("table_cumrisk_all_log")
ReadTable("table_plot_TPAEN_duration_all_log")
ReadTable("table_maxrisk_all_log_nobg")
ReadTable("table_cumrisk_all_log_nobg")
ReadTable("table_plot_TPAEN_duration_all_log_nobg")
table_maxrisk_joint <- bind_rows(
  table_maxrisk_log %>% mutate(name="with background risk"),
  table_maxrisk_log_nobg %>% mutate(name="ML background risk correction"),
  table_maxrisk_log_nofullbg %>% mutate(name="full background risk correction")) %>%
  mutate(name=factor(name,ordered = TRUE,levels = c("with background risk","ML background risk correction","full background risk correction")))
SaveTable(table_maxrisk_joint)
table_cumrisk_joint <- bind_rows(
  table_cumrisk_log %>% mutate(name="with background risk"),
  table_cumrisk_log_nobg %>% mutate(name="ML background risk correction"),
  table_cumrisk_log_nofullbg %>% mutate(name="full background risk correction")) %>%
  mutate(name=factor(name,ordered = TRUE,levels = c("with background risk","ML background risk correction","full background risk correction")))
SaveTable(table_cumrisk_joint)
table_duration_joint <- bind_rows(
  table_plot_TPAEN_duration_log %>% mutate(name="with background risk"),
  table_plot_TPAEN_duration_log_nobg %>% mutate(name="ML background risk correction"),
  table_plot_TPAEN_duration_log_nofullbg %>% mutate(name="full background risk correction")) %>%
  mutate(name=factor(name,ordered = TRUE,levels = c("with background risk","ML background risk correction","full background risk correction")))
SaveTable(table_duration_joint)
table_all_joint <- bind_rows(
  table_maxrisk_log %>% mutate(name="maximum risk score"),
  table_cumrisk_log %>% mutate(name="cumulative risk score"),
  table_plot_TPAEN_duration_log %>% mutate(name="duration (hours)")) %>%
  mutate(name=factor(name,ordered = TRUE,levels = c("maximum risk score","cumulative risk score","duration (hours)")))
SaveTable(table_all_joint)
table_overall_joint <- bind_rows(
  table_maxrisk_all_log %>% mutate(name="maximum risk score"),
  table_cumrisk_all_log %>% mutate(name="cumulative risk score"),
  table_plot_TPAEN_duration_all_log %>% mutate(name="duration (hours)")) %>%
  mutate(name=factor(name,ordered = TRUE,levels = c("maximum risk score","cumulative risk score","duration (hours)")))
SaveTable(table_overall_joint)
table_overall_joint_nobg <- bind_rows(
  table_maxrisk_all_log_nobg %>% mutate(name="maximum risk score"),
  table_cumrisk_all_log_nobg %>% mutate(name="cumulative risk score"),
  table_plot_TPAEN_duration_all_log_nobg %>% mutate(name="duration (hours)")) %>%
  mutate(name=factor(name,ordered = TRUE,levels = c("maximum risk score","cumulative risk score","duration (hours)")))
SaveTable(table_overall_joint_nobg)

table_summary_distributions<-bind_rows(table_outcome_riskScores %>% 
                                         filter(getmonthhq(date) %in% monthshq) %>%
                                         summarise(across="all contacts",
                                                   mean=mean(duration/3600),
                                                   median=quantile(duration/3600,1/2),
                                                   q25=quantile(duration/3600,1/4),
                                                   q75=quantile(duration/3600,3/4),
                                                   m1h=sum(duration>3600)/sum(duration>0),
                                                   m4h=sum(duration>4*3600)/sum(duration>0)
                                         ),
                                       table_outcome_riskScores %>% 
                                         filter(getmonthhq(date) %in% monthshq) %>%
                                         arrange(duration) %>% 
                                         mutate(transmissions=positive-(1-(1-bg_rate_cases_app)^rho_ml),sum_transmissions=cumsum(transmissions),tot_transmissions=sum(transmissions)) %>% 
                                         summarise(across="transmissions",
                                                   mean=sum(transmissions*duration/3600)/sum(transmissions),
                                                   q25=duration[min(which(sum_transmissions>0.25*tot_transmissions))]/3600,
                                                   median=duration[min(which(sum_transmissions>0.5*tot_transmissions))]/3600,
                                                   q75=duration[min(which(sum_transmissions>0.75*tot_transmissions))]/3600,
                                                   m1h=sum(transmissions[duration>3600])/sum(transmissions),
                                                   m4h=sum(transmissions[duration>4*3600])/sum(transmissions)
                                         )
)
SaveTable(table_summary_distributions)

# classification of exposures

table_classification <- table_outcome_riskScores %>% 
  filter(getmonthhq(date) %in% monthshq) %>%
  mutate(setting=factor(
    if_else(classification=='household',classification,if_else(number_exposures==1,"fleeting",if_else(peak_duration<number_exposures,"recurring","single day"))),
    levels=c('household','recurring','single day','fleeting'))) %>%
  group_by(setting) %>% 
  summarise(contacts=n(),
            cumulative_duration=sum(duration/3600),
            cumulative_risk=sum(cumRiskScore/90),
            #background_infections=sum((1-(1-bg_rate_cases_app)^rho_ml)),
            transmissions=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml)),
            duration=mean(duration/3600),
            risk=mean(cumRiskScore/90),
            SAR=transmissions/contacts,
            mean_risk_per_hour=cumulative_risk/cumulative_duration) %>%
  mutate(across(c(contacts,cumulative_duration,cumulative_risk,transmissions),function(x){x/sum(x)},.names="fraction_of_{.col}"))
SaveTable(table_classification)

