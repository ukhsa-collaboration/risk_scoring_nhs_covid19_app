# generate plots
using("ggpubr")
using("cowplot")
using("ggpmisc")


setwd(dir_ltla)

table_outcome_riskScores <- read_csv(file=paste0("table_riskScores.",extension_csv)) %>% 
  mutate(month=getmonthhq(date)) %>%
  filter(month %in% monthshq) 

# load info basic risk + background risk
ReadTable("MLrisk_figX1.1")
MLres<-t(data.matrix(MLrisk_figX1.1[1:27,1]))
rho_ml<-MLres[1,9]
rho_ml_lower<-MLres[1,18]
rho_ml_upper<-MLres[1,27]

# reload rural/urban scores if needed
ruc <- read.table(ruc_file,sep="\t",header=T,quote="|")
rownames(ruc) <- ruc$LAD11CD
ruc_from_ltla <- as.character(ruc[ltlas,"RUC11CD"])
ruc_from_ltla_main <- as.character(ruc[ltlas,"RUC11CD"])
ruc_from_ltla_main[ruc_from_ltla_main %in% c('1','2')]<-"rural"
ruc_from_ltla_main[ruc_from_ltla_main %in% c('3','4')]<-"urban"
ruc_from_ltla_main[ruc_from_ltla_main %in% c('5','6')]<-"conurbation"
ruc_from_ltla_simple <- as.numeric(ruc_from_ltla)
ruc_from_ltla_simple[region_from_ltla[ltlas]=="London"]<-0
ruc_from_ltla_simple[region_from_ltla[ltlas]=="Wales"]<-7
ruc_from_ltla_simple[is.na(ruc_from_ltla_simple)]<-3.5
ruc_from_ltla[region_from_ltla[ltlas]=="Wales"]<-"Wales"
ruc_from_ltla[region_from_ltla[ltlas]=="London"]<-"London"
ruc_from_ltla[ruc_from_ltla %in% c('5','6')]<-"5/6"
names(ruc_from_ltla)<-ltlas
names(ruc_from_ltla_simple)<-ltlas
names(ruc_from_ltla_main)<-ltlas
list_locations <- sort(unique(ruc_from_ltla[ruc_from_ltla!='NA']))

# compute transmission risk for risk=1 (weighted quadratic extrapolation of log-scale probability for risk<3)
ReadTable("table_relation_riskscore_transmissionrisk")
transmissionProbability_2m15min <- exp(predict(lm(I(log(transmissionProbability)) ~ x1+x2,
                                                  data=table_relation_riskscore_transmissionrisk %>% 
                                                    filter(x<3) %>% 
                                                    mutate(x1=x-1,x2=(x-1)^2),
                                                  weights=transmissionProbability^2*invwidth),
                                               newdata=data.frame(x1=0,x2=0)))
# compute infection risk for risk=1
ReadTable("table_maxrisk_1exposure_all_log")
TPAEN_2m15min <- table_maxrisk_1exposure_all_log %>%
  mutate(x=score/90) %>% #OLD VERSION: mutate(riskScore=cut(score,breaks = bucket_boundaries,labels = bucket_labels)) %>% left_join(table_relation_riskscore_transmissionrisk,.,by=c("riskScore")) %>%
  mutate(invwidth=1/((upperCI-TPAEN)^2+(TPAEN-lowerCI)^2)) %>%
  filter(x<3) %>% 
  mutate(x1=x-1,x2=(x-1)^2) %>% 
  lm(I(log(TPAEN)) ~ x1 + x2,data =.,weights=TPAEN^2*invwidth) %>% 
  predict(.,newdata=data.frame(x1=0,x2=0)) %>% 
  exp(.)

# figure showing contact tracing guidelines in distance/duration space
d1<-data.frame(x=seq(0,5,0.01),y=15*10/9*pmax(1,seq(0,5,0.01))^2/4,type="Risky contact for\n  NHS COVID-19 app,\n  standard infectiousness")
d2<-data.frame(x=seq(0,5,0.01),y=15*10/9*0.4*pmax(1,seq(0,5,0.01))^2/4,type="Risky contact for\n  NHS COVID-19 app,\n  high infectiousness")
d3<-data.frame(x=seq(0,5,0.01),y=15+(seq(0,5,0.01)>2)*(seq(0,5,0.01)-2)*1000,type="Widely used guidelines\n  to trace close contacts\n  (>15 min at <2m)")
d<-bind_rows(d1,d2,d3)
plot_Figure_X0 <- ggplot(data=d) + xlab("distance (meters)") + ylab("duration (minutes)") +
  geom_ribbon(aes(x=x,ymin=y,ymax=Inf,fill=type),alpha=0.3) + 
  theme_classic() + 
  scale_x_continuous(expand = expansion(mult = 0),limits = c(0,4.5)) + 
  scale_y_continuous(expand = expansion(mult = 0),limits = c(0,30),n.breaks = 7) +
  labs(fill='Notified contact:') 
SaveFigure(plot_Figure_X0,width=15,height=10)

# main figure of relation between risk score and transmission risk
ReadTable("table_relation_riskscore_transmissionrisk")
using("MASS")
select<-dplyr::select
plot_relation_riskscore_transmissionrisk<-ggplot(
  data=table_relation_riskscore_transmissionrisk,
  aes(x=x,y=transmissionProbability,group=1)) + 
  geom_smooth(method='rlm', formula= y~x+0, mapping = aes(weight = invwidth),fullrange=TRUE,fill="lightblue") + 
  geom_errorbar(aes(ymin=lowerCI,ymax=upperCI),width=0.1) +
  geom_point() +
  theme_classic() + 
  geom_point(size=5,col="gray",data = tibble(x=1,transmissionProbability=transmissionProbability_2m15min)) +
  geom_text(size=3,hjust = 1, nudge_x = 0.06,nudge_y = 0.0009,
            aes(x=x,y=transmissionProbability,label="2m,15min"),col="darkgray",
            data = tibble(x=1,transmissionProbability=transmissionProbability_2m15min)) +
  scale_y_continuous(limits = c(0,NA)) + 
  scale_x_continuous(limits = c(0,NA)) +
  coord_cartesian( xlim = c(0,NA),ylim = c(0,NA), expand = FALSE) +
  xlab("app-measured risk score") + 
  ylab("probability of reported transmission\nper exposure window")
SaveFigure(plot_relation_riskscore_transmissionrisk,width=15,height=12)
#logscale
plot_relation_riskscore_transmissionrisk_log<-ggplot(
  data=table_relation_riskscore_transmissionrisk,
  aes(x=x,y=transmissionProbability,group=1)) + 
  geom_smooth(method='rlm', formula= y~x, mapping = aes(weight = invwidth),fullrange=TRUE,fill="lightblue") + 
  geom_errorbar(aes(ymin=lowerCI,ymax=upperCI),width=0.1) +
  geom_point() +
  theme_classic() + 
  geom_point(size=5,col="gray",data = tibble(x=1,transmissionProbability=transmissionProbability_2m15min)) +
  geom_text(size=3,hjust = 1, nudge_x = 0.06,nudge_y = 0.0009,
            aes(x=x,y=transmissionProbability,label="2m,15min"),col="darkgray",
            data = tibble(x=1,transmissionProbability=transmissionProbability_2m15min)) +
  scale_y_continuous(limits = c(0.0025,NA),trans="log10") + 
  scale_x_continuous(limits = c(0.75,NA),trans="log10") +
  coord_cartesian( expand = FALSE) +
  xlab("app-measured risk score") + 
  ylab("probability of reported transmission\nper exposure window")
SaveFigure(plot_relation_riskscore_transmissionrisk_log,width=15,height=12)

# supplementary figure of relation between risk score and transmission risk
# with nonmonotonicity, heterogeneities, ascertainment

ReadTable("table_ML_with_het_asc")
ReadTable("table_MLrisk_figX1.1")
plot_MLrisk_figX1b <- ggplot(data=table_MLrisk_figX1.1 %>% filter(riskScore!="background") %>% left_join(.,table_ML_with_het_asc %>% select(risk_score_bin,average_risk_score),by=c(riskScore="risk_score_bin")) %>% mutate(model="monotonic_ML"),aes(group=model)) +
  geom_line(aes(x=average_risk_score,y=transmissionProbability),col="darkgray",lwd=2) +
  geom_ribbon(aes(x=average_risk_score,ymin=lowerCI,ymax=upperCI),alpha=0.3,fill="lightgray") +
  geom_line(data = table_ML_with_het_asc %>% pivot_longer(cols = c(nonmonotonic_ML,starts_with("with_")),names_to="model", values_to="transmissionProbability"),
            aes(x=average_risk_score,y=transmissionProbability,col=model)) +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  #scale_y_log10() +
  coord_cartesian(ylim=c(0,NA),xlim=c(0,NA),expand = FALSE) + 
  #ylab("risk of detectable transmission\n from a single exposure window") +
  ylab("probability of reported transmission\nper exposure window") +
  xlab("risk score") +
  scale_color_discrete(labels = ~ gsub("_", " ", .x))
SaveFigure(plot_MLrisk_figX1b,width=20,height=12)

plot_MLrisk_figX1b_log <- ggplot(data=table_MLrisk_figX1.1 %>% filter(riskScore!="background") %>% left_join(.,table_ML_with_het_asc %>% select(risk_score_bin,average_risk_score),by=c(riskScore="risk_score_bin")) %>% mutate(model="monotonic_ML"),aes(group=model)) +
  geom_line(aes(x=average_risk_score,y=transmissionProbability),col="darkgray",lwd=2) +
  geom_ribbon(aes(x=average_risk_score,ymin=lowerCI,ymax=upperCI),alpha=0.3,fill="lightgray") +
  geom_line(data = table_ML_with_het_asc %>% pivot_longer(cols = c(nonmonotonic_ML,starts_with("with_")),names_to="model", values_to="transmissionProbability"),
            aes(x=average_risk_score,y=transmissionProbability,col=model)) +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_y_log10() +
  scale_x_log10() +
  ylab("probability of reported transmission\nper exposure window") +
  xlab("risk score") +
  scale_color_discrete(labels = ~ gsub("_", " ", .x))
SaveFigure(plot_MLrisk_figX1b_log,width=20,height=12)

# individual plots on TPAEN vs maximum risk score
ReadTable("table_maxrisk_all_log")
plot_maxrisk_X1aa <- ggplot(data=table_maxrisk_all_log,
                            aes(x=meanscore/90,y=TPAEN)) + 
  geom_line(lwd=1) + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lowerCI,ymax=upperCI),width=0.02)+
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI),alpha=0.5,fill="lightgray") + 
  xlab("maximum risk score") +
  ylab("probability of\nreported infection") +
  theme_classic() +
  geom_point(size=5,col="gray",aes(x=x,y=TPAEN),data = tibble(x=1,TPAEN=TPAEN_2m15min)) +
  geom_text(size=4,hjust = 0, nudge_x = -0.08,nudge_y = 0.007,
            aes(x=x,y=TPAEN,label="2m,15min"),col="darkgray",
            data = tibble(x=1,TPAEN=TPAEN_2m15min)) +
  scale_y_continuous(limits=c(0,NA)) +
  scale_x_continuous(limits=c(0.8,NA),breaks = c(1,2,3,4,5,10),trans = "log10") +
  coord_cartesian(expand=FALSE)
SaveFigure(plot_maxrisk_X1aa,width = 15,height = 12)
ReadTable("table_maxrisk_log")
plot_maxrisk_X1a <- ggplot(data=table_maxrisk_log %>% mutate(month=factor(month,levels=monthshq)),
                           aes(x=meanscore/90,y=TPAEN,group=month)) + 
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) + 
  xlab("maximum risk score") +
  ylab("probability of reported infection") +
  theme_classic() + 
  scale_y_continuous(limits=c(0,NA)) +
  coord_cartesian(expand = F) +
  scale_x_continuous(breaks = c(1,2,3,4,5,10),trans = "log10") 
SaveFigure(plot_maxrisk_X1a,width = 15,height = 12)
ReadTable("table_maxrisk_log_nobg")
plot_maxrisk_log_nobg <- ggplot(data=table_maxrisk_log_nobg %>% mutate(month=factor(month,levels=monthshq)),
                                aes(x=meanscore/90,y=TPAEN,group=month)) + 
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) + 
  xlab("maximum risk score") +
  ylab("probability of reported infection") +
  theme_classic() + 
  coord_cartesian(expand = F) +
  scale_x_continuous(breaks = c(1,2,3,4,5,10),trans = "log10") 
SaveFigure(plot_maxrisk_log_nobg,width = 15,height = 12)
ReadTable("table_maxrisk_log_nofullbg")
plot_maxrisk_log_nofullbg <- ggplot(data=table_maxrisk_log_nofullbg %>% mutate(month=factor(month,levels=monthshq)),
                                    aes(x=meanscore/90,y=TPAEN,group=month)) + 
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) + 
  xlab("maximum risk score") +
  ylab("probability of reported infection") +
  theme_classic() + 
  coord_cartesian(expand = F) +
  scale_x_continuous(breaks = c(1,2,3,4,5,10),trans = "log10") 
SaveFigure(plot_maxrisk_log_nofullbg,width = 15,height = 12)



# individual plots on TPAEN vs cumulative risk score
ReadTable("table_cumrisk_all_log")
plot_cumrisk_X4ba <- ggplot(data=table_cumrisk_all_log,
                            aes(x=meanscore/90,y=TPAEN,group=1)) + 
  geom_line(lwd=1) + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lowerCI,ymax=upperCI),width=0.02)+
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI),alpha=0.5,fill="lightgray") + 
  xlab("cumulative risk score") +
  ylab("probability of reported infection") +
  theme_classic() + 
  geom_point(size=5,col="gray",aes(x=x,y=TPAEN),data = tibble(x=1,TPAEN=TPAEN_2m15min)) +
  geom_text(size=4,hjust = 0, nudge_x = -0.15,nudge_y = 0.01,
            aes(x=x,y=TPAEN,label="2m,15min"),col="darkgray",
            data = tibble(x=1,TPAEN=TPAEN_2m15min)) +
  scale_y_continuous(limits=c(0,NA)) +
  scale_x_log10() +
  coord_cartesian(expand=FALSE)
SaveFigure(plot_cumrisk_X4ba,width = 15,height = 12)
ReadTable("table_cumrisk_log")
plot_cumrisk_X4b <- ggplot(data=table_cumrisk_log %>% mutate(month=factor(month,levels=monthshq)),
                           aes(x=meanscore/90,y=TPAEN,group=month)) + 
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) + 
  xlab("cumulative risk score") +
  ylab("probability of reported infection") +
  theme_classic() + 
  scale_y_continuous(trans="log10") +
  scale_x_log10()
SaveFigure(plot_cumrisk_X4b,width = 15,height = 12)
ReadTable("table_cumrisk_log_nobg")
plot_cumrisk_log_nobg <- ggplot(data=table_cumrisk_log_nobg %>% mutate(month=factor(month,levels=monthshq)),
                                aes(x=meanscore/90,y=TPAEN,group=month)) + 
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) + 
  xlab("cumulative risk score") +
  ylab("probability of reported infection") +
  theme_classic() + 
  scale_x_log10() +
  coord_cartesian(expand=F)
SaveFigure(plot_cumrisk_log_nobg,width = 15,height = 12)
ReadTable("table_cumrisk_log_nofullbg")
plot_cumrisk_log_nofullbg <- ggplot(data=table_cumrisk_log_nofullbg %>% mutate(month=factor(month,levels=monthshq)),
                                    aes(x=meanscore/90,y=TPAEN,group=month)) + 
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) + 
  xlab("cumulative risk score") +
  ylab("probability of reported infection") +
  theme_classic() + 
  scale_x_log10() +
  coord_cartesian(expand=F)
SaveFigure(plot_cumrisk_log_nofullbg,width = 15,height = 12)



# individual plots on TPAEN vs duration
ReadTable("table_plot_TPAEN_duration_all_log")
plot_duration_all <- ggplot(data=table_plot_TPAEN_duration_all_log,
                            aes(x=meanscore,y=TPAEN,group=1)) + 
  geom_line(lwd=1) + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lowerCI,ymax=upperCI),width=0.02)+
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI),alpha=0.5,fill="lightgray") +
  xlab("duration of exposure (hours)") +
  ylab("probability of reported infection") +
  theme_classic() + 
  geom_point(size=5,col="gray",aes(x=x,y=TPAEN),data = tibble(x=15/60,TPAEN=TPAEN_2m15min)) +
  geom_text(size=4,hjust = 0, nudge_x = -0.06,nudge_y = 0.01,
            aes(x=x,y=TPAEN,label="2m,15min"),col="darkgray",
            data = tibble(x=15/60,TPAEN=TPAEN_2m15min)) +
  scale_y_continuous(limits=c(0,NA)) +
  scale_x_log10(limits=c(0.2,NA)) +
  coord_cartesian(expand=FALSE)
SaveFigure(plot_duration_all,width = 15,height = 12)
ReadTable("table_plot_TPAEN_duration_log")
plot_TPAEN_duration_log_X4c <- ggplot( table_plot_TPAEN_duration_log %>% mutate(month=factor(month,levels=monthshq)),
                                       aes(x=meanscore,y=TPAEN,group=month)) + 
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) +
  xlab("duration of exposure (hours)") +
  ylab("probability of reported infection") +
  theme_classic() +
  scale_x_log10(breaks=c(0.5,1,2,5,10,20)) +
  coord_cartesian(expand=F)
SaveFigure(plot_TPAEN_duration_log_X4c,width = 15,height = 12)
ReadTable("table_plot_TPAEN_duration_log_nobg")
plot_TPAEN_duration_log_nobg <- ggplot( table_plot_TPAEN_duration_log_nobg %>% mutate(month=factor(month,levels=monthshq)),
                                        aes(x=score,y=TPAEN,group=month)) + 
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) +
  xlab("duration of exposure (hours)") +
  ylab("probability of reported infection") +
  theme_classic() +
  scale_x_log10(breaks=c(0.5,1,2,5,10,20)) +
  coord_cartesian(expand=F)
SaveFigure(plot_TPAEN_duration_log_nobg,width = 15,height = 12)
ReadTable("table_plot_TPAEN_duration_log_nofullbg")
plot_TPAEN_duration_log_nofullbg <- ggplot( table_plot_TPAEN_duration_log_nofullbg %>% mutate(month=factor(month,levels=monthshq)),
                                            aes(x=score,y=TPAEN,group=month)) + 
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) +
  xlab("duration of exposure (hours)") +
  ylab("probability of reported infection") +
  theme_classic() +
  scale_x_log10(breaks=c(0.5,1,2,5,10,20)) +
  coord_cartesian(expand=F)
SaveFigure(plot_TPAEN_duration_log_nofullbg,width = 15,height = 12)

ReadTable("table_maxrisk_joint")
ReadTable("table_cumrisk_joint")
ReadTable("table_duration_joint")
#joint plot on TPAEN vs maximum risk score
plot_maxrisk_joint <- ggplot(data=table_maxrisk_joint %>% 
                               mutate(month=factor(month,levels=monthshq,ordered=FALSE)) %>%
                               mutate(name=factor(name,levels=c('with background risk','ML background risk correction','full background risk correction'))),
                             aes(x=meanscore/90,y=TPAEN,group=month)) + 
  geom_hline(yintercept=0,lty=2,col="lightgray") +
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) + 
  xlab("maximum risk score") +
  ylab("probability of reported infection") +
  theme_classic() + 
  scale_colour_viridis_d(option="H") + scale_fill_viridis_d(option="H") + 
  coord_cartesian(expand = F) +
  scale_x_continuous(breaks = c(1,2,3,4,5,10),trans = "log10") + facet_wrap(~name,scales="free_x") 
SaveFigure(plot_maxrisk_joint,width = 20,height = 12)

#joint plot on TPAEN vs cumulative risk score
plot_cumrisk_joint <- ggplot(data=table_cumrisk_joint %>% 
                               mutate(month=factor(month,levels=monthshq,ordered=FALSE)) %>%
                               mutate(name=factor(name,levels=c('with background risk','ML background risk correction','full background risk correction'))),
                             aes(x=meanscore/90,y=TPAEN,group=month)) + 
  geom_hline(yintercept=0,lty=2,col="lightgray") +
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) + 
  xlab("cumulative risk score") +
  ylab("probability of reported infection") +
  theme_classic() + 
  scale_colour_viridis_d(option="H") + scale_fill_viridis_d(option="H") + 
  coord_cartesian(expand = F) +
  scale_x_continuous(trans = "log10") + facet_wrap(~name,scales="free_x") 
SaveFigure(plot_cumrisk_joint,width = 20,height = 12)

#joint plot on TPAEN vs duration
plot_duration_joint <- ggplot(data=table_duration_joint %>% 
                                mutate(month=factor(month,levels=monthshq,ordered=FALSE)) %>%
                                mutate(name=factor(name,levels=c('with background risk','ML background risk correction','full background risk correction'))),
                              aes(x=meanscore,y=TPAEN,group=month)) + 
  geom_hline(yintercept=0,lty=2,col="lightgray") +
  geom_line(aes(col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=month),alpha=0.3) + 
  xlab("total duration of exposure (hours)") +
  ylab("probability of reported infection") +
  theme_classic() + 
  scale_colour_viridis_d(option="H") + scale_fill_viridis_d(option="H") + 
  coord_cartesian(expand = F) +
  scale_x_continuous(breaks=c(0.5,1,3,10,30), trans = "log10") + facet_wrap(~name,scales="free_x") 
SaveFigure(plot_duration_joint,width = 20,height = 12)

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  plot_maxrisk_joint + theme(legend.box.margin = margin(0, 0, 0, 3))
)
plot_bgremoval_all_joint<-plot_grid(plot_grid(
  plot_maxrisk_joint + theme(legend.position="none"),
  plot_cumrisk_joint + theme(legend.position="none"),
  plot_duration_joint + theme(legend.position="none"),ncol=1),
  legend,ncol=2,rel_widths=c(3,0.5))
SaveFigure(plot_bgremoval_all_joint,width = 20,height = 20)


# joint plot on TPAEN vs everything
ReadTable("table_overall_joint")
table_overall_joint <- table_overall_joint %>%
  mutate(name=factor(name,ordered = TRUE,levels = c("maximum risk score","cumulative risk score","duration (hours)")))
plot_TPAEN_combined <- ggplot(data=table_overall_joint  %>%
                                mutate(meanscore=if_else(name!="duration (hours)",meanscore/90,meanscore)),
                              aes(x=meanscore,y=TPAEN)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI),alpha=0.75,fill="lightgray") +
  geom_linerange(aes(ymin=lowerCI,ymax=upperCI))+
  geom_point(size=0.5) +
  geom_line(lwd=0.25) + 
  xlab("score") +
  ylab("probability of reported infection") +
  theme_classic() +
  geom_point(size=3,col="gray",aes(x=x,y=TPAEN),data = tibble(x=c(1,1,15/60),TPAEN=rep(TPAEN_2m15min,3),
                                                              name=factor(c("maximum risk score","cumulative risk score","duration (hours)"),ordered = TRUE,levels = c("maximum risk score","cumulative risk score","duration (hours)")))) +
  geom_text(size=4,hjust = 0, nudge_x = c(-0.03,0.1,0.1),nudge_y = 0.1-0.01,
            aes(x=x,y=TPAEN,label="2m,15min"),col=c("darkgray","white","white"),
            data = tibble(x=c(1,1,15/60),TPAEN=rep(TPAEN_2m15min,3),
                          name=factor(c("maximum risk score","cumulative risk score","duration (hours)"),ordered = TRUE,levels = c("maximum risk score","cumulative risk score","duration (hours)")))) +
  geom_segment( data = tibble(x=c(1,1,15/60),TPAEN=rep(TPAEN_2m15min,3),name=factor(c("maximum risk score","cumulative risk score","duration (hours)"),ordered = TRUE,levels = c("maximum risk score","cumulative risk score","duration (hours)"))),
             aes(xend=x,yend=TPAEN+0.015,x=x,y=TPAEN+0.1-0.025),arrow=arrow(length = unit(0.05, "npc")),lineend="round",linejoin="mitre",
             col=c("darkgray","white","white"),lwd=0.5) +
  scale_y_continuous(expand=expansion(add=0),limits=c(0,NA)) +
  scale_x_continuous(expand=expansion(mult=0.05),trans = "log10") +
  coord_cartesian(expand=TRUE) + facet_wrap(~name,scale="free_x") 
SaveFigure(plot_TPAEN_combined,width = 15,height = 8)

#joint plot on TPAEN vs everything by month
ReadTable("table_all_joint")
table_all_joint <- table_all_joint %>%
  mutate(name=factor(name,ordered = TRUE,levels = c("maximum risk score","cumulative risk score","duration (hours)")))
plot_all_joint <- ggplot(data=table_all_joint %>% mutate(month=factor(month,levels=monthshq,ordered=FALSE),meanscore=if_else(name!="duration (hours)",meanscore/90,meanscore)),
                         aes(x=meanscore,y=TPAEN)) + 
  geom_line(aes(group=month,col=month)) + 
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,group=month,fill=month),alpha=0.3) + 
  xlab("score") +
  ylab("probability of reported infection") +
  theme_classic() + 
  scale_colour_viridis_d(option="H") + scale_fill_viridis_d(option="H") + 
  geom_point(size=3,col="white",aes(x=x,y=TPAEN),data = tibble(x=c(1,1,15/60),TPAEN=rep(TPAEN_2m15min,3),
                                                              name=factor(c("maximum risk score","cumulative risk score","duration (hours)"),ordered = TRUE,levels = c("maximum risk score","cumulative risk score","duration (hours)")))) +
  scale_y_continuous(expand=expansion(add=0),limits=c(0,NA)) +
  scale_x_continuous(expand=expansion(mult=0.05),trans = "log10") + facet_wrap(~name,scales="free_x") 
SaveFigure(plot_all_joint,width = 20,height = 12)

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  plot_all_joint + theme(legend.box.margin = margin(0, 0, 0, 6))
)

#final plot figure 1 (combination of previous 2 plots)
combined_figure1 <- plot_grid(
  plot_grid(plot_TPAEN_combined #+ theme(axis.title.x = element_blank()) 
            ,plot_all_joint + theme(legend.position="none"), ncol = 1,rel_heights = c(0.9,1), align = "hv", axis="lr"),
  legend,ncol=2,rel_widths=c(3,0.5))
SaveFigure(combined_figure1,height=12,width=20)

# plots on heterogeneities

# day of week
table_MLrisk_fig4a<-ReadTable("table_MLrisk_fig4a") %>% 
  mutate(weekday=weekdays(ISOdate(1, 1, 1:7),abbreviate = TRUE)[weekday])
table_MLrisk_fig4a2<-ReadTable("table_MLrisk_fig4a2")
table_MLrisk_fig4a3<-bind_rows(table_MLrisk_fig4a,table_MLrisk_fig4a2) %>% 
  mutate(weekday=factor(weekday,ordered=TRUE,levels=c("Mon-Thu",weekdays(ISOdate(1, 1, 1:7),abbreviate = TRUE)[5:7])))  %>% 
  filter(!is.na(weekday)) 

# plot day of week

# plot all weekend days
plot_MLrisk_fig_dayweekend <- ggplot(data=full_join(table_relation_riskscore_transmissionrisk %>% select(x,riskScore) ,table_MLrisk_fig4a3,multiple="all") %>% filter(riskScore!="background") %>% 
                                       mutate(weekday=factor(weekday)),
                                     aes(x=x,y=transmissionProbability,group=weekday)) +
  geom_line(aes(x=x,y=transmissionProbability,col=weekday)) +
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=weekday),alpha=0.15) +
  theme_classic() +
  ylab("probability of reported transmission\nper exposure window") +
  xlab("risk score") +
  scale_colour_manual(
    values = c(`Mon-Thu`="red",`Fri`="brown",`Sat`="darkgreen",`Sun`="blue"),
    aesthetics = c("colour", "fill")
  ) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,0.026)#,trans='log10'
  ) + scale_x_continuous(expand = c(0,0),limits = c(0,NA))
SaveFigure(plot_MLrisk_fig_dayweekend,width=15,height=12)

# plot all days of week
plot_MLrisk_fig_alldayweek <- ggplot(data=full_join(table_relation_riskscore_transmissionrisk %>% select(x,riskScore) ,table_MLrisk_fig4a,multiple="all") %>% filter(riskScore!="background") %>% 
                                       mutate(weekday=factor(weekday,levels=weekdays(ISOdate(1, 1, 1:7),abbreviate = TRUE))),
                                     aes(x=x,y=transmissionProbability,group=weekday)) +
  geom_line(aes(x=x,y=transmissionProbability,col=weekday)) +
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=weekday),alpha=0.15) +
  theme_classic() +
  ylab("probability of reported transmission\nper exposure window") +
  xlab("risk score") +
  scale_colour_manual(
    values = c(`Mon`="firebrick",`Tue`="red",`Wed`="tomato",`Thu`="orange",`Fri`="brown",`Sat`="darkgreen",`Sun`="blue"),
    aesthetics = c("colour", "fill")
  ) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,0.026)#,trans='log10'
  ) + scale_x_continuous(expand = c(0,0),limits = c(0,NA))
SaveFigure(plot_MLrisk_fig_alldayweek,width=15,height=12)


# final plot day of week
ReadTable("table_MLrisk_fig4a2")
plot_MLrisk_fig_dayweek <- ggplot(data=full_join(table_relation_riskscore_transmissionrisk %>% select(x,riskScore) ,table_MLrisk_fig4a2,multiple="all") %>% filter(riskScore!="background") %>% 
                                    mutate(weekday=factor(weekday)),
                                  aes(x=x,y=transmissionProbability,group=weekday)) +
  geom_line(aes(x=x,y=transmissionProbability,col=weekday)) +
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=weekday),alpha=0.3) +
  theme_classic() +
  ylab("probability of reported transmission\nper exposure window") +
  xlab("risk score") +
  scale_colour_manual(
    values = c(`Mon-Thu`="blue",`Fri-Sun`="red"),
    aesthetics = c("colour", "fill")
  ) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,0.026)#,trans='log10'
  ) + scale_x_continuous(expand = c(0,0),limits = c(0,NA))
SaveFigure(plot_MLrisk_fig_dayweek,width=15,height=12)

# final plot risk profile by rural/urban score:
ReadTable("table_MLrisk_fig_urban3")
plot_MLrisk_fig_ruralurban <- ggplot(data=full_join(table_relation_riskscore_transmissionrisk %>% select(x,riskScore) ,table_MLrisk_fig_urban3,multiple="all") %>% filter(riskScore!="background") %>% 
                                       #filter(location %in% list_locations[1:5]) %>%
                                       mutate(area=factor(location,levels=c("rural","urban","conurbation"))),
                                     aes(x=x,y=transmissionProbability,group=area)) +
  geom_line(aes(x=x,y=transmissionProbability,col=area)) +
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=area),alpha=0.3) +
  theme_classic() +
  ylab("probability of reported transmission\nper exposure window") +
  xlab("risk score") +
  scale_color_viridis_d(option="E",direction = -1) + scale_fill_viridis_d(option="E",direction = -1) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,0.026)#,trans='log10'
                     ) + scale_x_continuous(expand = c(0,0),limits = c(0,NA))
SaveFigure(plot_MLrisk_fig_ruralurban,width=15,height=12)

# plot risk profile by rural/urban score:
table_MLrisk_fig_urban <- ReadTable("table_MLrisk_fig_urban")
plot_MLrisk_figX2c <- ggplot(data=full_join(table_relation_riskscore_transmissionrisk %>% select(x,riskScore) ,table_MLrisk_fig_urban,multiple="all") %>% filter(riskScore!="background") %>% 
                               #filter(location %in% list_locations[1:5]) %>%
                               mutate(location=factor(location,levels=c("Wales",1:4,"5/6","London"))),
                             aes(x=x,y=transmissionProbability,group=location)) +
  geom_line(aes(x=x,y=transmissionProbability,col=location)) +
  geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=location),alpha=0.1) +
  theme_classic() +
  ylab("probability of reported transmission\nper exposure window") +
  xlab("risk score") +
  scale_color_viridis_d() + scale_fill_viridis_d() + 
  scale_y_continuous(expand = c(0,0),limits = c(0,0.026)#,trans='log10'
  ) + scale_x_continuous(expand = c(0,0),limits = c(0,NA))
SaveFigure(plot_MLrisk_figX2c,width=15,height=12)


# final joint plot heterogeneities
joint_plot_weekday_ruralurban <- plot_grid(plot_MLrisk_fig_dayweek + 
                                             theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                                   legend.position = c(0.2, 0.8)),
                                           plot_MLrisk_fig_ruralurban + 
                                             theme(legend.position = c(0.2, 0.8), axis.text.y = element_blank(),
                                                   axis.ticks.y = element_blank(),
                                                   axis.title.y = element_blank(), 
                                                   axis.line.y = element_blank(),
                                                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)),
                                           nrow=1, labels = "auto", align = "h", vjust = +1,hjust=-1
)
SaveFigure(joint_plot_weekday_ruralurban,width=20,height=12)

p1<-plot_relation_riskscore_transmissionrisk
p2<-plot_MLrisk_fig_dayweek + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = c(0.2, 0.8))
p3<-plot_MLrisk_fig_ruralurban + 
  theme(legend.position = c(0.25, 0.8), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plots <- align_plots(p1, p2, align = 'v', axis = 'l')
# then build the bottom row
bottom_row <- plot_grid(plots[[2]], p3, align = "h", labels = c('b', 'c'), label_size = 18,hjust=-2,vjust=1.0,rel_widths=c(12,10))
# then combine with the top row for final plot
joint_plot_perwindow <- plot_grid(plots[[1]], bottom_row, labels = c('a', ''), label_size = 18,hjust=-2 , ncol = 1,  rel_heights = c(2,1))
SaveFigure(joint_plot_perwindow,width=20,height=24)


# population-level distributions of risk score, duration, etc

# heterogeneities in weekday, month, area
distribution_duration_bydayweek <- ggplot(table_outcome_riskScores %>% mutate(weekday=dayweek2short(peak_day)), aes(x=duration/3600,group=weekday))+
  geom_density(aes(color=weekday),bw=0.1,alpha=0.3) + scale_x_log10() + xlab("duration (hours)") + 
  theme_classic() + coord_cartesian(xlim=c(0.3,NA),expand = FALSE)
SaveFigure(distribution_duration_bydayweek,width=12,height=10)

distribution_duration_bymonth <- ggplot(table_outcome_riskScores %>% mutate(month=factor(month,levels = monthshq)), aes(x=duration/3600,group=month))+
  geom_density(aes(color=month),bw=0.1,alpha=0.3) + scale_x_log10() + xlab("duration (hours)") + 
  theme_classic() + coord_cartesian(xlim=c(0.3,NA),expand = FALSE)
SaveFigure(distribution_duration_bymonth,width=12,height=10)

distribution_duration_byarea <- ggplot(table_outcome_riskScores %>% mutate(area=factor(c("London",as.character(1:6),"Wales")[RuralUrbanScore+1],levels=c("Wales",as.character(1:6),"London"))), aes(x=duration/3600,group=area))+
  geom_density(aes(color=area),bw=0.1,alpha=0.3) + scale_x_log10() + xlab("duration (hours)") + 
  theme_classic() + coord_cartesian(xlim=c(0.3,NA),expand = FALSE)
SaveFigure(distribution_duration_byarea,width=12,height=10)

distribution_maxrisk_bydayweek <- ggplot(table_outcome_riskScores %>% mutate(weekday=dayweek2short(peak_day)), aes(x=maxRiskScore/90,group=weekday))+
  geom_density(aes(color=weekday),bw=0.01,alpha=0.3) + scale_x_continuous(breaks = c(1.5,2,3,5,10),trans="log10") + xlab("maximum risk score") + 
  theme_classic() + coord_cartesian(expand = FALSE)
SaveFigure(distribution_maxrisk_bydayweek,width=12,height=10)

distribution_maxrisk_bymonth <- ggplot(table_outcome_riskScores %>% mutate(month=factor(month,levels = monthshq)), aes(x=maxRiskScore/90,group=month))+
  geom_density(aes(color=month),bw=0.01,alpha=0.3) + scale_x_continuous(breaks = c(1.5,2,3,5,10),trans="log10") + xlab("maximum risk score") + 
  theme_classic() + coord_cartesian(expand = FALSE)
SaveFigure(distribution_maxrisk_bymonth,width=12,height=10)

distribution_maxrisk_byarea <- ggplot(table_outcome_riskScores %>% mutate(area=factor(c("London",as.character(1:6),"Wales")[RuralUrbanScore+1],levels=c("Wales",as.character(1:6),"London"))), aes(x=maxRiskScore/90,group=area))+
  geom_density(aes(color=area),bw=0.01,alpha=0.3) + scale_x_continuous(breaks = c(1.5,2,3,5,10),trans="log10") + xlab("maximum risk score") + 
  theme_classic() + coord_cartesian(expand = FALSE)
SaveFigure(distribution_maxrisk_byarea,width=12,height=10)

distribution_cumrisk_bydayweek <- ggplot(table_outcome_riskScores %>% mutate(weekday=dayweek2short(peak_day)), aes(x=cumRiskScore/90,group=weekday))+
  geom_density(aes(color=weekday),bw=0.1,alpha=0.3) + scale_x_continuous(trans="log10") + xlab("cumulative risk score") + 
  theme_classic() + coord_cartesian(expand = FALSE)
SaveFigure(distribution_cumrisk_bydayweek,width=12,height=10)

distribution_cumrisk_bymonth <- ggplot(table_outcome_riskScores %>% mutate(month=factor(month,levels = monthshq)), aes(x=cumRiskScore/90,group=month))+
  geom_density(aes(color=month),bw=0.1,alpha=0.3) + scale_x_continuous(trans="log10") + xlab("cumulative risk score") + 
  theme_classic() + coord_cartesian(expand = FALSE)
SaveFigure(distribution_cumrisk_bymonth,width=12,height=10)

distribution_cumrisk_byarea <- ggplot(table_outcome_riskScores %>% mutate(area=factor(c("London",as.character(1:6),"Wales")[RuralUrbanScore+1],levels=c("Wales",as.character(1:6),"London"))), aes(x=cumRiskScore/90,group=area))+
  geom_density(aes(color=area),bw=0.1,alpha=0.3) + scale_x_continuous(trans="log10") + xlab("cumulative risk score") + 
  theme_classic() + coord_cartesian(expand = FALSE)
SaveFigure(distribution_cumrisk_byarea,width=12,height=10)

# binned distribution duration for transmissions
discrete_distribution_transmission_duration_corrections <- ggplot(table_outcome_riskScores %>% mutate(duration_hours=cut(duration/3600,breaks=2^c(-5,-1:6,10),labels=c("<0.5","0.5-1","1-2","2-4","4-8","8-16","16-32","32-64",">64")
)) %>% group_by(duration_hours) %>% 
  summarise(counts_none=sum(positive-bg_rate_cases_app),counts_fair=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml)),counts_all=sum(positive)), aes(x=duration_hours))+
  geom_col(aes(y=counts_all,fill="no correction")) + 
  geom_col(aes(y=counts_fair,fill="ML background")) + 
  geom_col(aes(y=counts_none,fill="full background")) + 
  xlab("duration (hours)") + ylab("transmissions") +
  theme_classic() + 
  scale_x_discrete(expand = c(0,0)) +
  labs(fill = "Background correction") +
  scale_fill_manual(values = c(`no correction`=alpha("orange",alpha=0.1),`ML background`=alpha("red",alpha=0.3),`full background`=alpha("darkred",alpha=0.5))) +
  theme(legend.position = c(0.6,0.8))
SaveFigure(discrete_distribution_transmission_duration_corrections,width=15,height=12)

discrete_distribution_transmission_duration <- ggplot(table_outcome_riskScores %>% mutate(duration_hours=cut(duration/3600,breaks=2^c(-5,-1:6,10),labels=c("<0.5","0.5-1","1-2","2-4","4-8","8-16","16-32","32-64",">64")
)) %>% group_by(duration_hours) %>% 
  summarise(counts_mid=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml)),counts_lower=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml_upper)),counts_upper=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml_lower))), 
aes(x=duration_hours))+
  geom_col(aes(y=counts_mid),alpha=0.3,fill="red") + 
  geom_errorbar(aes(ymin=counts_lower, ymax=counts_upper), width=.2,
                position=position_dodge(.9)) +
  xlab("duration (hours)") + ylab("transmissions") +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))

SaveFigure(discrete_distribution_transmission_duration,width=15,height=12)
discrete_distribution_transmission_duration_joint <- plot_grid(
  discrete_distribution_transmission_duration_corrections,
  discrete_distribution_transmission_duration,
  nrow=1, labels = "auto", align = "h", vjust = +1,hjust=-1
)
SaveFigure(discrete_distribution_transmission_duration_joint,width=25,height=12)

# all distributions
table_population_duration_corr <- table_outcome_riskScores %>% mutate(duration_hours=cut(duration/3600,breaks=2^c(-7,-1:6,10),labels=c("<0.5","0.5-1","1-2","2-4","4-8","8-16","16-32","32-64",">64")
)) %>% group_by(duration_hours) %>% 
  summarise(counts_all=n(),counts_pos=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml),na.rm=T)) %>%
  mutate(`all contacts`=counts_all/sum(counts_all),`contacts testing positive`=counts_pos/sum(counts_pos))
SaveTable(table_population_duration_corr)
table_population_duration_raw <- table_outcome_riskScores %>% mutate(duration_hours=cut(duration/3600,breaks=2^c(-7,-1:6,10),labels=c("<0.5","0.5-1","1-2","2-4","4-8","8-16","16-32","32-64",">64")
)) %>% group_by(duration_hours) %>% 
  summarise(counts_all=n(),counts_pos=sum(positive)) %>%
  mutate(`all contacts`=counts_all/sum(counts_all),`contacts testing positive`=counts_pos/sum(counts_pos))
SaveTable(table_population_duration_raw)
mycols<-c(`all`="blue",`transmissions`="red")
plot_distr_dur_corr <- ggplot(table_population_duration_corr, 
aes(x=duration_hours))+
  geom_col(aes(y=`all contacts`,fill="all"),alpha=0.3) + 
  geom_col(aes(y=`contacts testing positive`,fill="transmissions"),alpha=0.3) +
  scale_fill_manual(values=mycols) +
  labs(x="duration (hours)",y="fraction of contacts",fill="Contacts:") +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
mycols<-c(`all`="blue",`testing positive`="darkred")
plot_distr_dur_raw <- ggplot(table_population_duration_raw,
  aes(x=duration_hours))+
  geom_col(aes(y=`all contacts`,fill="all"),alpha=0.3) + 
  geom_col(aes(y=`contacts testing positive`,fill="testing positive"),alpha=0.3) +
  scale_fill_manual(values=mycols) +
  labs(x="duration (hours)",y="fraction of contacts",fill="Contacts:") +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))

mycols<-c(`all`="blue",`transmissions`="red")
plot_distr_cum_corr <- ggplot(table_outcome_riskScores %>% mutate(duration_hours=cut(cumRiskScore/90,breaks=2^c(0:8,15),labels=c("1-2","2-4","4-8","8-16","16-32","32-64","64-128","128-256",">256")
)) %>% group_by(duration_hours) %>% 
  summarise(counts_all=n(),counts_pos=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml))) %>%
  mutate(`all contacts`=counts_all/sum(counts_all),`contacts testing positive`=counts_pos/sum(counts_pos)), 
aes(x=duration_hours))+
  geom_col(aes(y=`all contacts`,fill="all"),alpha=0.3) + 
  geom_col(aes(y=`contacts testing positive`,fill="transmissions"),alpha=0.3) +
  scale_fill_manual(values=mycols) +
  labs(x="cumulative risk score",y="fraction of contacts",fill="Contacts:") +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
mycols<-c(`all`="blue",`testing positive`="darkred")
plot_distr_cum_raw <- ggplot(table_outcome_riskScores %>% mutate(duration_hours=cut(cumRiskScore/90,breaks=2^c(0:8,15),labels=c("1-2","2-4","4-8","8-16","16-32","32-64","64-128","128-256",">256")
)) %>% group_by(duration_hours) %>% 
  summarise(counts_all=n(),counts_pos=sum(positive)) %>%
  mutate(`all contacts`=counts_all/sum(counts_all),`contacts testing positive`=counts_pos/sum(counts_pos)), 
aes(x=duration_hours))+
  geom_col(aes(y=`all contacts`,fill="all"),alpha=0.3) + 
  geom_col(aes(y=`contacts testing positive`,fill="testing positive"),alpha=0.3) +
  scale_fill_manual(values=mycols) +
  labs(x="cumulative risk score",y="fraction of contacts",fill="Contacts:") +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))

mycols<-c(`all`="blue",`transmissions`="red")
plot_distr_max_corr <- ggplot(table_outcome_riskScores %>% mutate(duration_hours=cut(maxRiskScore,breaks=c(100*1.25^c(0:12)),ordered_result = TRUE)
) %>% group_by(duration_hours) %>% 
  summarise(counts_all=n(),counts_pos=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml))) %>%
  mutate(`all contacts`=counts_all/sum(counts_all),`contacts testing positive`=counts_pos/sum(counts_pos)), 
aes(x=duration_hours))+
  geom_col(aes(y=`all contacts`,fill="all"),alpha=0.3) + 
  geom_col(aes(y=`contacts testing positive`,fill="transmissions"),alpha=0.3) +
  scale_fill_manual(values=mycols) +
  labs(x="maximum risk score",y="fraction of contacts",fill="Contacts:") +
  scale_x_discrete(
    labels=c("1.25","","2","","3","","","6","","","","12","",""),expand = c(0,0)) +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0))
mycols<-c(`all`="blue",`testing positive`="darkred")
plot_distr_max_raw <- ggplot(table_outcome_riskScores %>% mutate(duration_hours=cut(maxRiskScore,breaks=c(100*1.25^c(0:12)),ordered_result = TRUE)
) %>% group_by(duration_hours) %>% 
  summarise(counts_all=n(),counts_pos=sum(positive)) %>%
  mutate(`all contacts`=counts_all/sum(counts_all),`contacts testing positive`=counts_pos/sum(counts_pos)), 
aes(x=duration_hours))+
  geom_col(aes(y=`all contacts`,fill="all"),alpha=0.3) + 
  geom_col(aes(y=`contacts testing positive`,fill="testing positive"),alpha=0.3) +
  scale_fill_manual(values=mycols) +
  labs(x="maximum risk score",y="fraction of contacts",fill="Contacts:") +
  scale_x_discrete(
    labels=c("1.25","","2","","3","","","6","","","","12","",""),expand = c(0,0)) +
  theme_classic() + 
  scale_y_continuous(expand = c(0,0))

plot_distr_joint <- plot_grid(
  plot_distr_max_raw + theme(legend.position="none"),
  plot_distr_cum_raw + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position="none"),
  plot_distr_dur_raw + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position=c(0.6,0.7)),
  plot_distr_max_corr + theme(legend.position="none"),
  plot_distr_cum_corr + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position="none"),
  plot_distr_dur_corr + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position=c(0.6,0.7)),
  nrow=2, align = "h"#, labels = "auto",  vjust = +1,hjust=0
)
SaveFigure(plot_distr_joint,width=20,height=12)


# main figure with distributions
# and effectiveness of intervention by threshold on duration
mybw<-0.025
overall_density <- density(table_outcome_riskScores %>% mutate(x=log(maxRiskScore/90)) %>% pull(x), weights =rep(1/dim(table_outcome_riskScores)[1],dim(table_outcome_riskScores)[1]), bw=mybw, cut=0)
bg_density <- density(table_outcome_riskScores %>% mutate(x=log(maxRiskScore/90)) %>% pull(x),weights = 1-(1-table_outcome_riskScores$bg_rate_cases_app)^rho_ml, bw=mybw, cut=0)
positive_density <- density(table_outcome_riskScores %>% mutate(x=log(maxRiskScore/90)) %>% pull(x),weights = table_outcome_riskScores$positive, bw=mybw, cut=0)
positive_density_bg_removed <- positive_density
positive_density_bg_removed$y <- positive_density$y - bg_density$y
positive_density_bg_removed$y <- positive_density_bg_removed$y/sum(positive_density_bg_removed$y)*sum(overall_density$y)

mycols<-c(`contacts`="blue",`transmissions`="red")
plot_density_maxrisk <-ggplot(data=as_tibble(overall_density[c("x","y")]), aes(x=exp(x)))+
  geom_ribbon(aes(ymin=0,ymax=y,fill="contacts"),col="blue",alpha=0.3) + 
  geom_ribbon(data=as_tibble(positive_density_bg_removed[c("x","y")]),aes(ymin=0,ymax=y,fill="transmissions"),col="red",alpha=0.3) + 
  scale_x_continuous(breaks = c(1.5,2,3,5,10),trans="log10",labels = prettyNum) + 
  xlab("maximum risk score") + 
  ylab("density") +
  theme_classic() +
  scale_fill_manual(values=mycols) +
  guides(fill = guide_legend(
    override.aes = list(color = "white")
  )) +
  coord_cartesian(expand = FALSE)

mylog<-function(x){x*(x<=10)+10*(1+log(x/10))*(x>10)}
myderiv<-function(x){x*(x<=10)+10*(x>10)}
myexp<-function(x){x*(x<=10)+10*exp(x/10-1)*(x>10)}

mybw<-0.6
overall_density <- density(table_outcome_riskScores %>% mutate(x=mylog(cumRiskScore/90)) %>% pull(x), n=10000, weights =rep(1/dim(table_outcome_riskScores)[1],dim(table_outcome_riskScores)[1]), bw=mybw, cut=0)
bg_density <- density(table_outcome_riskScores %>% mutate(x=mylog(cumRiskScore/90)) %>% pull(x), n=10000,weights = 1-(1-table_outcome_riskScores$bg_rate_cases_app)^rho_ml, bw=mybw, cut=0)
positive_density <- density(table_outcome_riskScores %>% mutate(x=mylog(cumRiskScore/90)) %>% pull(x), n=10000,weights = table_outcome_riskScores$positive, bw=mybw, cut=0)
positive_density_bg_removed <- positive_density
positive_density_bg_removed$y <- positive_density$y - bg_density$y
positive_density_bg_removed$y <- positive_density_bg_removed$y/sum(positive_density_bg_removed$y)*sum(overall_density$y)
overall_density$x <- overall_density$x + mybw/2
positive_density_bg_removed$x <- positive_density_bg_removed$x + mybw/2
overall_density$y <- overall_density$y*myderiv(overall_density$x)
positive_density_bg_removed$y <- positive_density_bg_removed$y*myderiv(positive_density_bg_removed$x)
overall_density$y <- overall_density$y/sum((overall_density$y[-1]+overall_density$y[-length(overall_density$y)])/2*diff(log10(myexp(overall_density$x))))
positive_density_bg_removed$y <- positive_density_bg_removed$y/sum((positive_density_bg_removed$y[-1]+positive_density_bg_removed$y[-length(positive_density_bg_removed$y)])/2*diff(log10(myexp(positive_density_bg_removed$x))))


plot_density_cumrisk <-ggplot(data=as_tibble(overall_density[c("x","y")]), aes(x=myexp(x)))+
  geom_ribbon(aes(ymin=0,ymax=y),fill="blue",col="blue",alpha=0.3) + 
  geom_ribbon(data=as_tibble(positive_density_bg_removed[c("x","y")]),aes(ymin=0,ymax=y),fill="red",col="red",alpha=0.3) + 
  scale_x_continuous(breaks = c(2,5,10,20,50,100,200,1000),trans="log10",labels = prettyNum) + 
  xlab("cumulative risk score") + 
  ylab("density") +
  theme_classic() +
  coord_cartesian(expand = FALSE)

mybw<-0.3
overall_density <- density(table_outcome_riskScores %>% mutate(x=mylog(pmax(duration/3600,0.5))) %>% pull(x), n=10000, weights =rep(1/dim(table_outcome_riskScores)[1],dim(table_outcome_riskScores)[1]), bw=mybw, cut=0)
bg_density <- density(table_outcome_riskScores %>% mutate(x=mylog(pmax(duration/3600,0.5))) %>% pull(x), n=10000, weights = 1-(1-table_outcome_riskScores$bg_rate_cases_app)^rho_ml, bw=mybw, cut=0)
positive_density <- density(table_outcome_riskScores %>% mutate(x=mylog(pmax(duration/3600,0.5))) %>% pull(x), n=10000, weights = table_outcome_riskScores$positive, bw=mybw, cut=0)
positive_density_bg_removed <- positive_density
positive_density_bg_removed$y <- positive_density$y - bg_density$y
positive_density_bg_removed$y <- positive_density_bg_removed$y/sum(positive_density_bg_removed$y)*sum(overall_density$y)
overall_density$x <- overall_density$x + mybw/2
positive_density_bg_removed$x <- positive_density_bg_removed$x + mybw/2
overall_density$y <- overall_density$y*myderiv(overall_density$x)
positive_density_bg_removed$y <- positive_density_bg_removed$y*myderiv(positive_density_bg_removed$x)
overall_density$y <- overall_density$y/sum((overall_density$y[-1]+overall_density$y[-length(overall_density$y)])/2*diff(log10(myexp(overall_density$x))))
positive_density_bg_removed$y <- positive_density_bg_removed$y/sum((positive_density_bg_removed$y[-1]+positive_density_bg_removed$y[-length(positive_density_bg_removed$y)])/2*diff(log10(myexp(positive_density_bg_removed$x))))

mycols<-c(`contacts`="blue",`transmissions`="red")
plot_density_duration <- ggplot(data=as_tibble(overall_density[c("x","y")]), aes(x=myexp(x)))+
  geom_ribbon(aes(ymin=0,ymax=y),fill="blue",col="blue",alpha=0.3) + 
  geom_ribbon(data=as_tibble(positive_density_bg_removed[c("x","y")]),aes(ymin=0,ymax=y),fill="red",col="red",alpha=0.3) + 
  scale_x_continuous(breaks = c(0.1,0.2,0.5,1,2,5,10,20,50),trans="log10",labels = prettyNum) + 
  xlab("duration of exposure\n(hours)") +
  ylab("density") +
  theme_classic() +
  coord_cartesian(expand = FALSE)

plot_outcome_duration <- ggplot(table_outcome_riskScores %>% mutate(duration_hours=cut(duration/3600,breaks=2^c(-5,-1:6,10),labels=c("none",">1/2h",">1h",">2h",">4h",">8h",">16h",">32h",">64h"),ordered_result = TRUE
)) %>% group_by(duration_hours) %>% 
  summarise(counts_mid=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml)),counts_lower=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml_upper)),counts_upper=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml_lower))) %>% 
  arrange(desc(duration_hours)) %>% mutate(across(starts_with("counts"),~ cumsum(.x)/sum(.x))), 
aes(x=duration_hours))+
  geom_col(aes(y=counts_upper*100),alpha=0.3,fill="darkred") + 
  geom_col(aes(y=counts_mid*100),alpha=0.3,fill="red") + 
  geom_col(aes(y=counts_lower*100),alpha=0.3,fill="orange") + 
  xlab("threshold on duration") + ylab("tracing effectiveness \n (% of transmissions detected)") +
  theme_classic() + 
  coord_cartesian(expand = FALSE)

plot_SAR_duration <- ggplot(table_outcome_riskScores %>% mutate(duration_hours=cut(duration/3600,breaks=2^c(-5,-1:6,10),labels=c("none",">1/2h",">1h",">2h",">4h",">8h",">16h",">32h",">64h"),ordered_result = TRUE
)) %>% group_by(duration_hours) %>% 
  summarise(counts_mid=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml)),counts_lower=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml_lower)),counts_upper=sum(positive-(1-(1-bg_rate_cases_app)^rho_ml_upper)),counts_all=n()) %>% 
  arrange(desc(duration_hours)) %>% mutate(across(starts_with("counts"),~ cumsum(.x))) %>%
  mutate(sar_mid=counts_mid/counts_all,sar_lower=counts_lower/counts_all,sar_upper=counts_upper/counts_all), 
aes(x=duration_hours))+
  geom_col(aes(y=sar_upper*100),alpha=0.3,fill="darkred") + 
  geom_col(aes(y=sar_mid*100),alpha=0.3,fill="red") + 
  geom_col(aes(y=sar_lower*100),alpha=0.3,fill="orange") + 
  xlab("threshold on duration") + ylab("% of contacts testing positive\ndue to transmission") +
  theme_classic() + 
  coord_cartesian(expand = FALSE)

# plot classification
ReadTable("table_classification")
plot_fraction_classification <- table_classification %>% 
  select(setting,starts_with("fraction")) %>% 
  pivot_longer(cols = starts_with("fraction"),names_to = "statistic", values_to = "value") %>%
  mutate(statistic=str_replace(str_replace(str_replace(statistic,"fraction_of_",""),"_"," ")," ","\n")) %>%
  mutate(`type of\ncontact`=factor(setting,levels=c("household","recurring","single day","fleeting"))) %>%
  ggplot(aes(x=statistic,y=value,group=`type of\ncontact`)) +
  geom_col(aes(fill=`type of\ncontact`)) +
  xlab("") +
  ylab("fraction") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
SaveFigure(plot_fraction_classification,height=10,width=10)

plot_distr_joint_main <- plot_grid(
  plot_density_maxrisk + theme(legend.title=element_blank(),legend.position=c(0.7,0.8)),
  plot_density_duration + theme(legend.position="none"),
  plot_density_cumrisk + theme(legend.position="none"),
  plot_fraction_classification,
  plot_outcome_duration,
  plot_SAR_duration,
  nrow=3, labels = "auto", align = "hv", axis="lb", vjust = 2,hjust=0
)
SaveFigure(plot_distr_joint_main,width=20,height=30)



# plots for TPAEN vs duration and proximity
results_maxproximity_duration <- table_outcome_riskScores %>% 
  select(month,duration,maxProximityScore,avgProximityScore,bg_rate_cases_app,positive) %>%
  filter(month %in% monthshq) %>% 
  mutate(binduration=cut(duration/3600,
                         breaks = c(0,0.5,1,2,3,10,Inf),
                         labels=paste(c(0,0.5,1,2,3,">"),c(rep("-",5),""),c(0.5,1,2,3,10,10),sep="")),
         binmaxprox=cut((maxProximityScore),
                        breaks = rev(0.01*round(1/c(0,2.5,2.75,3,3.25,3.5,Inf)^2*100)),
                        labels=paste(c("",rep(".",4),""),rev(c(">",4*round(1/c(2.75,3,3.25,3.5,Inf)^2*100))),c(rep("-",5),""),".",rev(4*round(1/c(2.5,2.5,2.75,3,3.25,3.5)^2*100)),sep="")),
         binavgprox=cut((avgProximityScore),
                        breaks = rev(0.01*round(1/c(0,2.5,2.75,3,3.25,3.5,Inf)^2*100)),
                        labels=paste(c("",rep(".",4),""),rev(c(">",4*round(1/c(2.75,3,3.25,3.5,Inf)^2*100))),c(rep("-",5),""),".",rev(4*round(1/c(2.5,2.5,2.75,3,3.25,3.5)^2*100)),sep=""))) %>%
  group_by(binduration,binmaxprox) %>%
  summarise(n=n(),TPAEN=mean(positive-(1-(1-bg_rate_cases_app)^rho_ml))) %>%
  mutate(lowerCI=qbeta(0.025,round(TPAEN*n),n-round(TPAEN*n)+1),upperCI=qbeta(1-0.025,round(TPAEN*n)+1,n-round(TPAEN*n)))
results_avgproximity_duration <- table_outcome_riskScores %>% 
  select(month,duration,maxProximityScore,avgProximityScore,bg_rate_cases_app,positive) %>%
  filter(month %in% monthshq) %>% 
  mutate(binduration=cut(duration/3600,
                         breaks = c(0,0.5,1,2,3,10,Inf),
                         labels=paste(c(0,0.5,1,2,3,">"),c(rep("-",5),""),c(0.5,1,2,3,10,10),sep="")),
         binmaxprox=cut((maxProximityScore),
                        breaks = rev(0.01*round(1/c(0,2.5,2.75,3,3.25,3.5,Inf)^2*100)),
                        labels=paste(c("",rep(".",4),""),rev(c(">",4*round(1/c(2.75,3,3.25,3.5,Inf)^2*100))),c(rep("-",5),""),".",rev(4*round(1/c(2.5,2.5,2.75,3,3.25,3.5)^2*100)),sep="")),
         binavgprox=cut((avgProximityScore),
                        breaks = rev(0.01*round(1/c(0,2.5,2.75,3,3.25,3.5,Inf)^2*100)),
                        labels=paste(c("",rep(".",4),""),rev(c(">",4*round(1/c(2.75,3,3.25,3.5,Inf)^2*100))),c(rep("-",5),""),".",rev(4*round(1/c(2.5,2.5,2.75,3,3.25,3.5)^2*100)),sep=""))) %>%
  group_by(binduration,binavgprox) %>%
  summarise(n=n(),TPAEN=mean(positive-(1-(1-bg_rate_cases_app)^rho_ml))) %>%
  mutate(lowerCI=qbeta(0.025,round(TPAEN*n),n-round(TPAEN*n)+1),upperCI=qbeta(1-0.025,round(TPAEN*n)+1,n-round(TPAEN*n)))

miny<-min(c(results_avgproximity_duration$lowerCI,results_maxproximity_duration$lowerCI))/1.025
maxy<-max(c(results_avgproximity_duration$upperCI,results_maxproximity_duration$upperCI))*1.025
plot_duration_avgproximity <- results_avgproximity_duration %>% ggplot(aes(x=binduration,y=TPAEN,group=binavgprox)) + geom_line(aes(col=binavgprox)) + geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=binavgprox),alpha=0.3) + ylab("probability of\nreported transmission") + # ylab("probability of detected transmission") +
  scale_y_continuous(limits=c(miny,maxy)) + coord_cartesian(expand=FALSE) + scale_fill_viridis_d(direction=-1,option="magma") + scale_colour_viridis_d(direction=-1,option="magma") + theme_classic() + xlab("duration (hours)") + facet_grid(. ~ "by average proximity score")
plot_duration_maxproximity <- results_maxproximity_duration %>% ggplot(aes(x=binduration,y=TPAEN,group=binmaxprox)) + geom_line(aes(col=binmaxprox)) + geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=binmaxprox),alpha=0.3) + ylab("probability of\nreported transmission") + # ylab("probability of detected transmission") +
  scale_y_continuous(limits=c(miny,maxy)) + coord_cartesian(expand=FALSE) + scale_fill_viridis_d(direction=-1,option="magma") + scale_colour_viridis_d(direction=-1,option="magma") + theme_classic() + xlab("duration (hours)") + facet_grid(. ~ "by maximum proximity score")
plot_avgproximity_duration <- results_avgproximity_duration %>% ggplot(aes(x=binavgprox,y=TPAEN,group=binduration)) + geom_line(aes(col=binduration)) + geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=binduration),alpha=0.3) + ylab("probability of\nreported transmission") + # ylab("probability of detected transmission") + 
  scale_y_continuous(limits=c(miny,maxy)) + coord_cartesian(expand=FALSE) + scale_fill_viridis_d(direction=-1) + scale_colour_viridis_d(direction=-1) + theme_classic() + xlab("average proximity score") + facet_grid(. ~ "by duration")
plot_maxproximity_duration <- results_maxproximity_duration %>% ggplot(aes(x=binmaxprox,y=TPAEN,group=binduration)) + geom_line(aes(col=binduration)) + geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=binduration),alpha=0.3) + ylab("probability of\nreported transmission") + # ylab("probability of detected transmission") + 
  scale_y_continuous(limits=c(miny,maxy)) + coord_cartesian(expand=FALSE) + scale_fill_viridis_d(direction=-1) + scale_colour_viridis_d(direction=-1) + theme_classic() + xlab("maximum proximity score") + facet_grid(. ~ "by duration")
legend_proximity <- cowplot::get_legend(
  # create some space to the left of the legend
  plot_duration_avgproximity +
    guides(fill=guide_legend(title="proximity\n score"),color=guide_legend(title="proximity\n score")) + 
    theme(legend.box.margin = margin(0, 0, 0, 0))
)
legend_duration <- cowplot::get_legend(
  # create some space to the left of the legend
  plot_avgproximity_duration + 
    guides(fill=guide_legend(title="duration\n (hours)"),color=guide_legend(title="duration\n (hours)")) + 
    theme(legend.box.margin = margin(0, 0, 0, 0))
)

plot_main_proximity_duration <- plot_grid(plot_duration_avgproximity + theme(legend.position="none",plot.margin =margin(t=5,r=15,b=5,l=5)) + scale_y_log10(limits=c(miny,maxy)),
                                          plot_duration_maxproximity + scale_y_log10(limits=c(miny,maxy)) + 
                                            theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5),
                                                  axis.title.y = element_blank()
                                            ),
                                          legend_proximity,
                                          plot_avgproximity_duration + theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5)) + scale_y_log10(limits=c(miny,maxy)),
                                          plot_maxproximity_duration + scale_y_log10(limits=c(miny,maxy)) +
                                            theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5),
                                                  axis.title.y = element_blank()
                                            ),
                                          legend_duration,
                                          nrow=2,align="hv",axis="lb",rel_widths=c(1,1,0.4))
SaveFigure(plot_main_proximity_duration,width=20,height=16)
plot_main_proximity_duration_log <- plot_grid(plot_duration_avgproximity + theme(legend.position="none",plot.margin =margin(t=5,r=15,b=5,l=5)) ,#+ scale_y_log10(limits=c(miny,maxy)),
                                              plot_duration_maxproximity + #scale_y_log10(limits=c(miny,maxy)) + 
                                                theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5),
                                                      axis.title.y = element_blank()
                                                ),
                                              legend_proximity,
                                              plot_avgproximity_duration + theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5)) ,#+ scale_y_log10(limits=c(miny,maxy)),
                                              plot_maxproximity_duration + #scale_y_log10(limits=c(miny,maxy)) +
                                                theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5),
                                                      axis.title.y = element_blank()
                                                ),
                                              legend_duration,
                                              nrow=2,align="hv",axis="lb",rel_widths=c(1,1,0.4))
SaveFigure(plot_main_proximity_duration_log,width=20,height=16)


# plots for TPAEN vs duration and riskscore
results_maxriskscore_duration <- table_outcome_riskScores %>% mutate(avgRiskScore=cumRiskScore/(duration/1800)) %>% 
  filter(avgRiskScore>100) %>%
  select(month,duration,maxRiskScore,avgRiskScore,bg_rate_cases_app,positive) %>%
  filter(month %in% monthshq) %>% 
  mutate(binduration=cut(duration/3600,
                         breaks = c(0,0.5,1,2,3,10,Inf),
                         labels=paste(c(0,0.5,1,2,3,">"),c(rep("-",5),""),c(0.5,1,2,3,10,10),sep="")),
         binmaxprox=cut((maxRiskScore),
                        breaks = c(100,120,140,160,200,300,Inf),
                        labels=paste(c(round(c(100,120,140,160,200)/90,1),">"),c(rep("-",5),""),round(c(120,140,160,200,300,300)/90,1),sep="")),
         binavgprox=cut((avgRiskScore),
                        breaks = c(100,120,140,160,200,300,Inf),
                        labels=paste(c(round(c(100,120,140,160,200)/90,1),">"),c(rep("-",5),""),round(c(120,140,160,200,300,300)/90,1),sep=""))) %>%
  group_by(binduration,binmaxprox) %>%
  summarise(n=n(),TPAEN=mean(positive-(1-(1-bg_rate_cases_app)^rho_ml))) %>%
  mutate(lowerCI=qbeta(0.025,round(TPAEN*n),n-round(TPAEN*n)+1),upperCI=qbeta(1-0.025,round(TPAEN*n)+1,n-round(TPAEN*n)))
results_avgriskscore_duration <- table_outcome_riskScores %>% mutate(avgRiskScore=cumRiskScore/(duration/1800)) %>%
  filter(avgRiskScore>100) %>%
  select(month,duration,maxRiskScore,avgRiskScore,bg_rate_cases_app,positive) %>%
  filter(month %in% monthshq) %>% 
  mutate(binduration=cut(duration/3600,
                         breaks = c(0,0.5,1,2,3,10,Inf),
                         labels=paste(c(0,0.5,1,2,3,">"),c(rep("-",5),""),c(0.5,1,2,3,10,10),sep="")),
         binmaxprox=cut((maxRiskScore),
                        breaks = c(100,120,140,160,200,300,Inf),
                        labels=paste(c(round(c(100,120,140,160,200)/90,1),">"),c(rep("-",5),""),round(c(120,140,160,200,300,300)/90,1),sep="")),
         binavgprox=cut((avgRiskScore),
                        breaks = c(100,120,140,160,200,300,Inf),
                        labels=paste(c(round(c(100,120,140,160,200)/90,1),">"),c(rep("-",5),""),round(c(120,140,160,200,300,300)/90,1),sep=""))) %>%
  group_by(binduration,binavgprox) %>%
  summarise(n=n(),TPAEN=mean(positive-(1-(1-bg_rate_cases_app)^rho_ml))) %>%
  mutate(lowerCI=qbeta(0.025,round(TPAEN*n),n-round(TPAEN*n)+1),upperCI=qbeta(1-0.025,round(TPAEN*n)+1,n-round(TPAEN*n)))
SaveTable(results_avgriskscore_duration)

miny<-min(c(results_avgriskscore_duration$lowerCI,results_maxriskscore_duration$lowerCI))/1.025
maxy<-max(c(results_avgriskscore_duration$upperCI,results_maxriskscore_duration$upperCI))*1.025
plot_duration_avgriskscore <- results_avgriskscore_duration %>% ggplot(aes(x=binduration,y=TPAEN,group=binavgprox)) + geom_line(aes(col=binavgprox)) + geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=binavgprox),alpha=0.3) + ylab("probability of\nreported transmission") + # ylab("probability of detected transmission") + 
  scale_y_continuous(limits=c(0,maxy)) + coord_cartesian(expand=FALSE) + scale_fill_viridis_d(direction=-1,option="magma") + scale_colour_viridis_d(direction=-1,option="magma") + theme_classic() + xlab("duration (hours)") + facet_grid(. ~ "by risk score per half-hour")
plot_duration_maxriskscore <- results_maxriskscore_duration %>% ggplot(aes(x=binduration,y=TPAEN,group=binmaxprox)) + geom_line(aes(col=binmaxprox)) + geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=binmaxprox),alpha=0.3) + ylab("probability of\nreported transmission") + # ylab("probability of detected transmission") + 
  scale_y_continuous(limits=c(0,maxy)) + coord_cartesian(expand=FALSE) + scale_fill_viridis_d(direction=-1,option="magma") + scale_colour_viridis_d(direction=-1,option="magma") + theme_classic() + xlab("duration (hours)") + facet_grid(. ~ "by maximum risk score")
plot_avgriskscore_duration <- results_avgriskscore_duration %>% ggplot(aes(x=binavgprox,y=TPAEN,group=binduration)) + geom_line(aes(col=binduration)) + geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=binduration),alpha=0.3) + ylab("probability of\nreported transmission") + # ylab("probability of detected transmission") + 
  scale_y_continuous(limits=c(0,maxy)) + coord_cartesian(expand=FALSE) + scale_fill_viridis_d(direction=-1) + scale_colour_viridis_d(direction=-1) + theme_classic() + xlab("risk score per half-hour") + facet_grid(. ~ "by duration")
plot_maxriskscore_duration <- results_maxriskscore_duration %>% ggplot(aes(x=binmaxprox,y=TPAEN,group=binduration)) + geom_line(aes(col=binduration)) + geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=binduration),alpha=0.3) + ylab("probability of\nreported transmission") + # ylab("probability of detected transmission") + 
  scale_y_continuous(limits=c(0,maxy)) + coord_cartesian(expand=FALSE) + scale_fill_viridis_d(direction=-1) + scale_colour_viridis_d(direction=-1) + theme_classic() + xlab("maximum risk score") + facet_grid(. ~ "by duration")
legend_riskscore <- cowplot::get_legend(
  # create some space to the left of the legend
  plot_duration_avgriskscore +
    guides(fill=guide_legend(title="risk score"),color=guide_legend(title="risk score")) + 
    theme(legend.box.margin = margin(0, 0, 0, 0))
)
legend_duration <- cowplot::get_legend(
  # create some space to the left of the legend
  plot_avgriskscore_duration + 
    guides(fill=guide_legend(title="duration\n (hours)"),color=guide_legend(title="duration\n (hours)")) + 
    theme(legend.box.margin = margin(0, 0, 0, 0))
)

plot_full_riskscore_duration <- plot_grid(plot_duration_avgriskscore + theme(legend.position="none",plot.margin =margin(t=5,r=15,b=5,l=5)) ,#+ scale_y_log10(limits=c(miny,maxy)),
                                          plot_duration_maxriskscore + #scale_y_log10(limits=c(miny,maxy)) + 
                                            theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5),
                                                  axis.title.y = element_blank()
                                            ),
                                          legend_riskscore,
                                          plot_avgriskscore_duration + theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5)) ,#+ scale_y_log10(limits=c(miny,maxy)),
                                          plot_maxriskscore_duration + #scale_y_log10(limits=c(miny,maxy)) +
                                            theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5),
                                                  axis.title.y = element_blank()
                                            ),
                                          legend_duration,
                                          nrow=2,align="hv",axis="lb",rel_widths=c(1,1,0.4))
SaveFigure(plot_full_riskscore_duration,width=20,height=16)

results_avgriskscoreperhour_duration <- table_outcome_riskScores %>% mutate(avgRiskScore=cumRiskScore/(duration/3600)) %>%
  filter(avgRiskScore>200) %>%
  select(month,duration,maxRiskScore,avgRiskScore,bg_rate_cases_app,positive) %>%
  filter(month %in% monthshq) %>% 
  mutate(binduration=cut(duration/3600,
                         breaks = c(0,0.5,1,2,3,10,Inf),
                         labels=paste(c(0,0.5,1,2,3,">"),c(rep("-",5),""),c(0.5,1,2,3,10,10),sep="")),
         binmaxprox=cut((maxRiskScore),
                        breaks = c(100,120,140,160,200,300,Inf),
                        labels=paste(c(round(c(100,120,140,160,200)/90,1),">"),c(rep("-",5),""),round(c(120,140,160,200,300,300)/90,1),sep="")),
         binavgprox=cut((avgRiskScore),
                        breaks = c(100,120,140,160,200,300,Inf)*2,
                        labels=paste(c(round(2*c(100,120,140,160,200)/90,1),">"),c(rep("-",5),""),round(2*c(120,140,160,200,300,300)/90,1),sep=""))) %>%
  group_by(binduration,binavgprox) %>%
  summarise(n=n(),TPAEN=mean(positive-(1-(1-bg_rate_cases_app)^rho_ml)),meanduration=mean(duration/3600),meanavgprox=mean(avgRiskScore)) %>%
  mutate(lowerCI=qbeta(0.025,round(TPAEN*n),n-round(TPAEN*n)+1),upperCI=qbeta(1-0.025,round(TPAEN*n)+1,n-round(TPAEN*n)))

plot_main_riskscore_duration <- results_avgriskscoreperhour_duration %>% 
  filter(!is.na(binavgprox)) %>% 
  ggplot(aes(x=meanduration,y=TPAEN,group=binavgprox)) + 
    geom_line(aes(col=binavgprox)) + 
    geom_ribbon(aes(ymin=lowerCI,ymax=upperCI,fill=binavgprox),alpha=0.3) + 
    ylab("probability of reported transmission") + # ylab("probability of detected transmission") + 
    scale_x_log10() + 
    scale_y_log10() + 
    coord_cartesian(expand=FALSE) + 
    scale_fill_viridis_d(direction=-1,option="magma") + 
    scale_colour_viridis_d(direction=-1,option="magma") + theme_classic() + xlab("duration (hours)")  + 
    guides(fill=guide_legend(title="risk score\n per hour"),color=guide_legend(title="risk score\n per hour"))
SaveFigure(plot_main_riskscore_duration,width=12,height=10)


plot_full_riskscore_duration_log <- plot_grid(plot_duration_avgriskscore + theme(legend.position="none",plot.margin =margin(t=5,r=15,b=5,l=5)) + scale_y_log10(limits=c(miny,maxy)),
                                              plot_duration_maxriskscore + scale_y_log10(limits=c(miny,maxy)) + 
                                                theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5),
                                                      axis.title.y = element_blank()
                                                ),
                                              legend_riskscore,
                                              plot_avgriskscore_duration + theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5)) + scale_y_log10(limits=c(miny,maxy)),
                                              plot_maxriskscore_duration + scale_y_log10(limits=c(miny,maxy)) +
                                                theme(legend.position="none",plot.margin=margin(t=5,r=15,b=5,l=5),
                                                      axis.title.y = element_blank()
                                                ),
                                              legend_duration,
                                              nrow=2,align="hv",axis="lb",rel_widths=c(1,1,0.4))
SaveFigure(plot_full_riskscore_duration_log,width=20,height=16)

# joint plot on TPAEN vs cumrisk,duration with regressions
ReadTable("table_overall_joint")
ReadTable("table_overall_joint_nobg")
table_overall_joint_regr <- bind_rows(
  table_overall_joint %>%
    mutate(type="...infection",name=factor(name,ordered = TRUE,levels = c("maximum risk score","cumulative risk score","duration (hours)"))),
  table_overall_joint_nobg %>%
    mutate(type="...transmission",name=factor(name,ordered = TRUE,levels = c("maximum risk score","cumulative risk score","duration (hours)")))
) %>%
  mutate(type=factor(type,levels=c("...infection","...transmission")))
plot_TPAEN_loglogregression <- ggplot(data=table_overall_joint_regr  %>%
                                        filter(name!="maximum risk score") %>%
                                        mutate(meanscore=if_else(name!="duration (hours)",meanscore/90,meanscore)),
                                      aes(x=meanscore,y=TPAEN)) + 
  geom_errorbar(aes(ymin=lowerCI,ymax=upperCI),col="blue") + 
  geom_point(col="blue") + 
  xlab("score") +
  ylab("probability of reported...") +
  theme_classic() +
  scale_y_continuous(expand=expansion(add=0),trans = "log10") +
  scale_x_continuous(expand=expansion(mult=0.05),trans = "log10") +
  coord_cartesian(expand=TRUE) + facet_wrap(~name,scale="free_x") +
  geom_smooth(method="lm",col=alpha("brown",0.3),fill=alpha("brown",0.1)) +
  stat_regline_equation(
    aes(label =  ..eq.label..),label.y.npc = 1.05,label.x.npc = 0,
    col="brown") +
  geom_smooth(data = table_overall_joint_regr  %>%
                mutate(meanscore=if_else(name!="duration (hours)",meanscore/90,meanscore)) %>%
                filter(type=="...transmission") %>%
                filter((name=="duration (hours)" & meanscore<3) | (name=="cumulative risk score" & meanscore<20) )
              ,method="lm",col=alpha("orange",0.5),fill=alpha("orange",0.3)) +
  stat_regline_equation(
    data = table_overall_joint_regr  %>%
      mutate(meanscore=if_else(name!="duration (hours)",meanscore/90,meanscore)) %>%
      filter(type=="...transmission") %>%
      filter((name=="duration (hours)" & meanscore<3) | (name=="cumulative risk score" & meanscore<20) ),
    aes(label =  ..eq.label..),label.y.npc = 0.05,label.x.npc = 0.3,
    col="orange"
    ) +
  facet_grid(~ type + name)
SaveFigure(plot_TPAEN_loglogregression,width = 20,height = 8)
# regression coefficients
table_overall_joint %>% filter(name=="cumulative risk score") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.) %>% .$coefficients
table_overall_joint %>% filter(name=="duration (hours)") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.) %>% .$coefficients
table_overall_joint_nobg %>% filter(name=="cumulative risk score") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.) %>% .$coefficients
table_overall_joint_nobg %>% filter(name=="duration (hours)") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.) %>% .$coefficients
table_overall_joint_nobg %>% filter(meanscore/90<20,name=="cumulative risk score") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.) %>% .$coefficients
table_overall_joint_nobg %>% filter(meanscore<3,name=="duration (hours)") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.) %>% .$coefficients
#
fit_cum<-table_overall_joint %>% filter(name=="cumulative risk score") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.)
summary(fit_cum)
confint(fit_cum,"I(log(meanscore))",0.95)
fit_dur<-table_overall_joint %>% filter(name=="duration (hours)") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.) 
summary(fit_dur)
confint(fit_dur,"I(log(meanscore))",0.95)
fit_cum_tr<-table_overall_joint_nobg %>% filter(name=="cumulative risk score") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.) 
summary(fit_cum_tr)
confint(fit_cum_tr,"I(log(meanscore))",0.95)
fit_dur_tr<-table_overall_joint_nobg %>% filter(name=="duration (hours)") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.) 
summary(fit_dur_tr)
confint(fit_dur_tr,"I(log(meanscore))",0.95)
fit_cum_low<-table_overall_joint_nobg %>% filter(meanscore/90<20,name=="cumulative risk score") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.) 
summary(fit_cum_low)
confint(fit_cum_low,"I(log(meanscore))",0.95)
fit_dur_low<-table_overall_joint_nobg %>% filter(meanscore<3,name=="duration (hours)") %>% lm(I(log(TPAEN)) ~ I(log(meanscore)),data=.) 
summary(fit_dur_low)
confint(fit_dur_low,"I(log(meanscore))",0.95)

