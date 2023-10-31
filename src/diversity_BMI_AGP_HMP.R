#################### BMI and Diversity for HMP and AGP study ########################
pkg <- c("tidyverse","coin", "cowplot","ggthemes")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}

meta <- read.csv("meta_BMI_obesity.csv")
metric <- read.csv("combine_metrics_rdp_lg2fc.txt", sep='\t')
metricTotal <- metric %>%
  mutate(Diversity=if_else(Metric=="FBratio", log2(Diversity), Diversity))

submetawestern <- meta%>%filter(!is.na(BMI))%>%filter(BMI<80)

submetawestern <- submetawestern%>%mutate(BMI_interval=case_when(
  BMI<16.5 ~ 16,
  BMI>42 ~ 42,
  TRUE ~ round(BMI,digits=0)
))
submetawestern <- submetawestern%>%filter(StudyID %in% c("AGP","HMP"))
submetawestern <- submetawestern[,c("Run","BMI_interval","race")]

tablewestern <- plyr::join(metricTotal, submetawestern, by="Run",type="inner") %>% 
  filter(race !="Other") %>%
  filter(race !="Hispanic") %>%
  group_by(StudyID, BMI_interval, Metric, race, group)%>%
  dplyr::summarize(mean=mean(Diversity), sd=sd(Diversity)) %>%
  ungroup()

BMI_western <- tibble()
for (j in unique(tablewestern$race)) {
  for (i in unique(tablewestern$Metric)){
    datalmm <- subset(tablewestern,Metric==i) %>% filter(race==j)
    fit<-lmerTest:::lmer(mean~BMI_interval+(1|StudyID), data=datalmm)
    beta <- summary(fit)$coefficients["BMI_interval", "Estimate"]
    Pval <- summary(fit)$coefficients["BMI_interval", "Pr(>|t|)"]
    OR <- exp(summary(fit)$coefficients["BMI_interval", "Estimate"])
    cf<-lme4::confint.merMod(object = fit,method ="Wald")
    ORCI <- exp(cf)
    lowOR <- ORCI["BMI_interval",1]
    highOR <- ORCI["BMI_interval",2]
    low <- cf["BMI_interval",1]
    high <- cf["BMI_interval",2]
    so <- summary(fit)
    corr <- so$vcov@factors$correlation["BMI_interval",1]
    BMI_western <- bind_rows(BMI_western, tibble(
      Study=j,Metric=i,OR=OR,low_OR=lowOR,high_OR=highOR,beta=beta,low_CI=low, high_CI=high, Pval= Pval,correlation=corr
    ))
  }
}

BMI_western <- mutate(BMI_western, Significance=if_else(Pval<0.05, "sig","not sig")) 
BMI_western[BMI_western$Metric=="FBratio", ]$Metric <- "log(F/B ratio)"
BMI_western$Metric <- factor(BMI_western$Metric, levels = c("Shannon","Evenness","OBS","Chao1","log(F/B ratio)"))

forestci <- ggplot(aes(x=OR, y=Study,color=Study), data=BMI_western) + 
  geom_vline(xintercept =1, linetype="dashed", color="grey50") +
  geom_linerange(aes(xmin=low_OR, xmax=high_OR), colour='darkgrey') + 
  geom_point(aes(shape=Significance)) +
  facet_grid(~Metric, scales="free") +
  theme_few() + 
  scale_shape_manual(values = c(1,16))+
  scale_color_manual(values=c('#66c2a5','#fc8d62','#8da0cb'))+
  theme(axis.text.x=element_text(angle=45, hjust=1,size=9), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank()
  )

tablewestern$low <- tablewestern$mean - tablewestern$sd
tablewestern$high <- tablewestern$mean + tablewestern$sd
tablewestern[is.na(tablewestern$low),]$low <- 0
tablewestern[is.na(tablewestern$high),]$high <- 0
tablewestern[tablewestern$Metric=="FBratio", ]$Metric <- "log(F/B ratio)"

tablewestern$Metric <- factor(tablewestern$Metric, levels = c("Shannon","Evenness","OBS","Chao1","log(F/B ratio)"))

regression <- ggplot(tablewestern,aes(x=BMI_interval, y=mean, ymin=mean-sd, ymax=mean+sd, color=race)) +
  geom_smooth(se=F,method = "lm") +
  geom_errorbar(width=0) +
  geom_point(aes(shape = race)) +
  facet_wrap(~Metric, scales="free",nrow = 1) +
  theme_few() +
  theme(axis.text.y = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1,size=9),
        axis.ticks.y=element_blank())+
  xlab("BMI") +
  ylab("Diversity¡ÀSD") +
  scale_color_manual(values=c('#66c2a5','#fc8d62','#8da0cb'))

forest_reg <- plot_grid(regression, forestci,  align = 'h', ncol = 1, rel_heights = c(1,0.3))

save_plot("F:/Rproject/pic/meta_fig/Forest_HMP_AGP_diversity.pdf", forest_reg, base_height=6)


