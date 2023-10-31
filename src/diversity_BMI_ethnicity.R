#################### BMI and Diversity stratified by ethnicity#######################
pkg <- c("tidyverse","coin", "cowplot","ggthemes")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}

meta <- read.csv("meta_BMI_obesity.csv")
metric <- read.csv("combine_metrics_rdp_lg2fc.txt", sep='\t')

submeta <- meta %>% filter(!is.na(BMI)) %>% filter(BMI<100)
submeta <- submeta %>% mutate(BMI_interval=round(BMI,digits=0))
submeta <- submeta[, c("Run","race","BMI_interval","BMI")]

metricTotal <- metric %>%
  mutate(Diversity=if_else(Metric=="FBratio", log2(Diversity), Diversity))
tablemeta <- plyr::join(metricTotal, submeta, by="Run", type="inner")
tablemetasub <- tablemeta[, c("StudyID","Run","group","Metric","Diversity","race","BMI","BMI_interval")]

corr_BMI_div <- function(dftable, studytype){
  Int_BMItotal <- tibble()
  for (i in unique(dftable$Metric)){
    datalmm <- subset(dftable, Metric==i)
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
    Int_BMItotal <- bind_rows(Int_BMItotal, tibble(
      Study=studytype, Metric=i, OR=OR, low_OR=lowOR, high_OR=highOR, beta=beta, low_CI=low, high_CI=high, Pval= Pval, correlation=corr
    ))
  }
  return(Int_BMItotal)
}

tabletotal <- tablemetasub%>%
  group_by(StudyID, BMI_interval,Metric,group)%>%
  dplyr::summarize(mean=mean(Diversity), sd=sd(Diversity)) %>%
  ungroup()

BMItotal <- corr_BMI_div(tabletotal, "Combined")

tablecacusian <- tablemetasub %>% 
  filter(race=="Caucasian") %>%
  group_by(StudyID, BMI_interval,Metric,group)%>%
  dplyr::summarize(mean=mean(Diversity), sd=sd(Diversity)) %>%
  ungroup()
BMIwhite <- corr_BMI_div(tablecacusian, "Caucasian")

tableafrican <- tablemetasub %>% 
  filter(race=="African") %>%
  group_by(StudyID, BMI_interval,Metric,group)%>%
  dplyr::summarize(mean=mean(Diversity), sd=sd(Diversity)) %>%
  ungroup()
BMIafrican <- corr_BMI_div(tableafrican, "African")

tableaasian <- tablemetasub %>% 
  filter(race=="Asian") %>%
  group_by(StudyID, BMI_interval,Metric,group) %>%
  dplyr::summarize(mean=mean(Diversity), sd=sd(Diversity)) %>%
  ungroup()
BMIasian <- corr_BMI_div(tableaasian, "Asian")

bmitable <- rbind(BMItotal, BMIasian, BMIwhite, BMIafrican)

# write.table(bmitable,file="race_diversity_BMIint.csv",quote=F, sep=',',row.names=F)

bmitable <- mutate(bmitable, Significance=if_else(Pval<0.05, "sig","not sig")) 
bmitable[bmitable$Metric=="FBratio", ]$Metric <- "log(F/B ratio)"
bmitable$Metric <- factor(bmitable$Metric, levels = c("Shannon", "Evenness", "OBS", "Chao1", "log(F/B ratio)"))

div_BMI <- ggplot(aes(x=OR, y=Study,color=Study), data=bmitable) + 
  geom_vline(xintercept =1, linetype="dashed", color="grey50") +
  geom_linerange(aes(xmin=low_OR, xmax=high_OR), colour='darkgrey') + 
  geom_point(aes(shape=Significance)) +
  facet_wrap(~Metric, scales="free") +
  theme_few() + 
  scale_shape_manual(values = c(1, 16))+
  scale_color_manual(values=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3'))+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9), 
        legend.position = "left",
        legend.title = element_blank(),
  )

ggsave("F:/Rproject/pic/meta_fig/diversity_forest_plot.pdf", div_BMI, width = 9, height =6)
