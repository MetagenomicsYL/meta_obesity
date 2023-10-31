######################### calculate for the combined alpha diversity ######################
### 
# some codes referred https://jbisanz.github.io/MetaDiet/analysis/HFD_diversity.html
###

pkg <- c("tidyverse", "Matrix", "lme4", "ggthemes")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}

setwd("F:/Rproject/pic/meta_dataset")
metric <- read.csv("combined_metrics_rdp.txt",sep = "\t")

log2r <- list()
for (y in unique(metric$StudyID)){
  metricdata <- subset(metric, StudyID==y)
  eachmetric <- rename(metricdata, Chao1=chao1,Shannon = shannon_entropy, OBS=observed_features,FBratio=FB,Evenness=pielou_evenness) %>%
    na.omit()
  metric2 <- eachmetric[,-which(colnames(eachmetric) %in% c("b","f"))]%>%
    select("StudyID","Run","group",everything())
  metric_gather <- gather(metric2, Metric, Diversity, OBS:FBratio)
  metric_smr <- metric_gather %>%
    group_by(Metric, group) %>% 
    summarize(mean=mean(log2(Diversity))) %>%
    spread(key=group, value=mean) %>%
    rename(mean_log2_Healthy=Healthy, mean_log2_Obese=Obese) 
  metric_join <- left_join(metric_gather, metric_smr)
  metric_join$log2FC <- log2(metric_join$Diversity)-(metric_join$mean_log2_Healthy)
  log2r[[y]] <- metric_join
}

metricTotal <- log2r%>%
  do.call(bind_rows, .)
write.table(metricTotal, file="combine_metrics_rdp_lg2fc.txt", quote=F, sep='\t', row.names=F)

singledata<-list()
for (y in unique(metric$StudyID)){
  eachst <- list()
  eachst$AlphaDiversity <- log2r[[y]]
  eachst$AlphaDiversity_Stats <-
    eachst$AlphaDiversity %>%
    group_by(Metric) %>%  
    do(
      broom::tidy(t.test(log2FC~group, data=., conf.int=TRUE, conf.level=0.95)) 
    ) %>%
    mutate(StudyID=y) %>%
    select(StudyID, Metric, log2FC=estimate, Pvalue=p.value, mean_Healthy=estimate1, mean_Obese=estimate2, CI_low=conf.high, CI_high=conf.low)
  eachst$AlphaDiversity_Stats$log2FC <- eachst$AlphaDiversity_Stats$log2FC
  eachst$AlphaDiversity_Stats$CI_low <- eachst$AlphaDiversity_Stats$CI_low
  eachst$AlphaDiversity_Stats$CI_high <- eachst$AlphaDiversity_Stats$CI_high
  singledata[[y]] <- eachst
}


Alphasum<-lapply(singledata, function(x) x$AlphaDiversity)
Alphasum <- Alphasum%>%
  do.call(bind_rows, .) %>%
  mutate(group=factor(group, levels=c("Healthy","Obese")))

AlphaCombined<-tibble()
for(i in unique(Alphasum$Metric)){
  fit<-lmerTest:::lmer(log2FC~group+(1|StudyID), data=subset(Alphasum, Metric==i))
  cf<-confint(fit,level = 0.95)
  temp <- tibble(
    StudyID="Combined", 
    Metric=i, 
    log2FC=summary(fit)$coefficients["groupObese", "Estimate"], 
    Pvalue=anova(fit)$`Pr(>F)`, 
    mean_Healthy=NA, 
    mean_Obese=NA, 
    CI_low=cf["groupObese",1], 
    CI_high=cf["groupObese",2]
  )
  AlphaCombined<-bind_rows(AlphaCombined, temp)
}

Aisa <- c("He", "HKGutMicMap", "Gao", "Zeevi", "Ahmad")
Europe <- c("Goodrich", "Olsson", "Lippert")
America <- c("AGP", "HMP", "Zupancic", "Stanislawski", "Wu", "Barengolts ", "Ross", " Baxter", "Chavez", "Escobar")

SampleSort <- lapply(names(singledata), function(x) tibble(StudyID=x, Nsamples=length(unique(singledata[[x]]$AlphaDiversity$Run)))) %>% 
  do.call(bind_rows, .)%>% 
  mutate(Geography=case_when(
    StudyID %in% Aisa ~"Asia",
    StudyID %in% Europe ~"Europe",
    TRUE ~ "America"
  ))

SampleSort <- SampleSort[order(SampleSort$Geography, SampleSort$Nsamples,decreasing = T),]%>% 
  mutate(Study=StudyID)%>%
  bind_rows(tibble(StudyID="Combined", Study="Combined",Geography="Combined")) %>%
  mutate(Study=factor(Study, levels=rev(Study))) %>%
  mutate(StudyID=factor(StudyID, levels=rev(StudyID))) %>%
  mutate(Geography=factor(Geography, levels=c("America","Asia","Europe","Combined")))

alphatotal <- lapply(singledata, function(x) x$AlphaDiversity_Stats)%>% 
  do.call(bind_rows, .) %>% 
  bind_rows(AlphaCombined) %>%
  mutate(Significance=case_when(Pvalue<0.05 & log2FC>0 ~ "up", 
                                Pvalue<0.05 & log2FC<0 ~ "down",
                                TRUE~"not sig"))%>%
  ungroup()%>%
  left_join(SampleSort)
alpha_fig <- alphatotal %>%
  mutate(type=case_when(
    StudyID!="Combined" ~ "one study",
    TRUE~"total effect"
  ))
alpha_fig <- alpha_fig%>%
  mutate(Metric=case_when(
    Metric=="FBratio"~"F/B ratio",
    TRUE~as.character(Metric)))%>%
  mutate(Metric=factor(Metric, levels=c("Shannon", "Evenness","OBS","Chao1","F/B ratio")))

alpha_plot <- ggplot(alpha_fig, aes(x=log2FC, y=Study, color=Geography)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_linerange(aes(xmin=CI_low, xmax=CI_high), colour='darkgrey') + 
  geom_point(aes(shape = Significance)) +
  facet_grid(~Metric, scales="free_x") +
  theme_few() +
  scale_color_manual(values=rev(c('#034B61','#E3A6A1','#19B3B1','#BC5F6A')))+
  scale_shape_manual(values=c(17,1,16))+
  theme(axis.text.x=element_text(angle=45, hjust=1,size=9), 
        axis.text.y = element_text(size=11),
        legend.title = element_text(size=11),
        axis.title.x=element_text(size=12,face = "bold"))+
  ylab("")

ggsave("F:/Rproject/pic/meta_fig/combined_alpha_geo.pdf", alpha_plot, width = 9, height = 4.5)

