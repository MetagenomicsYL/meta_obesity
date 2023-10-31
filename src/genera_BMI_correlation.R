##########################association of genera and BMI #############################
pkg <- c("vegan", "tidyverse", "cowplot", "coin", "boot", "ggsci", "ggpubr", "ggthemes")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}
source('F:/Rproject/pic/meta_code/function_library.R')
setwd("F:/Rproject/pic/meta_dataset")

meta <- read.csv("meta_BMI_obesity.csv")
tablegenus <- read.csv("genus_combined_rel_rdp.csv")
tablegenus$Genus <- gsub("g__","",tablegenus$Genus)
tablegenus <- tablegenus%>%column_to_rownames("Genus")
submeta <- meta%>%filter(!is.na(BMI))%>%filter(BMI<100)
tablegenus <- tablegenus[,submeta$Run]

markerdown <- read.csv("REM_FC_P_count12.csv")%>% 
  filter(REM<0.05) %>% 
  filter(FDR<0.001) %>% 
  filter(FC<0)

markerup <- read.csv("REM_FC_P_count12.csv")%>% 
  filter(REM<0.05) %>% 
  filter(FDR<0.001) %>% 
  filter(FC>0)

###################### negative association of BMI and genera ################### 
deplete <- markerdown$genuslevel
tabledown <- tablegenus[deplete,]

submetadown <- submeta[,c("StudyID","Run","BMI","group","race")]

tabledownmeta <- tabledown%>%t()%>%as.data.frame()%>%
  rownames_to_column("Run")
tabledownmeta <- plyr::join(tabledownmeta,submetadown,by="Run",type="inner")

BMIdowntotal <- REMfun(tabledownmeta,"Combined")

tabledown_white <- tabledownmeta %>% filter(race=="Caucasian")
BMIresultwhite <- REMfun(tabledown_white,"Caucasian")

tabledown_asian <- tabledownmeta %>% filter(race=="Asian")
BMIresultasian <- REMfun(tabledown_asian,"Asian")

tabledown_african <- tabledownmeta %>% filter(race=="African")
BMIresultafrican <- REMfun(tabledown_african,"African")

totaldown <- rbind(BMIdowntotal,BMIresultwhite,BMIresultasian,BMIresultafrican)

totaldown$logP <--log10(totaldown$Pval)

totaldown <- mutate(totaldown,significance=case_when(
  Pval<0.05 ~ "true",
  TRUE ~"false"
))

genusorder <-subset(totaldown,Study=="Asian")%>%arrange(logP) %>% pull(Genus)

totaldown$Genus <- factor(totaldown$Genus,levels = genusorder)
totaldown$REM_beta <- abs(totaldown$beta)
totaldown <- mutate(totaldown,Effect=case_when(
  beta>0 ~ "positve",
  TRUE ~ "negative"
))

write.table(totaldown,file="bubble_depleted_genera.csv",quote=F, sep=',',row.names=F)

totaldown.plot <- totaldown %>% filter(beta<0)
bubble_deplete <-  ggplot(totaldown.plot,aes(logP,Genus,color=Study))+ 
  geom_point(aes(size=-beta,shape=significance))+
  scale_color_manual(values=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3')) + 
  scale_shape_manual(values=c(21,19)) +
  geom_vline(xintercept = totaldown.plot %>% 
               filter(Pval <0.05) %>% 
               pull(logP) %>% min, col='#525252',linetype="dashed") +
  theme_few()+
  xlab('REM -log10(Pvalue)')

ggsave("F:/Rproject/pic/meta_fig/BMI_negative_beta_bubble.pdf", bubble_deplete, width = 9, height =5)


###################### positive association of BMI and genera ################### 

enrich <- markerup$genuslevel
tableup <- tablegenus[enrich,]
tableup <- tableup%>%t()%>%as.data.frame()%>%
  rownames_to_column("Run")

tableup_meta <- plyr::join(tableup,submetadown,by="Run",type="inner")
BMIuptotal <- REMfun(tableup_meta,"Combined")

tableup_white <- tableup_meta %>% filter(race=="Caucasian")
BMIupwhite <- REMfun(tableup_white,"Caucasian")

tableup_asian <- tableup_meta %>% filter(race=="Asian")
BMIupasian <- REMfun(tableup_asian,"Asian")

tableup_african <- tableup_meta %>% filter(race=="African")
BMIupafrican <- REMfun(tableup_african,"African")

totalup <- rbind(BMIuptotal,BMIupwhite,BMIupasian,BMIupafrican)
totalup$logP <--log10(totalup$Pval)

totalup <- mutate(totalup,significance=case_when(
  Pval<0.05 ~ "true",
  TRUE ~"false"
))

genusorder <- subset(totalup, Study=="Asian") %>% arrange(logP) %>% pull(Genus)

totalup$Genus <- factor(totalup$Genus,levels = genusorder)
totalup$REM_beta <- abs(totalup$beta)
totalup <- mutate(totalup,Effect=case_when(
  beta>0 ~ "positve",
  TRUE ~ "negative"
))

write.table(totalup,file="bubble_enriched_genera.csv",quote=F, sep=',',row.names=F)

totalup.plot <- totalup %>% filter (beta > 0)

bubble_enrich <-  ggplot(totalup.plot,aes(logP,Genus,color=Study))+ 
  geom_point(aes(size=beta,shape=significance))+
  scale_color_manual(values=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3')) + 
  scale_shape_manual(values=c(21,19)) +
  geom_vline(xintercept = totalup.plot %>% 
               filter(Pval <0.05) %>% 
               pull(logP) %>% min, col='#525252',linetype="dashed") +
  theme_few()+
  xlab('REM -log10(Pvalue)')

ggsave("F:/Rproject/pic/meta_fig/BMI_positive_beta_bubble.pdf", bubble_enrich, width = 9, height =5)

