######################select significant microbial genera##############################
pkg <- c("qiime2R","tidyverse","broom","Matrix","lme4")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}
source('F:/Rproject/pic/meta_code/function_library.R')

setwd("F:/Rproject/pic/meta_dataset")

genusOTU <- read.csv("genus_combined_rel_rdp.csv") %>%
  mutate(genusOTU, Genus=gsub("g__","", Genus)) %>%
  column_to_rownames("Genus")

metagroup <- read.csv("meta_BMI_obesity.csv") %>%
  dplyr::rename(Samplename=SampleID, SampleID=Run)

genusOTU <- genusOTU[,metagroup$SampleID]

#############################for all marker_candi_OTU#################################
OTU_czm<-zCompositions::cmultRepl(t(genusOTU),
                                 label=0,
                                 method="CZM",
                                 output="prop") %>%
                                 as.data.frame() %>%
                                 rownames_to_column("SampleID")
upsetOTU <- read.csv("genus_count12_rdp.csv") %>%
  mutate(Genus=gsub("g__","", Genus))
upsetgenus <- upsetOTU$Genus


selectOTU <- OTU_czm[,c("SampleID",upsetgenus)]%>%as.data.frame()%>%
  left_join(metagroup)

metacol <- colnames(metagroup)

OTUgather <- gather(selectOTU, key=genuslevel, value=reads, -all_of(metacol))

OTUlog <- cal_logfun(metagroup, OTUgather, "genuslevel")

singledata <- cal_log2fc(metagroup, OTUlog, "genuslevel")

OTUlogcom <- combine_log(singledata)

otuCombined <- combine_logfc(OTUlogcom, "genuslevel")

allstdata <- lapply(singledata, function(x) x$logstats)%>% 
  do.call(bind_rows, .) %>% 
  bind_rows(otuCombined) %>%
  mutate(Significance=case_when(Pvalue<0.05 & log2FC>0 ~ "up", 
                                Pvalue<0.05 & log2FC<0 ~ "down",
                                TRUE~"not sig"))%>%
  ungroup()%>%
  mutate(genuslevel=factor(genuslevel, levels=c(upsetgenus)))

write.table(allstdata,file="OTU_forest_plot_count12.csv",quote=F, sep=',',row.names=F)

############# save random effect model FC value ##################################
coefexpo<- allstdata[,c("StudyID","genuslevel","log2FC")]%>%
  spread(key=StudyID,value=log2FC)
studyname <- unique(metagroup$StudyID)
coefexpo <- coefexpo[,c("genuslevel",studyname,"Combined")]
write.table(coefexpo,file="REM_FC_rdp_study_count12.csv",quote=F, sep=',',row.names=F)

############# save random effect model P value ###################################

pexpo<- allstdata[,c("StudyID","genuslevel","Pvalue")]%>%
  spread(key=StudyID,value=Pvalue)
pexpo <- pexpo[,c("genuslevel",studyname,"Combined")]
write.table(pexpo,file="REM_P_rdp_study_count12.csv",quote=F, sep=',',row.names=F)

###############save each OTU count in each study and wilcoxon test result###########
pcount <- read.csv("genus_count12_rdp.csv")%>%
  rename(genuslevel=Genus)
pvalotu <-read.csv("p_val.csv") %>%
  rename(genuslevel=Genus)
pvalotu <- pvalotu[,c("genuslevel","all")]
pcount <- plyr::join(pcount,pvalotu, by="genuslevel",type="inner")

pcount$genuslevel <- gsub("g__","",pcount$genuslevel)
write.table(pcount,file="marker_wilpval_count12.csv",quote=F, sep=',',row.names=F)
