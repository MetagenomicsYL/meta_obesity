###############calculated stratified wilcoxon test#################################
### 
# some codes referred https://github.com/zellerlab/crc_meta
###
pkg <- c("tidyverse", "coin")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}
source('F:/Rproject/pic/meta_code/function_library.R')
setwd("F:/Rproject/pic/meta_dataset")

######################### prepare for the genus level OTU table ##################
OTUtable <- read.csv("combined_uniqueOTU_rdp.csv") %>%
  column_to_rownames("Taxon")
OTUtableflt <- Confidence.Filter(OTUtable, 3, 10, TRUE) %>% 
  as.data.frame() %>%
  rownames_to_column("Taxon")

write.table(OTUtableflt, file="combined_uniqueOTU_filter_ab_rdp.csv", quote=FALSE, sep=',',
            row.names=FALSE, col.names=TRUE)

OTUsep <- OTUtableflt %>% separate(Taxon, 
                                    c("Kingdom",
                                      "Phylum",
                                      "Class",
                                      "Order",
                                      "Family",
                                      "Genus",
                                      "Species"
                                    ), 
                                    sep="; ", 
                                    remove=FALSE)
OTUgenus <- dplyr::select(OTUsep,-c("Taxon","Kingdom", "Phylum", "Class", 
                                    "Order", "Family", "Species")) %>%
  na.omit() %>%
  filter(Genus!="g__")
Genusagg <- aggregate(OTUgenus[,-1], by = list(Genus=OTUgenus$Genus), FUN = sum)

write.table(Genusagg, file="genus_combinedOTU_rdp.csv", quote=F, sep=',', row.names=F)


filter_name <- OTUtableflt %>% column_to_rownames("Taxon")
OTUrel <- prop.table(as.matrix(OTUtable), 2) %>% as.data.frame()
OTUrel <- OTUrel[row.names(filter_name),] %>% rownames_to_column("Taxon")
OTUrelsep <- OTUrel %>% separate(Taxon, 
                                  c("Kingdom",
                                    "Phylum",
                                    "Class",
                                    "Order",
                                    "Family",
                                    "Genus",
                                    "Species"
                                  ), 
                                  sep="; ", 
                                  remove=FALSE)
tem <- dplyr::select(OTUrelsep,-c("Taxon","Kingdom", "Phylum", "Class", "Order", "Family", "Species"))
OTUrelgenus <- subset(tem,!is.na(Genus))%>%
  na.omit() %>%
  filter(Genus!="g__")

Genusrelagg <- aggregate(OTUrelgenus[,-1], by = list(Genus=OTUrelgenus$Genus), FUN = sum)
write.table(Genusrelagg,file="genus_combined_rel_rdp.csv",quote=F, sep=',',row.names=F)

#################### calculate Stratified Wilcoxon test ##########################

meta <- read.csv("meta_BMI_obesity.csv")
studies <- meta %>% pull(StudyID) %>% unique
Genusrelagg <- Genusrelagg %>% column_to_rownames("Genus")
p_value <- matrix(NA, nrow=nrow(Genusrelagg), ncol=length(studies)+1, 
                dimnames=list(row.names(Genusrelagg), c(studies, 'all')))
bcnames <- Genusagg$Genus

for (i in studies){
  for (t in bcnames){
    metasty <- subset(meta,StudyID==i)
    otusty <- Genusrelagg[,metasty$Run]
    totusty <- t(otusty)%>%as.data.frame()%>%rownames_to_column("Run")
    metagroup <- metasty[,c("Run","group")]
    merge <- totusty%>%left_join(metagroup)
    w <- wilcox.test(merge[,t][merge[,"group"]=="Obese"], merge[,t][merge[,"group"]=="Healthy"])
    p_value[t,i] <- w$p.value
  }
}

for (t in bcnames){
  tfeatall=t(Genusrelagg)%>%as.data.frame()%>%rownames_to_column("Run")
  metagroup <- meta[,c("Run","group","StudyID")]
  metagroup$StudyID <- factor(metagroup$StudyID)
  metagroup$group<- factor(metagroup$group)
  merge <- tfeatall%>%left_join(metagroup)
  d <- data.frame(y=merge[,t],x=merge$group, block=merge$StudyID)
  p_value[t,"all"]  <-pvalue(wilcox_test(y ~ x | block, data=d))
}


fdrvalue <- data.frame(apply(p_value, MARGIN=2, FUN=p.adjust, method='fdr'),
                    check.names = FALSE)
fdrvalue <- fdrvalue%>%rownames_to_column("Genus")
write.table(fdrvalue, file = "p_adj.csv", 
            quote=FALSE, sep=',',row.names = F)
p_value<- p_value%>%as.data.frame()%>%rownames_to_column("Genus")
write.table(p_value, file = "p_val.csv", 
            quote=FALSE, sep=',',row.names = F)


         
