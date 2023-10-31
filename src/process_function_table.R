############# Stratified Wilcoxon test and REM for functional table ####################
pkg <- c("tidyverse", "coin")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}
setwd("F:/Rproject/pic/meta_dataset")
meta <- read.csv("meta_BMI_obesity.csv")
pathwaytable <- read.csv("combined_pathway.txt",sep = "\t") %>%
  column_to_rownames("pathway")
pathwaytable <- pathwaytable[, meta$Run]

p_value <- matrix(NA, nrow=nrow(pathwaytable), ncol=length(unique(meta$StudyID))+1, 
                dimnames=list(row.names(pathwaytable), c(unique(meta$StudyID), 'all')))

for (i in rownames(pathwaytable)){
  for (s in unique(meta$StudyID)){
    x <- pathwaytable[i, meta %>% filter(StudyID==s) %>% 
                        filter(group=='Healthy') %>% pull(Run)]
    x <- as.numeric(x)
    y <- pathwaytable[i, meta %>% filter(StudyID==s) %>% 
                        filter(group=='Obese') %>% pull(Run)]
    y <- as.numeric(y)
    p_value[i,s] <- wilcox.test(x, y, exact=FALSE)$p.value
  }
}

for (i in rownames(pathwaytable)){
  y=as.numeric(pathwaytable[i,])
  x=as.factor(meta$group)
  block=as.factor(meta$StudyID)
  df <- data.frame(y,x,block)
  p_value[i,'all'] <- pvalue(wilcox_test(y ~ x | block, data=df))
}

fdrvalue <- data.frame(apply(p_value, MARGIN=2, FUN=p.adjust, method='fdr'),
                    check.names = FALSE)

p_value <- p_value %>% as.data.frame() %>%  rownames_to_column("pathid")
write.table(p_value,file="pathway_wilcoxon_combined.txt",quote=F, sep='\t',row.names=F)

fdrvalue <- fdrvalue %>% as.data.frame() %>%  rownames_to_column("pathid")
write.table(fdrvalue,file="pathway_fdr_combined.txt",quote=F, sep='\t',row.names=F)

############################ the prevalance of each function #####################
upsetdf<-
  pathwaytable %>%
  rownames_to_column("pathid")%>%
  gather(-pathid, key=Run, value=Count) %>% 
  mutate(Count=if_else(Count==0, 0, 1)) %>% 
  #if an OTU is not zero, then it will be labeled 1
  left_join(meta[,c("Run","StudyID")])%>%
  select(-Run)

uplot <- upsetdf%>% group_by(StudyID, pathid)%>%
  dplyr::summarise(Count=max(Count))%>%
  spread(key=StudyID, value=Count, fill=0) %>%
  as.data.frame()%>%
  column_to_rownames("pathid")

uplot$count <- apply(uplot,1,sum)
uplot <- as.data.frame(uplot)%>%rownames_to_column("pathid")
write.table(uplot,file="path_upsetRDP.csv",quote=F, sep=',',row.names=F)
