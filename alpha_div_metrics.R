################# calculate alpha diversity metrics #################################
pkg <- c("tidyverse", "coin","vegan")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}
setwd("F:/Rproject/pic/meta_dataset")

BFratio_fun <- function(metadata, subOTU){
  tabletaxa <- subOTU %>% rownames_to_column("Taxon")
  tabletaxa <- tabletaxa %>% separate(Taxon, 
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
  n_seqs <- apply(subOTU, 2, sum)
  b_otusum <- tabletaxa[grep("Bacteroidetes", tabletaxa$Phylum),]
  
  b_OTUtable <- b_otusum[, colnames(b_otusum) %in% c(metadata$Run)]%>%
    as.data.frame()
  rownames(b_OTUtable) <- b_otusum$Taxon
  b_countOTU <- apply(b_OTUtable, 2, sum)
  b_relab <- b_countOTU / n_seqs
  
  f_otusum <- tabletaxa[grep("Firmicutes", tabletaxa$Phylum),]
  f_OTUtable <- f_otusum[, colnames(f_otusum) %in% c(metadata$Run)]%>%
    as.data.frame()
  rownames(f_OTUtable) <- f_otusum$Taxon
  f_countOTU <- apply(f_OTUtable, 2, sum)
  f_relab <- f_countOTU / n_seqs
  bf_ratio <- f_relab / b_relab
  bf_ratio[!is.finite(bf_ratio)] <- 1e6
  bf_relab <- data.frame(cbind(b = b_relab, f = f_relab, FB = bf_ratio)) %>%
    rownames_to_column("Run")
  return(bf_relab)
}

metametric <- list()
meta <- read.csv("meta_BMI_obesity.csv")
otutable <- read.csv("combined_uniqueOTU_filter_ab_rdp.csv")
otutable <- otutable%>%
  column_to_rownames("Taxon")

for (y in unique(meta$StudyID)){
  metadata <- subset(meta, StudyID==y)
  subOTU <- otutable[, metadata$Run]
  colsums <- apply(subOTU, 2, sum)
  size <- min(colsums)
  oturare = t(rrarefy(floor(t(subOTU)), sample=size))
  oturare2 <- as.data.frame(oturare)
  
  df <- data.frame(
    observed_features = specnumber(t(oturare2)),
    shannon_entropy = vegan::diversity(oturare2, index="shannon", MARGIN=2),
    chao1 = vegan::estimateR(t(oturare2)) %>% t() %>% as.data.frame %>% rename(Chao1=S.chao1) %>% pull(Chao1)
  ) 
  
  df <- df%>%rownames_to_column("Run")
  df$pielou_evenness <- df$shannon_entropy/log(df$observed_features)
  df$StudyID <- y
  df <- df[,c("StudyID","Run","observed_features","pielou_evenness", "shannon_entropy","chao1")]
  bftable <- BFratio_fun(metadata,subOTU)
  metagroup <- metadata[,c("Run","group")]
  df <- df %>% left_join(metagroup) %>% left_join(bftable)
  metametric[[y]] <- df
}

tablemetric <- metametric%>%
  do.call(bind_rows, .)

write.table(tablemetric,file="combined_metrics_rdp.txt",quote=F, sep='\t',row.names=F)
