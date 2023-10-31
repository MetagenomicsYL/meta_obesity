############################# calculate for the phyla abundance ####################
pkg <- c("tidyverse", "broom", "ggthemes", "UpSetR", "RColorBrewer")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}
source('F:/Rproject/pic/meta_code/function_library.R')

setwd("F:/Rproject/pic/meta_dataset")
OTUtable <- read.csv("combined_uniqueOTU_rdp.csv")
meta <- read.csv("meta_BMI_obesity.csv")
meta <- meta[, -which(colnames(meta) %in% c("SampleID"))]
meta$SampleID <- meta$Run

OTUtable <- OTUtable[, c("Taxon", meta$Run)]
OTUsep <- OTUtable %>% separate(Taxon, 
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

OTUPhylum <- select(OTUsep,-c("Taxon","Kingdom", "Class", "Order", "Family", "Genus", "Species"))
OTUPhylum$Phylum <- factor(OTUPhylum$Phylum)

aggdata<-aggregate(OTUPhylum[, -1], by=list(Phylum=OTUPhylum$Phylum),FUN=sum) %>%
  as.data.frame() %>%
  column_to_rownames(var="Phylum")

metasub <- meta[, c("SampleID", "StudyID","group")]

tablephyla <- prop.table(as.matrix(aggdata), 2) %>%
  as.data.frame() %>%
  rownames_to_column("Phylum")%>%
  filter(Phylum!="p__")%>%
  mutate(Phylum=gsub("p__","", Phylum))

tablephylagather <- tablephyla%>%
  gather(-Phylum, key="SampleID", value="Abundance")

tp<- tablephylagather %>% group_by(Phylum) %>% summarize(mean=mean(Abundance))%>%
  arrange(desc(mean))%>%
  top_n(9, mean)%>% 
  bind_rows(., tibble(Phylum="Other", mean=0))

table_temp <- tablephyla %>%
  as.data.frame() %>%
  column_to_rownames('Phylum') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID")

phyla_study <- plyr::join(table_temp, metasub, by="SampleID",type = "inner") %>% 
  gather(-SampleID, -StudyID, -group, key="Phylum", value="Abundance") %>% 
  mutate(Phylum=if_else(Phylum %in% tp$Phylum, Phylum, "Other")) %>%
  mutate(Phylum=factor(Phylum, levels = rev(tp$Phylum))) %>%
  group_by(StudyID,group,Phylum) %>% 
  summarise(prop=sum(Abundance))

phyla_fig <- tibble()
for (i in unique(phyla_study$StudyID)) {
  hdf <- meta %>% filter(StudyID == i) %>% filter(group == "Healthy")
  hcount <- length(hdf$SampleID)
  odf <- meta %>% filter(StudyID == i) %>% filter(group == "Obese")
  ocount <- length(odf$SampleID)
  htentative <- subset(phyla_study,StudyID==i) %>% filter(group == "Healthy")
  otentative <- subset(phyla_study,StudyID==i) %>% filter(group == "Obese")
  htentative$prop <- htentative$prop/hcount
  otentative$prop <- otentative$prop/ocount
  phyla_fig <- rbind(phyla_fig,htentative,otentative)
}

studyorder <- tablephylagather %>%
  mutate(Phylum=if_else(Phylum %in% tp$Phylum, Phylum, "Other")) %>%
  mutate(Phylum=factor(Phylum, levels = rev(tp$Phylum)))%>%
  left_join(metasub) %>% 
  group_by(Phylum,StudyID) %>% 
  summarise(prop=mean(Abundance)) %>% 
  filter(Phylum=="Firmicutes") %>% 
  arrange(desc(prop)) %>% 
  pull("StudyID")

phyla_fig$StudyID <- factor(phyla_fig$StudyID, levels = studyorder)  
phyla_fig[phyla_fig$group=="Obese",]$group <- "OB"
phyla_fig[phyla_fig$group=="Healthy",]$group <- "HT"
phyla_fig$group <- factor(phyla_fig$group, levels = c("OB","HT"))

########################################draw phyla plot#############################

phyla_plot <- ggplot(aes(x=group, y=prop, fill=Phylum),data = phyla_fig) +
  geom_bar(stat="identity") +
  facet_wrap(~StudyID, scales="free_x", nrow=1) +
  theme_few() +
  theme(
    text=element_text(family="Helvetica"),
    axis.text.x=element_text(angle=45, hjust=1,size=7),
    axis.title.y = element_text(size = 13),
    strip.text.x = element_text(angle=90,size=10))+
  xlab("") +
  ylab("Phylum Abundance (%)") +
  scale_fill_manual(values=rev(c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd')))+
  coord_cartesian(expand=F)

ggsave("F:/Rproject/pic/meta_fig/combined_phyla.pdf", phyla_plot, width = 9, height = 4.5)

########################################draw upset plot#############################
genusOTU <- read.csv("genus_combinedOTU_rdp.csv")%>%
  filter(Genus!="g__")

genusOTU <- genusOTU[,c("Genus",meta$Run)]

uplot<-
  genusOTU %>%
  as.data.frame() %>%
  dplyr::rename(OTU=Genus)%>%
  as_tibble() %>%
  gather(-OTU, key=SampleID, value=Count) %>% 
  mutate(Count=if_else(Count==0, 0, 1))%>% 
  #if an OTU is not zero, then it will be labeled 1
  left_join(meta[,c("SampleID","StudyID")])%>%
  select(-SampleID)%>% 
  group_by(StudyID, OTU)%>%
  dplyr::summarise(Count=max(Count)) %>%
  spread(key=StudyID, value=Count, fill=0) %>%
  as.data.frame()


pdf(file ="F:/Rproject/pic/meta_fig/genus_uplot.pdf", width=12, height=8, family="Helvetica",pointsize=8)
upset(uplot, nsets=length(uplot), nintersects=30, order.by="freq", show.numbers=TRUE,
      mb.ratio = c(0.55, 0.45), 
      matrix.dot.alpha=0,
      main.bar.color="#7fcdbb",
      sets.bar.color='#edf8b1',
      matrix.color = "#737373",
      text.scale=c(2, 2, 1.5, 1.5, 2, 1.3),
      queries = list(list(query = intersects, params = unique(meta$StudyID), active = T)))
dev.off()

####### genera that were detected more than 12 studies and wilcoxon p< 0.05########

uplot_select <- uplot%>%column_to_rownames("OTU")
uplot_select$count <- apply(uplot_select,1,sum)
uplot_select <- as.data.frame(uplot_select)%>%rownames_to_column("Genus")
write.table(uplot_select,file="genus_OTUupsetRDP.csv",quote=F, sep=',',row.names=F)

selectgenus <- uplot_select%>%subset(count>=12)
write.table(selectgenus,file="genus_count12_rdp.csv",quote=F, sep=',',row.names=F)

wilpval<- read.csv("p_val.csv")
selectp <- wilpval%>%subset(all<=0.05)
selectp <- selectp[,c("Genus","all")]
select_final <- plyr::join(selectgenus,selectp,by="Genus",type="inner")
write.table(select_final,file="marker_pval_12_rdp.csv",quote=F, sep=',',row.names=F)

