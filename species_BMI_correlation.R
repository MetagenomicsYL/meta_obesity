######################### significant species from HKgutmicmap #######################
pkg <- c("tidyverse", "coin", "boot", "cowplot", "ggsci", "ggpubr", "ggthemes")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}

setwd("F:/Rproject/pic/meta_dataset")

meta <- read.csv("meta_BMI_obesity.csv")

HKGut <- read.csv("metaphlan3_relabun_L7.txt",sep = "\t")
mapid <- read.csv("map_id_table.csv")
mapid <- rename(mapid,sample=SampleID)
metahk <- meta %>% filter(StudyID=="HKGutMicMap")
HKGut <- HKGut[,-1]
HKGut <- HKGut %>% separate(taxonomy, 
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
deplete_genera <- read.csv("bubble_depleted_genera.csv")

nameasdown <- subset(deplete_genera, Study=="Asian") %>% filter(beta<0) %>% filter(Pval<0.05)
downname <- paste0("g__",nameasdown$Genus)

###################### species BMI netagive association #########################

hksptable <- subset(HKGut, Genus %in% downname)
hksp <- hksptable[,-which(colnames(hksptable) %in% c("taxonomy","Kingdom",
                                                     "Phylum",
                                                     "Class",
                                                     "Order",
                                                     "Family",
                                                     "Genus"))]
rownames(hksp) <- hksp$Species
hksp <- hksp[,-length(hksp)]
hksp <- hksp%>% t() %>% 
  as.data.frame() %>% rownames_to_column("sample")

mapid$Run <- gsub("-","",mapid$Run)
metasampleID <- plyr::join(metahk, mapid, by="Run",type="inner")
metasampleID <- metasampleID[,c("sample","group","BMI")]
hksptable_meta <- plyr::join(hksp,metasampleID,by="sample",type="inner")

hksptable.plot <- gather(hksptable_meta, key="species",value = "Abundance", -sample,-group,-BMI)

getcor <- function(x, ndx) {
  return(cor(x[ndx,1], x[ndx,2], method="spearman"))
}
calculate_spearman_CI <- function(hksptable.plot){
  spci <- tibble()
  hksptable.plot$Abundance <- as.numeric(hksptable.plot$Abundance)
  for (i in unique(hksptable.plot$species)){
    subtable <- subset(hksptable.plot, species==i)
    if(sum(subtable$Abundance)!=0){
      af <- subtable[,c("Abundance","BMI")]
      r <- cor(subtable$Abundance,subtable$BMI, method="spearman")
      result <- boot(af,getcor, R=1000, stype="i")
      CI <- boot.ci(result, type="bca")
      low_ci <- CI$bca[4]
      high_ci <- CI$bca[5]
      spci <- bind_rows(spci, tibble(
        Study="HKGutMicMap",species=i,low=low_ci, SpearmanR=r, high=high_ci
      ))
    }
  }
  return(spci)
}

spcidown <- calculate_spearman_CI(hksptable.plot)

spcinega <- spcidown %>% filter(SpearmanR<0)
sigsp <- spcidown %>% filter(SpearmanR<0) %>% pull("species")

calculate_wilcoxon <- function(hksptable.plot,hksptable_meta){
  wilhk <- tibble()
  for (i in unique(hksptable.plot$species)){
    t <- wilcox.test(as.numeric(hksptable_meta[,i])~hksptable_meta[,"group"])
    tem <- tibble(species=i,pval=t$p.value)
    wilhk <- rbind(wilhk,tem)
  }
  return(wilhk)
}

wilhkdown <- calculate_wilcoxon(hksptable.plot, hksptable_meta)

negative <- plyr::join(spcinega,wilhkdown,by="species", type="inner") %>% 
  mutate(Wilcoxon_sig=if_else(pval<0.05,"true","false")) %>% 
  mutate(Spearman_sig=if_else(high<0,"true","false"))
negative$species <- gsub("s__","",negative$species)

BMI_nega_sp <- ggplot(aes(x=SpearmanR, y=species),data=negative) + 
  geom_vline(xintercept =0, linetype="dashed", color="grey50") +
  geom_linerange(aes(xmin=low, xmax=high), colour='darkgrey') + 
  geom_point(aes(shape=Wilcoxon_sig,color=Spearman_sig)) +
  theme_classic() + 
  scale_shape_manual(values=c(1,16))+
  scale_color_manual(values=c('#BDCD00','#006165'))+
  theme(legend.position = "right")+
  guides(shape = guide_legend(override.aes = list(size = 1.5)),
         color = guide_legend(override.aes = list(size = 1.5)))+
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 9),
        axis.text.y = element_text(size = 8))+
  xlab("Spearman correlation")+
  ylab("")

hksptable_meta <- hksptable_meta%>%mutate(BMI_interval=case_when(
  BMI<15.5 ~ 15,
  BMI>42 ~ 42,
  TRUE ~ round(BMI,digits=0)
))
select_sp <- subset(negative, (Wilcoxon_sig == "true") & (Spearman_sig == "true"))
select_spname <- select_sp$species
select_spname <- paste0("s__",select_spname)

sp_plot <- gather(hksptable_meta, key = "Species",value = "Abundance", -sample, -BMI, -BMI_interval, -group) %>%
  subset(Species %in% select_spname)

sp_plot$Abundance <- as.numeric(sp_plot$Abundance)
sp_group <- sp_plot%>%
  group_by(Species,BMI_interval,group)%>%
  dplyr::summarize(Mean_Abandance=mean(Abundance), sd=sd(Abundance)) %>%
  ungroup()

sp_group$Species <- gsub("s__", "", sp_group$Species)
sigsp <- gsub("s__", "", select_spname)

df_anno <- data.frame(Species=sigsp,
                BMI_interval=c(26, 26, 26, 26),
                sd=c(0,0,0,0),
                Mean_Abandance=c( 0.065, 6.5, 0.125, 0.76),
                label=c("Correlation = -0.16, P = 0.02","Correlation = -0.16, P = 0.02","Correlation = -0.15, P = 0.03","Correlation = -0.14, P = 0.04"))

sp_asso <- ggplot(sp_group,aes(x=BMI_interval, y=Mean_Abandance,ymin=Mean_Abandance-sd, ymax=Mean_Abandance+sd))+
  geom_smooth(method = "glm",fill="#bababa",color="gray60") +
  geom_point(aes(color=group))+
  geom_errorbar(width=0,aes(color=group))+
  facet_wrap(~Species, scales="free", nrow = 1) +
  theme_few() +
  theme(strip.text =element_text(size = 9))+
  xlab("BMI") +
  ylab("Relative Abandance %") +
  scale_color_manual(values=c('#a6cee3','#fb9a99'))+
  geom_text(data=df_anno,aes(x=BMI_interval,y=Mean_Abandance,label=label),size=2.5)

ggsave("F:/Rproject/pic/meta_fig/BMI_sig_species_pointline.pdf", sp_asso, width = 11, height = 5)

###################### species BMI positive association #########################
enrich_genera <- read.csv("bubble_enriched_genera.csv")

nameasup<- subset(enrich_genera, Study=="Asian") %>% filter(Genus!="Roseburia")
upname <- nameasup$Genus
upname <- paste0("g__",upname)

hksptableup <- subset(HKGut, Genus %in% upname)

hkspup <- hksptableup[,-which(colnames(hksptableup) %in% c("taxonomy","Kingdom",
                                                           "Phylum",
                                                           "Class",
                                                           "Order",
                                                           "Family",
                                                           "Genus"))]
rownames(hkspup) <- hkspup$Species
hkspup <- hkspup[,-length(hksp)]
hkspup <- hkspup%>% t() %>% 
  as.data.frame() %>% rownames_to_column("sample")

hksptableup_meta <- plyr::join(hkspup,metasampleID,by="sample",type="inner")
hksptableup.plot <- gather(hksptableup_meta, key="species",value = "Abundance", -sample,-group,-BMI)

spciup <- calculate_spearman_CI(hksptableup.plot)

wilhkup <- calculate_wilcoxon(hksptableup.plot, hksptableup_meta)

positive <- plyr::join(spciup,wilhkup,by="species", type="inner") %>% 
  mutate(Wilcoxon_sig=if_else(pval<0.05,"true","false")) %>% 
  mutate(Spearman_sig=if_else(high<0,"true","false"))

positive$species <- gsub("s__","",positive$species)

BMI_posi_sp <- ggplot(aes(x=SpearmanR, y=species),data=positive) + 
  geom_vline(xintercept =0, linetype="dashed", color="grey50") +
  geom_linerange(aes(xmin=low, xmax=high), colour='darkgrey') + 
  geom_point(aes(shape=Wilcoxon_sig,color=Spearman_sig)) +
  theme_classic() + 
  scale_shape_manual(values=c(1,16))+
  scale_color_manual(values=c('#BDCD00','#006165'))+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 9),
        axis.text.y = element_text(size = 8))+
  xlab("Spearman correlation")+
  ylab("")

speciesforest <- plot_grid(BMI_nega_sp, BMI_posi_sp, align = 'v',nrow=1,rel_widths = c(0.6,0.4))
ggsave("F:/Rproject/pic/meta_fig/forest_BMI_species_all.pdf", speciesforest, width = 10, height = 6)