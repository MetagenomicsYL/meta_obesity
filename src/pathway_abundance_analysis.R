##################### fuctional differential abundance analysis #####################
pkg <- c("tidyverse", "lme4", "cowplot", "coin", "ggthemes", "ggpattern")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}
source('F:/Rproject/pic/meta_code/function_library.R')

setwd("F:/Rproject/pic/meta_dataset")

mapdata <- read.csv("MetaCyc_pathway_map_update.tsv",sep="\t")
funct_table <- read.csv("combined_pathway.txt",sep = "\t")
funct_table <- column_to_rownames(funct_table, "pathway")

metagroup <- read.csv("meta_BMI_obesity.csv") %>%
  dplyr::rename(Samplename=SampleID,SampleID=Run)
funct_table <- funct_table[,metagroup$SampleID]
dropid <- c("PWY-5634", "PWY-5679", "PWY-7022", "PWY-7401")
funct_update <- subset(funct_table, !(rownames(funct_table) %in% dropid))

fun_table<-zCompositions::cmultRepl(t(funct_update),
                                    label=0,
                                    method="CZM",
                                    output="p-counts")%>%
  as.data.frame() %>%
  rownames_to_column("SampleID")

fdrpathway <-read.csv("pathway_fdr_combined.txt",sep="\t") %>% 
  column_to_rownames("pathid")

pathwaysig <- subset(fdrpathway,all<0.05) %>% rownames()

selectfun <- fun_table[,c("SampleID",pathwaysig)]%>%as.data.frame()%>%
  left_join(metagroup)

############## calculate the random effect model ################################

fungather <- gather(selectfun, key=pathid, value=reads, -colnames(metagroup))

funlog <- cal_logfun(metagroup, fungather, "pathid")

eachst <- cal_log2fc(metagroup, funlog, "pathid")

totaldf <- combine_log(eachst)

funCombined <- combine_logfc(totaldf, "pathid")

uplot <- read.csv("path_upsetRDP.csv")
pathcount <- uplot[, c("pathid","count")]
funjoin <- plyr::join(funCombined, pathcount, by="pathid", type="left")

fdrvalue <- fdrpathway %>% rownames_to_column("pathid")
fdrvalue <- fdrvalue[,c("pathid","all")]
funjoinfdr <- plyr::join(funjoin, fdrvalue, by="pathid", type="inner")
combsub <- subset(funjoinfdr, Pvalue<0.05) %>% filter(count>9)%>% 
  filter(all<0.01)
pathwaysub <- subset(fdrpathway, rownames(fdrpathway) %in% combsub$pathid) %>% 
  rownames_to_column("mapid")

pathwaymap <- plyr::join(pathwaysub, mapdata, by="mapid", type="inner")

mapdata$pval <- pathwaymap[match(mapdata$mapid, pathwaymap$mapid), 'all']
mapdataex <- filter(mapdata, !is.na(pval))
module.name <- mapdata[,c("mapid","pathway")]%>%as.data.frame()

carbohydrate.degradation <- mapdataex %>% 
  filter(Superclass2=='Carbohydrate Degradation') %>% 
  pull(mapid)

amino.acid.modules <- mapdataex %>% 
  filter(Superclass2=='Amino Acid Degradation') %>% 
  pull(mapid)

Fatty.Acid <- mapdataex %>% 
  filter(Superclass2=='Fatty Acid and Lipid Biosynthesis') %>% 
  pull(mapid)

Vitamin.Biosynthesis <- mapdataex %>% 
  filter(Superclass2=='Cofactor, Carrier, and Vitamin Biosynthesis') %>% 
  pull(mapid)

df.category <- tibble(
  `Amino acid degradation` = 
    colMeans(log10(funct_update[amino.acid.modules,] + 1e-05), na.rm=TRUE),
  `Fatty Acid and Lipid Biosynthesis` = 
    colMeans(log10(funct_update[Fatty.Acid,] + 1e-05), 
             na.rm=TRUE),
  `Carbohydrate degradation` = 
    colMeans(log10(funct_update[carbohydrate.degradation,] + 1e-05), na.rm=TRUE),
  `Cofactor, Carrier, and Vitamin Biosynthesis` =
    colMeans(log10(funct_update[Vitamin.Biosynthesis,] + 1e-05), na.rm=TRUE),
  Group=metagroup$group,
  Study=metagroup$StudyID)

df.temp <- df.category %>% 
  gather(key=key, value=value, -Group, -Study)

totalpval <- tibble()
for (i in colnames(df.category)[-c(5,6)]){
  t <- wilcox_test(value~as.factor(Group)|as.factor(Study), data=df.temp %>% 
                     filter(key==i))
  totalpval <- rbind(totalpval, tibble(pathway=i,pval=signif(pvalue(t),digits = 2)))
}
totalpval$pval
totalpval$pathway

datalabel = data.frame(group = totalpval$pathway,
                        Group=c("Obese","Obese","Obese","Obese"),
                        y = c(4.5, 4.5, 4.5, 4.5),
                        lab = c("P=1.2e-07", "P=6.2e-06", "P=1.5e-13", "P=1.1e-08"))

superlevels <- c('Amino acid degradation', 'Carbohydrate degradation', 'Fatty Acid and Lipid Biosynthesis','Cofactor, Carrier, and Vitamin Biosynthesis')

df.figure <- df.category %>% 
  gather(key=group, value=value, -Group, -Study)%>% 
  mutate(Group=factor(Group, levels=c("Obese","Healthy"))) %>% 
  mutate(group=factor(group, levels=superlevels))

category.plot <- ggplot(df.figure,aes(x=group, fill=Group, y=value)) + 
  geom_boxplot(outlier.shape = NA)+
  theme_classic() + 
  xlab('') + ylab('Log10(Normalized Abundance)')+ 
  scale_fill_manual(values=c('#fbb4ae',"#b3cde3"))+
  scale_color_manual(values=c('#fbb4ae',"#b3cde3"))+
  geom_text(aes(x=group,y=y,label=lab),data=datalabel,size = 2.5)+
  scale_x_discrete(labels=function(x) str_wrap(x, width = 18))

ggsave("F:/Rproject/pic/meta_fig/pathway_boxplot_FDR0.01_count9_four.pdf", category.plot, width = 7, height = 4)


######################## histo plot #############################################
carbohydrate.degradation.name <- mapdataex %>% 
  filter(Superclass2=='Carbohydrate Degradation') %>% 
  pull(pathway)

amino.acid.modules.name <- mapdataex %>% 
  filter(Superclass2=='Amino Acid Degradation') %>% 
  pull(pathway)

Fatty.Acid.name <- mapdataex %>% 
  filter(Superclass2=='Fatty Acid and Lipid Biosynthesis') %>% 
  pull(pathway)

Vitamin.Biosynthesis.name <- mapdataex %>% 
  filter(Superclass2=='Cofactor, Carrier, and Vitamin Biosynthesis') %>% 
  pull(pathway)

order <- c(amino.acid.modules.name, carbohydrate.degradation.name, Fatty.Acid.name, Vitamin.Biosynthesis.name)

pathwayhisto <- pathwaymap %>% filter(pathway %in% order) %>% left_join(module.name)

fdrorder<- subset(fdrpathway,rownames(fdrpathway) %in% pathwayhisto$mapid)%>%arrange(all) %>%rownames_to_column("mapid") %>% 
  left_join(module.name)

histodf <- tibble(pathway=factor(fdrorder$pathway,levels = rev(order)), p.vals=-log10(fdrorder$all))
histo.plot <- histodf%>% 
  ggplot(aes(x=p.vals, y=pathway)) + 
  geom_bar(stat='identity') + 
  theme_classic() + 
  xlab('') +
  ylab('-log10(FDR)')+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        # axis.title.x=element_text(size=14),
        # title=element_text(size=14),
        panel.background = element_rect(fill=NULL, colour = 'black'))+
  # panel.border = element_rect(fill=NA,color="white", size=1, linetype="solid"))+
  scale_y_discrete(position='right')

##########################lollipop plot#############################################
fourmoduleid <- c(amino.acid.modules,carbohydrate.degradation,Fatty.Acid,Vitamin.Biosynthesis)

four_module <- funct_update[fourmoduleid,] %>%
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  left_join(metagroup[,c("SampleID","group")]) %>%
  gather(key = mapid,value=abundance,-SampleID,-group) 
four_modulelog <- four_module%>% mutate(log_abundance=log10(four_module[,"abundance"] + 1e-05))

four_modulegroup <- four_modulelog %>% 
  group_by(mapid,group) %>% 
  summarise(median_log_abundance=mean(log_abundance)) %>% 
  left_join(mapdataex)

four_modulegroup$group <- factor(four_modulegroup$group, levels=c("Obese","Healthy"))
four_modulegroup$pathway <- factor(four_modulegroup$pathway, levels = rev(order))
lollipop <- ggplot(four_modulegroup,aes(x=median_log_abundance, y=pathway,fill=Superclass2,shape=group)) + 
  geom_linerange(aes(xmin = -5, xmax = median_log_abundance, y = pathway,color = Superclass2),position = position_dodge(width = 0.75), size = 0.5)+
  geom_point(position = position_dodge(width=0.75),size=2)+
  scale_shape_manual(values = c(21,24))+
  theme_classic() + 
  xlab('Mean log10(Normalized Abundance)') + ylab('')+ 
  scale_color_manual(values=c('#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#fff2ae','#e6f5c9'))+
  scale_fill_manual(values=c('#b3e2cd','#fdcdac','#cbd5e8','#f4cae4','#fff2ae','#e6f5c9'))+
  theme(
    legend.position = "left",
    axis.title.x=element_text(size=9))


final.plot <- plot_grid(lollipop, histo.plot,  align = 'v', nrow=1, rel_widths = c(1,0.1))
ggsave("F:/Rproject/pic/meta_fig/pathway_lollipop_fourmodule.pdf", final.plot, width =11, height = 5)

################ save significant functional pathways ############################
allsigmodule <- funct_update[combsub$pathid,] %>% t() %>%
  as.data.frame() %>% rownames_to_column("SampleID") %>% 
  left_join(metagroup[,c("SampleID","group")]) %>% 
  gather(key = mapid,value=abundance,-SampleID,-group)

allsigmodule <- allsigmodule %>% 
  mutate(log10_abundance=log10(allsigmodule[,"abundance"] + 1e-05))%>% 
  group_by(mapid,group) %>% 
  summarise(mean_log10_abundance=mean(log10_abundance)) %>% 
  left_join(mapdataex)

allsigtable <- spread(allsigmodule,key = group, value = mean_log10_abundance)%>% 
  rename(FDR=pval,Healthy_logAbundance=Healthy,Obese_logAbundance=Obese,pathid=mapid) %>% 
  left_join(combsub[,c("pathid","Pvalue")]) %>% 
  rename(REM_Pval=Pvalue)

write.table(allsigtable,file="pathway_table_39modules.txt",quote=F, sep='\t',row.names=F)
