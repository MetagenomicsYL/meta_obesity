########################### volcano plot for REM and Wilcoxon #######################
pkg <- c("tidyverse","cowplot","ggthemes")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}

setwd("F:/Rproject/pic/meta_dataset")

remp <- read.csv("REM_P_rdp_study_count12.csv")
remfc <- read.csv("REM_FC_rdp_study_count12.csv")
wilfdr <- read.csv("p_adj.csv")
wilfdr$Genus <- gsub("g__","",wilfdr$Genus)
wilfdr <- wilfdr%>%rename(genuslevel=Genus)
fdr <- wilfdr[,c("genuslevel","all")]
fc <- remfc[,c("genuslevel","Combined")]%>%rename(FC=Combined)
rem <- plyr::join(remp,fc,by="genuslevel",type="inner")%>%left_join(fdr)
rem <- rem%>%rename(REM=Combined,FDR=all)

write.table(rem,file="REM_FC_fdr_count12.csv",quote=F, sep=',',row.names=F)

markerset <- read.csv("marker_wilpval_count12.csv")
markerset <- markerset%>%rename(Pvalue=all)
wilpval <- markerset[,c("genuslevel","Pvalue")]
wilpval$genuslevel <- gsub("g__","",wilpval$genuslevel)
remWilp <- plyr::join(rem,wilpval,by="genuslevel",type="inner")

write.table(remWilp,file="REM_FC_P_count12.csv",quote=F, sep=',',row.names=F)

valcano <- remWilp %>% mutate(logrem=(-log10(REM))) %>% 
  mutate(Significance=case_when(REM<0.05 ~ "true", 
                                TRUE ~ "false")) %>% 
  mutate(Direction=if_else(FC>0,"up","down"))

val.plot <- ggplot(aes(FC, logrem),data=valcano)+
  geom_point(aes(color =Significance,shape=Direction))+
  geom_vline(xintercept = 0, col='#525252',linetype="dashed")+
  geom_hline(yintercept = valcano %>% 
               filter(REM <=0.05) %>% 
               pull(logrem) %>% min, col='#525252',linetype="dashed")+
  ylab('REM -log10(P-value)') + 
  xlab('log2FC')+ 
  scale_shape_manual(values=c(17,16))+
  scale_colour_manual(values=c('grey','#BDCD00'))+
  theme_few()+
  theme(legend.position = "bottom",
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=9),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"))

marker.plot <- remWilp %>%
  mutate(type=case_when(FDR<0.001 & REM<0.05 ~ 'Both',
                        FDR>=0.001 & REM<0.05 ~'REM',
                        FDR<0.001 & REM>=0.05 ~'FDR',
                        TRUE~ 'Neither'))%>%
  mutate(type=factor(type, levels =
                       c('Both', 'REM', 'FDR','Neither')))


marker.plot$logP <- -log10(marker.plot$REM)
marker.plot$logfdr <- -log10(marker.plot$FDR)

corspearman <- cor.test(marker.plot$logP, marker.plot$logfdr,alternative = "two.sided",method = "spearman",conf.level = 0.95)

marker.fig <- ggplot(marker.plot,aes(x=logfdr, y=logP, col=type)) +
  geom_point()+
  geom_smooth(method = "glm",fill="grey80", color="grey50",linetype=2) +
  geom_hline(yintercept = marker.plot %>% 
               filter(REM <=0.05) %>% 
               pull(logP) %>% min, col='#525252',linetype="dashed")+
  geom_vline(xintercept = marker.plot %>% 
               filter(FDR <=0.001) %>% 
               pull(logfdr) %>% min, col='#525252',linetype="dashed") +
  ylab('REM -log10(P-value)') + 
  xlab('Wilcoxon -log10(FDR)')+ 
  scale_colour_manual(values=c('#006165','#BDCD00', 'grey'), 
                      name='Significance')+
  annotate('text', x=5, y=11, label=("Correlation=0.925  P=2.2e-16"), size= 3)+
  theme_few()+
  theme(legend.position = "bottom",
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=9),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"))

point.graph <- plot_grid(val.plot, marker.fig, align = 'v', nrow=1,rel_widths = c(0.5,0.5))

ggsave("F:/Rproject/pic/meta_fig/valcano_marker_fig.pdf",point.graph, width = 9, height =4)

############################# forest plot for significant genus ##################

meta <- read.csv("meta_BMI_obesity.csv")

fdrtable <- read.csv("p_adj.csv")%>%
  mutate(genuslevel=gsub("g__","", Genus))%>%
  select(-Genus) %>% column_to_rownames("genuslevel")

submarker <- subset(remWilp, REM<0.05)%>% filter(FDR<0.001)

order <- submarker%>%arrange(FDR)
genusorder <- order$genuslevel

fdr.plot <- fdrtable[genusorder,]
fdr.fig <- tibble(species=factor(rownames(fdr.plot),levels = rev(genusorder)), p.vals=-log10(fdr.plot$all))

fdr.histo <- fdr.fig %>% 
  ggplot(aes(x=species, y=p.vals)) + 
  geom_bar(stat='identity') + 
  theme_classic() + 
  xlab('significance') +
  ylab('-log10(FDR)') + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill=NULL, colour = 'black'))+
  scale_x_discrete(position='top')+
  coord_flip()


log2fcvalue <- read.csv("OTU_forest_plot_count12.csv")
combined <- subset(log2fcvalue, StudyID=="Combined")
subcom <- subset(combined, genuslevel %in% submarker$genuslevel)
subcom$genuslevel <- factor(subcom$genuslevel, levels=rev(genusorder))

forest <- ggplot(aes(x=log2FC, y=genuslevel,fill=Significance),data=subcom) + 
  geom_vline(xintercept =0, linetype="dashed", color="grey50") +
  geom_linerange(aes(xmin=CI_low, xmax=CI_high), colour='darkgrey') + 
  geom_point(pch=23) +
  theme_classic() + 
  scale_fill_manual(values=c('#b3cde3','#fbb4ae'),name='Direction')+
  theme(legend.position = "left")+
  ylab("")

forest_histo.plot <- plot_grid(forest, fdr.histo,  align = 'v', nrow=1, rel_widths = c(1,0.2))
ggsave("F:/Rproject/pic/meta_fig/Forest_fdr_genus_rel.pdf", forest_histo.plot, width = 8, height =4)
