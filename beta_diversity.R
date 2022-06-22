############################# calculate for the beta diversity ########################
pkg <- c("labdsv","coin","vegan","ggpubr","tidyverse", "cowplot","ggsci")
for(i in pkg){
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}

setwd("F:/Rproject/pic/meta_dataset")
genusOTU <- read.csv("genus_combinedOTU_rdp.csv")%>%
  column_to_rownames("Genus")
meta_obe <- read.csv("meta_BMI_obesity.csv")
genusOTU <- genusOTU[,meta_obe$Run]

#take a long time to calculate the distance

dist = vegdist(t(genusOTU), method = 'bray')
pco.results = pco(dist, k=2)
axis.1.title <- paste('PCo1 [', 
                      round((pco.results$eig[1]/sum(pco.results$eig))*100,1),
                      '%]', sep='')

axis.2.title <- paste('PCo2 [', 
                      round((pco.results$eig[2]/sum(pco.results$eig))*100,1),
                      '%]', sep='')
df.plot <- tibble(Axis1 = -1*pco.results$points[,1],
                  Axis2 = pco.results$points[,2],
                  Run = rownames(pco.results$points))
metasub <- meta[,c("Run","group","StudyID")]
df <- inner_join(df.plot,metasub,by="Run")


study.cols=c('#7bccc4','#ffff99','#33a02c','#fb9a99','#d4b9da','#fdbf6f','#ff7f00','#cab2d6','#a6bddb','#a6cee3','#bebada','#fb8072','#fdb462','#b3de69','#fccde5','#d9d9d9','#ccebc5','#d0d1e6')
w1 <- wilcox.test(df[,"Axis1"][df[,"group"]=="Obese"], df[,"Axis1"][df[,"group"]=="Healthy"])
summary(w1)
w2 <- wilcox.test(df[,"Axis2"][df[,"group"]=="Obese"], df[,"Axis2"][df[,"group"]=="Healthy"])
summary(w2)

#take a long time for the statistical test

anosim.result.group<-anosim(dist,df$group,permutations = 999, parallel=16)
summary(anosim.result.group)
anosim.result.study<-anosim(dist,df$StudyID,permutations = 999,parallel=16)
summary(anosim.result.study)

gs <- df %>% 
  ggplot(aes(x=Axis1, y=Axis2, shape=group, col=group)) +
  geom_point(alpha=0.6,shape=16) + 
  scale_colour_manual(values=c("#a6cee3", "#1f78b4")) + 
  scale_x_continuous(position='top') +
  xlab(axis.1.title) + ylab(axis.2.title) +
  annotate("text",x=0.25,y=0.5,label="ANOSIM R=0.042 P=0.001")+
  theme(panel.background = element_rect(fill='white', color = 'black'),
        legend.key = element_rect(fill="transparent"),
        legend.position = c(0.9,0.8),
        axis.ticks=element_blank(), 
        axis.text = element_blank(),
        panel.grid = element_blank())

g.g.1 <- 
  df %>% 
  ggplot(aes(x=group, y=Axis1, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#a6cee3", "#1f78b4"), guide=FALSE) + 
  ylab("PCo1")+
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank()) +
  coord_flip()+
  xlab("P=4.7e-6")

g.g.2 <- 
  df %>% 
  ggplot(aes(x=group, y=Axis2, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#a6cee3", "#1f78b4"), guide=FALSE) + 
  scale_x_discrete(position='top') + 
  scale_y_continuous(position = 'right') +
  ylab("PCo2")+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank())+
  xlab("P=3.2e-61")

g2 <- df %>% 
  ggplot(aes(x=Axis1, y=Axis2, shape=group, col=StudyID)) +
  geom_point(alpha=0.6,shape=20) + 
  scale_colour_manual(values=study.cols)+
  scale_x_continuous(position='top') +
  xlab(axis.1.title) + ylab("") +
  annotate("text",x=0.25,y=0.5,label="ANOSIM R=0.302 P=0.001")+
  theme(panel.background = element_rect(fill='white', color = 'black'),
        legend.text=element_text(size=6),
        legend.title=element_text(size=9),
        legend.key = element_rect(fill="transparent"),
        axis.ticks=element_blank(), axis.text = element_blank(),
        legend.position="none",
        panel.grid = element_blank())

combined <- plot_grid(gs, g2, nrow=1, rel_widths = c(0.5, 0.5))

plottotal <- plot_grid(combined, g.g.2, g.g.1,
                   nrow=2,
                   rel_widths = c(0.9, 0.1), rel_heights = c(0.85, 0.15))

ggsave("F:/Rproject/pic/meta_fig/PcoAbeta.pdf",plottotal, width = 13, height = 7)


