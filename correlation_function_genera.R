#########################correlation between function and genera#####################
reqpkg <- c("tidyverse", "coin", "ComplexHeatmap", "RColorBrewer", "circlize")
for(i in reqpkg)
{
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}

setwd("F:/Rproject/pic/meta_dataset")

meta <- read.csv("meta_BMI_obesity.csv")
dropsample <- c("ERR1072662", "ERR1074388", "ERR1077871", "ERR1090234", "ERR1250429", "ERR1076833", "SRR1796120", "ERR2993019") 
metagroup <- subset(meta, !(Run %in% dropsample))
pathwaytable <- read.csv("combined_pathway.txt",sep = "\t") %>%
  column_to_rownames("pathway")
pathwaytable <- pathwaytable[,metagroup$Run]

mapdata <- read.csv("MetaCyc_pathway_map_update.tsv",sep="\t")

genustable <- read.csv('genus_combined_rel_rdp.csv') %>%
  column_to_rownames("Genus")
genustable <- genustable[,metagroup$Run]
########################### significant genera ###################################
markerset <- read.csv("pFCfdr_REM_count.csv") %>%
  filter(count>12) %>%
  filter(all<0.001) %>%
  filter(Combined<0.05) %>% 
  pull(genuslevel)
fdrotu <-read.csv("p_adj.csv") %>%
  column_to_rownames("Genus")
siggenusname <- fdrotu[fdrotu$all < 0.001,] %>% rownames()
markerset <- paste0("g__",markerset)
siggenus <- genustable[markerset,]

############################### significant pathway ###############################

sigpath <- read.csv("pathway_table_39modules.txt", sep = "\t")

pathwaysig <- sigpath$pathid
pathwaysub <- pathwaytable[pathwaysig,]

siggenus_czm<-zCompositions::cmultRepl(t(siggenus),
                                    label=0,
                                    method="CZM",
                                    output="prop") %>% t()
loggenus <- log(siggenus_czm)
rownames(loggenus) <- gsub("g__","",rownames(loggenus))

pathwaysub_czm<-zCompositions::cmultRepl(t(pathwaysub),
                                      label=0,
                                      method="CZM",
                                      output="p-count") %>% t()
logpathway <- log(pathwaysub_czm)

corgenus_path <- cor(t(loggenus), t(logpathway), method="spearman")
rownames(corgenus_path) <- gsub("g__", "", rownames(corgenus_path))

pathwayname <- mapdata[,c("mapid","pathway")]
cordata <- t(corgenus_path) %>%
  as.data.frame() %>%
  rownames_to_column("mapid") %>%
  left_join(pathwayname) %>% 
  column_to_rownames("pathway")

cordata <- cordata[,-1]

####################### heatmap for the genus and function #######################

anndata <- data.frame(mapdata %>% select(mapid, Superclass2,pathway)) %>%
  filter(mapid %in% pathwaysig)
rownames(anndata) <- anndata$pathway
anndata$mapid <- NULL
anndata$pathway <- NULL

pathwaycols<-list(`Amino Acid Degradation`="#b3e2cd",`Carbohydrate Degradation`="#fdcdac",`Fatty Acid and Lipid Biosynthesis`="#f4cae4",`Cofactor, Carrier, and Vitamin Biosynthesis`="#cbd5e8")

anndata$Superclass2[!anndata$Superclass2 %in% names(pathwaycols)] <- NA
anndata <- rename(anndata, type=Superclass2)
cordata <- as.matrix(cordata)
order <- rownames(cordata)

anndataorder <- anndata[order,]
anndatadf <- data.frame(pathway=order,Superclass=anndataorder)
anndatadf <- plyr::arrange(anndatadf, Superclass)
anndatadf <- column_to_rownames(anndatadf, "pathway")

sidedf <- rowAnnotation(
  df = anndatadf, 
  col=list('Superclass' = unlist(pathwaycols)))
cordata <- cordata[rownames(anndatadf), ]

heatmp <- Heatmap(cordata, 
                name = "spearman correlation",
                show_column_names = TRUE,
                right_annotation  = sidedf,
                show_row_names = TRUE,
                row_names_max_width = max_text_width(
                  rownames(cordata), 
                  gp = gpar(fontsize = 14)),
                column_names_max_height =max_text_width(
                  colnames(cordata), 
                  gp = gpar(fontsize = 14)),
                cluster_rows = FALSE,
                # clustering_method_rows = 'mcquitty',
                clustering_method_columns = 'mcquitty', 
                col=colorRamp2(seq(from=-0.5, to=0.5, length.out = 19),
                               c(rev(brewer.pal(9, 'PuBu')), 'white', 
                                 brewer.pal(9, 'OrRd'))),
                row_names_gp = gpar(fontsize = 14),
                column_names_gp=gpar(fontsize = 16),
                show_row_dend = TRUE,
                heatmap_legend_param = list(direction = "horizontal"))

pdf(file="F:/Rproject/pic/meta_fig/correlation_plot_pathway_arrange.pdf", width=15,height=10)
draw(heatmp,heatmap_legend_side = "bottom")
dev.off()