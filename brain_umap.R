library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(cowplot)
library(viridis)
library(ggthemes)
library(ggalluvial)

### Origin dimplot
color_origin <- c("#E2D200", "#46ACC8", "#A6CEE3", "#B40F20", "#FB9A99", "#984EA3", "#fdbd10", "#0066b2")
umap_data <- DimPlot(brain.harmony, reduction = "umap", group.by = "origin")
umap_data <- umap_data$data

ggplot() + geom_point(data = umap_data, mapping = aes(UMAP_1, UMAP_2, color = origin), size = 0.3, stroke = 0.5)+
  scale_color_manual(values = color_origin)+
  theme_classic()+
  theme(legend.text = element_text(colour = "black",size=11,face="bold"),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(size = 1, color = "black")
  ) + guides(color = guide_legend(override.aes = list(size=6)))




### Sample dimplot
color_sample <- c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F", "#E41A1C","#E7298A", "#C3BC3F", "#FF7F00", "#F1788D","#223D6C","#D20A13","#FFD121","#088247","#B037C4", "#58CDD9","#7A142C","#5D90BA","#7CC767","#91612D", "#6E568C","#E0367A","#D8D155","#64495D")
umap_data <- DimPlot(brain.harmony, reduction = "umap", group.by = "orig.ident")
umap_data <- umap_data$data

ggplot() + geom_point(data = umap_data, mapping = aes(UMAP_1, UMAP_2, color = orig.ident), size = 0.3, stroke = 0.5)+
  scale_color_manual(values = color_sample)+
  theme_classic()+
  theme(legend.text = element_text(colour = "black",size=11,face="bold"),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(size = 1, color = "black")
  ) + guides(color = guide_legend(override.aes = list(size=6)))


### Cell dimplot
color_brain <- c("#0073C2FF", "#EEC000FF", "#CD534CFF", "#00CC33FF", "#984EA3", "#8F7700FF", "#FCCDE5", "#C46DA0FF", "#F0027F","#374E55FF", "#FD6467", "#0B775E", "#FF59FF", "#FDC086", "#5B1A18", "#69C8ECFF", "#3da8e0", "#CBDF9AFF","#802268FF")
cluster_order <- c("OLG", "MG", "ASC", "OPC", "PC", "InN_1", "InN_2", "InN_3", "InN_4", "InN_5", "ExN_1", "ExN_2", "DoN", "D1-MSN", "D2-MSN", "VLMC", "ABC", "EC", "CPC")
col_cluster <- setNames(color_brain,
                        c("OLG", "MG", "ASC", "OPC", "PC", "InN_1", "InN_2", "InN_3", "InN_4", "InN_5", "ExN_1", "ExN_2", "DoN", "D1-MSN", "D2-MSN", "VLMC", "ABC", "EC", "CPC")
)
df <- DimPlot(brain.harmony, reduction = "umap")$data
cell <- df %>% group_by(ident) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

rownames(cell) <- cell$ident
A <- cell[cluster_order, ]

p1 <- ggplot(df, aes(x = UMAP_1 , y = UMAP_2 ,col = ident)) +
      geom_point(size = 0.3, shape = 16) +
      scale_color_manual("", values = col_cluster) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(size = 1, color = "black"))+
      geom_text_repel(data = A, aes(label = ident),color = "black", size = 5, point.padding = 0.3)


#####barplot
count_df <- data.frame(Cell_Type = names(table(brain.harmony@meta.data$celltype)), Cell_Count = unname(table(brain.harmony@meta.data$celltype)))[, c(1,3)]
count_df$Cell_Type <- factor(count_df$Cell_Type, levels =  rev(cluster_order))
colnames(count_df) <- c("Cell_Type", "Cell_Count")


p2 <- ggplot(count_df, aes(x = Cell_Type, y = Cell_Count, fill = Cell_Type)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'none',
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(size = 1, color = "black")) +
      geom_text(aes(label = Cell_Count), hjust = -0.05) +
      coord_flip() + ylim(0,132000) +
      scale_fill_manual(values = rev(color_brain))


####legend
df <- data.frame(x = 0, y = rev(levels(brain.harmony)), stringsAsFactors = F )
df$y <- factor(df$y, levels = df$y )

p3 <- ggplot(df, aes(x, y, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = rev(color_brain)) +
    theme_classic() +
    scale_x_continuous(expand = c(0,0)) + scale_y_discrete(position = "right") +
    theme(
        plot.margin = margin(r=0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )


plot_grid(plotlist = list(p1, p2, p3), ncol = 3, rel_widths = c(3, 1.3, 0.6))


### DotPlot
gene <- c("Cldn11", "Mag", "Plp1", "C1qa", "Csf1r", "Cx3cr1", "Aqp4", "Gja1", "Slc1a3", "Pdgfra", "Cspg4", "Tnr", "Vtn", "Atp13a5", "Rgs5", "Gad1", "Gad2", "Slc32a1", "Pvalb", "Meis2", "Vip", "Sst", "Npy", "Slc17a6",
          "Slc17a7", "Slc18a2", "Slc6a3", "Th", "Drd1", "Tac1", "Drd2", "Penk", "Slc6a13", "Slc47a1", "Cldn5", "Flt1", "Prlr", "Ttr")
marker_df <- data.frame(gene = gene, cluster =  c(rep("OLG", 3), rep("MG", 3), rep("ASC", 3), rep("OPC", 3), rep("PC", 3), rep("InN_1", 4), rep("InN_2", 1), rep("InN_3", 1), rep("InN_4", 1), rep("InN_5", 1),
                                                  rep("ExN_1", 1), rep("ExN_2", 1),  rep("DoN", 3), rep("D1-MSN", 2), rep("D2-MSN", 2), rep("VLMC", 1), rep("ABC", 1), rep("EC", 2), rep("CPC", 2)))
marker_df$cluster <- factor(marker_df$cluster, levels = cluster_order)

p5 <- DotPlot(brain.harmony,
    features = marker_df$gene
    ) + scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    RotatedAxis() +
    theme(
         panel.border = element_rect(color = "black"),
         panel.spacing = unit(1, "mm"),
         axis.title = element_blank(),
         axis.text.y = element_blank(),
         axis.text.x = element_text(face = "italic"),
         strip.text.x = element_text(angle = 90)
     )


df <- data.frame(x = 0, y =rev(levels(brain.harmony)), stringsAsFactors = F )
df$y <- factor(df$y, levels = rev(df$y) )

p6 <- ggplot(df, aes(x, y, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = color_brain) +
    theme_classic() +
    scale_x_continuous(expand = c(0,0)) +
    theme(
        plot.margin = margin(r=0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

plot_grid(p6, p5, align ="h", axis="bt", rel_widths = c(1.5, 9))



### FeaturePlot
library(reshape2)
GeneExp <- FetchData(brain.harmony, vars = c("Plp1", "Mag", "Csf1r", "Aqp4", "Gja1","Pdgfra", "Cspg4", "Vtn", "Gad1", "Slc32a1", "Pvalb", "Meis2", "Vip", "Sst", "Npy",
          "Slc17a6", "Slc17a7", "Slc18a2", "Drd1", "Drd2", "Slc6a13", "Slc47a1", "Flt1", "Ttr")
)
pc <- Embeddings(brain.harmony, reduction = 'umap') %>% data.frame()
gbid <- cbind(pc, GeneExp)
gbidlong <- melt(gbid,id.vars = c('UMAP_1','UMAP_2'),
                 value.name = 'exp',
                 variable.name = 'gene')

gbidlong$exp <- ifelse(gbidlong$exp > 3, 3,gbidlong$exp)

ggplot(gbidlong,aes(x = UMAP_1,y = UMAP_2,color = exp)) +
  geom_point(size = 0.01, show.legend = T) +
  scale_color_gradientn(colors = c("#e3e0e2","#8c304d"),name = '')  +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA, size = 2),
        axis.text = element_blank(),
        legend.text = element_text(family = "sans", face = "bold"),
        legend.title = element_text(family = "sans", face = "bold"),
        text = element_text(family = "sans", face = "italic"), axis.title = element_text(family = "sans", face = "bold")
  ) + facet_wrap(~gene, ncol = 6)




### Heatmap_marker
library(Scillus)
library(dplyr)
markers <- read.table("DEGs_celltype_brain.txt", sep = '\t', row.names = 1, header = T)
marker_top30 <- data.frame(markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC))


library(scRNAtoolVis)
AverageHeatmap(object = brain.harmony,
               markerGene = marker_top30$gene,
               annoCol = TRUE, showRowNames = F,row_title = "Top30 marker genes",
               myanCol =  c("#0073C2FF", "#EEC000FF", "#CD534CFF", "#00CC33FF", "#984EA3", "#8F7700FF", "#FCCDE5", "#C46DA0FF", "#F0027F","#374E55FF", "#FD6467", "#0B775E", "#FF59FF", "#FDC086", "#5B1A18", "#69C8ECFF", "#3da8e0", "#CBDF9AFF","#802268FF"))


### Enrichment of top30 marker genes
library(clusterProfiler)
library(org.Mm.eg.db)
allcluster_go=data.frame()
for (i in unique(marker_top30$cluster)) {
  small_gene_group=marker_top30[marker_top30$cluster==i,]
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.1,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])

select_pathway <- allcluster_go[c(1,4,783,787,1733,1734,3245,3246,4128,4138,5057,5059,5809,5812,6458,6459,7768,7769,8378,8395,9492,9505,10605,10611,11718,11719,12610,12614,13526,13528,14558,14559,15494,15495,16166,16168,17167,17190), ]
select_pathway$cluster <- factor(select_pathway$cluster, levels = levels(brain.harmony))

ggplot(select_pathway,
       aes(x=reorder(Description,-log10(pvalue)),y=-log10(pvalue), fill=cluster)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = color_brain) +
  facet_grid(cluster~., scale = 'free_y', space = 'free_y')+
  coord_flip() +
  xlab("GO term") +
  ylab("-log10(Pvalue)") +
  labs(title = "GO Terms Enrich")+
  scale_x_discrete(position = "top") + theme_few() + NoLegend()


### Vlnplot of qc
Idents(brain.harmony) <- "orig.ident"
library(ggplot2)
color_sample <- c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F", "#E41A1C","#E7298A", "#C3BC3F", "#FF7F00", "#F1788D","#223D6C","#D20A13","#FFD121","#088247","#B037C4", "#58CDD9","#7A142C","#5D90BA","#7CC767","#91612D", "#6E568C","#E0367A","#D8D155","#64495D")

VlnPlot(brain.harmony, features = c("nFeature_RNA"), pt.size = 0, cols = color_sample, log = T) +
  theme(legend.text = element_text(colour = "black",size = 11,face = "bold"),
        legend.title = element_text(colour = "black",size = 15,face = "bold"),
        axis.title.x = element_text(size = 14,face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14,face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) + ylab("Genes") + xlab("") + labs("nFeature_RNA")

VlnPlot(brain.harmony, features = c("nCount_RNA"), pt.size = 0, cols = color_sample, log = T) +
  theme(legend.text = element_text(colour = "black",size = 11,face = "bold"),
        legend.title = element_text(colour = "black",size = 15,face = "bold"),
        axis.title.x = element_text(size = 14,face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14,face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) + ylab("Counts") + xlab("") + labs("nCount_RNA")

VlnPlot(brain.harmony, features = c("percent.mt"), pt.size = 0, cols = color_sample, log = T) +
  theme(legend.text = element_text(colour = "black",size = 11,face = "bold"),
        legend.title = element_text(colour = "black",size = 15,face = "bold"),
        axis.title.x = element_text(size = 14,face = "bold"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14,face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.background = element_rect(size = 1, color = "black")
  ) + ylab("MT-genes(%)") + xlab("") + labs("percent.mt")


### Proportion

data_Young <- brain.harmony[, brain.harmony@meta.data$origin =="Young"]
data_Old <- brain.harmony[, brain.harmony@meta.data$origin =="Old"]
data_MET <- brain.harmony[, brain.harmony@meta.data$origin =="MET"]
data_NR <- brain.harmony[, brain.harmony@meta.data$origin =="NR"]
data_DQ <- brain.harmony[, brain.harmony@meta.data$origin =="D+Q"]
data_SPD <- brain.harmony[, brain.harmony@meta.data$origin =="SPD"]
data_MN <- brain.harmony[, brain.harmony@meta.data$origin =="MET+NR"]
data_MS <- brain.harmony[, brain.harmony@meta.data$origin =="MET+SPD"]

fre_Young <- table(data_Young@meta.data$celltype)/sum(table(data_Young@meta.data$celltype))
fre_Old <- table(data_Old@meta.data$celltype)/sum(table(data_Old@meta.data$celltype))
fre_MET <- table(data_MET@meta.data$celltype)/sum(table(data_MET@meta.data$celltype))
fre_NR <- table(data_NR@meta.data$celltype)/sum(table(data_NR@meta.data$celltype))
fre_DQ <- table(data_DQ@meta.data$celltype)/sum(table(data_DQ@meta.data$celltype))
fre_SPD <- table(data_SPD@meta.data$celltype)/sum(table(data_SPD@meta.data$celltype))
fre_MN <- table(data_MN@meta.data$celltype)/sum(table(data_MN@meta.data$celltype))
fre_MS <- table(data_MS@meta.data$celltype)/sum(table(data_MS@meta.data$celltype))

fre_df <- data.frame(matrix(c(fre_Young, fre_Old, fre_MET, fre_NR, fre_DQ, fre_SPD, fre_MN, fre_MS), nrow  = 8, byrow = T), row.names = c("Young","Old","MET","NR","D+Q","SPD","MET+NR","MET+SPD"))
colnames(fre_df) <- names(fre_Young)

table <- data.frame(sample = rep(rownames(fre_df), each = ncol(fre_df)), cell = rep(colnames(fre_df), times = nrow(fre_df)), fre = c(t(fre_df)))
table$sample <- factor(table$sample, levels = unique(table$sample))
table$cell <- factor(table$cell, levels =  cluster_order)


p <- ggplot(table[table$sample %in% c("Young", "Old", "MET", "NR", "D+Q", "SPD"), ], aes(x = sample, y = 100 * fre, fill = cell, stratum = cell, alluvium = cell))
p <- p + scale_fill_manual(values = color_brain)


p1 <- p +
  geom_bar(position = "fill",
           stat="identity",
           alpha = 1, width = 0.5) + theme(panel.grid = element_blank(),
                                           panel.background = element_rect(color = 'black', fill = 'transparent', size = 1)) +
  geom_col(position = 'stack', width = 0.6) +
  scale_y_continuous(expand=c(0.01, 0.2)) +
  labs(fill = "Cell Types", x = "Samples", y = "Frequence(%)") +
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0)

p1


### MET+NR
p <- ggplot(table[table$sample %in% c("Young", "Old", "MET", "NR", "MET+NR"), ], aes(x = sample, y = 100*fre, fill = cell, stratum = cell, alluvium = cell))
p <- p + scale_fill_manual(values = color_brain)

p1 <- p +
  geom_bar(position = "fill",
           stat="identity",
           alpha = 1, width = 0.5) + theme(panel.grid = element_blank(),
                                           panel.background = element_rect(color = 'black', fill = 'transparent', size = 1)) +
  geom_col(position = 'stack', width = 0.6) +
 scale_y_continuous(expand=c(0.01, 0.2)) +
  labs(fill = "Cell Types", x = "Samples", y = "Frequence(%)") +
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0)

p1


### MET+SPD
p <- ggplot(table[table$sample %in% c("Young", "Old", "MET", "SPD", "MET+SPD"), ], aes(x = sample, y = 100*fre, fill = cell, stratum = cell, alluvium = cell))
p <- p + scale_fill_manual(values = color_brain)

p1 <- p +
  geom_bar(position = "fill",
           stat="identity",
           alpha = 1, width = 0.5) + theme(panel.grid = element_blank(),
                                           panel.background = element_rect(color = 'black', fill = 'transparent', size = 1)) +
  geom_col(position = 'stack', width = 0.6) +
 scale_y_continuous(expand=c(0.01, 0.2)) +
  labs(fill = "Cell Types", x = "Samples", y = "Frequence(%)") +
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0)

p1
