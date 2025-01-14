library(Seurat)
library(sceasy)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

### Brain
load("brain.harmony_singlet_round2.RData")

brain_O_Y <- subset(brain.harmony, subset = origin %in% c("Young", "Old"))
brain_O_Y@meta.data$origin <- as.character(brain_O_Y@meta.data$origin)
brain_O_Y@meta.data$origin <- factor(brain_O_Y@meta.data$origin, levels = c("Young", "Old"))

Idents(brain_O_Y) <- "origin"

### aging-related DEGs
### DEGs for each celltype
all_markers_O_Y <- c()
for(i in levels(brain_O_Y@meta.data$celltype)){
  brain_O_Y_i <- subset(brain_O_Y, subset = celltype == i)
  Idents(brain_O_Y_i) <- "origin"
  markers_O_Y <- FindAllMarkers(brain_O_Y_i)
  markers_O_Y$celltype <- i
  all_markers_O_Y <- rbind(all_markers_O_Y, markers_O_Y)
}
brain_O_Y_aging <- subset(brain_O_Y, features = unique(all_markers_O_Y$gene))


brain_aging <- subset(brain.harmony, features = unique(all_markers_O_Y$gene))
convertFormat(brain_aging, from="seurat", to="anndata", outFile='brain_aging.h5ad')


### Heart
load("heart.harmony_singlet_round2.RData")

heart_O_Y <- subset(heart.harmony, subset = origin %in% c("Young", "Old"))
heart_O_Y@meta.data$origin <- as.character(heart_O_Y@meta.data$origin)
heart_O_Y@meta.data$origin <- factor(heart_O_Y@meta.data$origin, levels = c("Young", "Old"))

Idents(heart_O_Y) <- "origin"

### aging-related DEGs
### DEGs for each celltype
all_markers_O_Y <- c()
for(i in levels(heart_O_Y@meta.data$celltype)){
  heart_O_Y_i <- subset(heart_O_Y, subset = celltype == i)
  Idents(heart_O_Y_i) <- "origin"
  markers_O_Y <- FindAllMarkers(heart_O_Y_i)
  markers_O_Y$celltype <- i
  all_markers_O_Y <- rbind(all_markers_O_Y, markers_O_Y)
}
heart_O_Y_aging <- subset(heart_O_Y, features = unique(all_markers_O_Y$gene))


heart_aging <- subset(heart.harmony, features = unique(all_markers_O_Y$gene))
convertFormat(heart_aging, from="seurat", to="anndata", outFile='heart_aging.h5ad')


### Kidney
load("kidney.harmony_singlet_round2.RData")

kidney_O_Y <- subset(kidney.harmony, subset = origin %in% c("Young", "Old"))
kidney_O_Y@meta.data$origin <- as.character(kidney_O_Y@meta.data$origin)
kidney_O_Y@meta.data$origin <- factor(kidney_O_Y@meta.data$origin, levels = c("Young", "Old"))

Idents(kidney_O_Y) <- "origin"


### aging-related DEGs
### DEGs for each celltype
all_markers_O_Y <- c()
for(i in levels(kidney_O_Y@meta.data$celltype)){
  kidney_O_Y_i <- subset(kidney_O_Y, subset = celltype == i)
  Idents(kidney_O_Y_i) <- "origin"
  markers_O_Y <- FindAllMarkers(kidney_O_Y_i)
  markers_O_Y$celltype <- i
  all_markers_O_Y <- rbind(all_markers_O_Y, markers_O_Y)
}
kidney_O_Y_aging <- subset(kidney_O_Y, features = unique(all_markers_O_Y$gene))


kidney_aging <- subset(kidney.harmony, features = unique(all_markers_O_Y$gene))
convertFormat(kidney_aging, from="seurat", to="anndata", outFile='kidney_aging.h5ad')


### Liver
load("liver.harmony_singlet_round2.RData")

liver_O_Y <- subset(liver.harmony, subset = origin %in% c("Young", "Old"))
liver_O_Y@meta.data$origin <- as.character(liver_O_Y@meta.data$origin)
liver_O_Y@meta.data$origin <- factor(liver_O_Y@meta.data$origin, levels = c("Young", "Old"))

Idents(liver_O_Y) <- "origin"

### aging-related DEGs
### DEGs for each celltype
all_markers_O_Y <- c()
for(i in levels(liver_O_Y@meta.data$celltype)){
  liver_O_Y_i <- subset(liver_O_Y, subset = celltype == i)
  Idents(liver_O_Y_i) <- "origin"
  markers_O_Y <- FindAllMarkers(liver_O_Y_i)
  markers_O_Y$celltype <- i
  all_markers_O_Y <- rbind(all_markers_O_Y, markers_O_Y)
}
liver_O_Y_aging <- subset(liver_O_Y, features = unique(all_markers_O_Y$gene))


liver_aging <- subset(liver.harmony, features = unique(all_markers_O_Y$gene))
convertFormat(liver_aging, from="seurat", to="anndata", outFile='liver_aging.h5ad')




### result
color_nmn <- c("#46ACC8", "#E2D200", "#A6CEE3", "#B40F20", "#FB9A99",
               "#984EA3", "#F1B412", "#0066B2")
names(color_nmn) <- c("Old", "Young", "MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")

### brain
brain_n_o_y <- subset(brain.harmony, subset = origin %in% c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD"))

### bys
result_brain <- read.csv("prediction_brain_aging_bys.csv")
all(result_brain$cell_name == rownames(brain_n_o_y@meta.data))
brain_n_o_y@meta.data$log_label <- result_brain$y_pred
x <- sort((as.matrix(table(result_brain$y_pred, brain_n_o_y@meta.data$origin))[1, 3:8]) / colSums(as.matrix(table(result_brain$y_pred, brain_n_o_y@meta.data$origin))[, 3:8]))
x <- data.frame(origin = names(x), Ratio = x)
x$origin <- factor(x$origin, levels = c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD"))

celltype_pre_age <- brain_n_o_y@meta.data %>% group_by(celltype) %>% reframe(Ratio = (as.matrix(table(log_label, origin))[1, 3:8]) / colSums(as.matrix(table(log_label, origin))[, 3:8])) %>% as.data.frame()
celltype_pre_age$origin <- c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
celltype_pre_age <- celltype_pre_age %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_pre_age) <- celltype_pre_age[, 1]
celltype_pre_age <- celltype_pre_age[, -1]

pheatmap(celltype_pre_age,
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         scale = "none",
         breaks = seq(from = 0, to = 1, length.out = 100),
         filename = "brain_ratio_by_cell_type_heatmap_all.pdf",
         width = 6, height = 3
)


### barplot
ggplot(data = x %>% filter(origin %in% c("MET", "NR", "D+Q", "SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Brain")
ggsave("brian_ratio.pdf", width = 4, height = 4)
ggsave("brian_ratio.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "NR", "MET+NR")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Brain")
ggsave("brian_ratio_MN.pdf", width = 4, height = 4)
ggsave("brian_ratio_MN.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "SPD", "MET+SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Brain")
ggsave("brian_ratio_MS.pdf", width = 4, height = 4)
ggsave("brian_ratio_MS.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "NR", "SPD", "MET+NR", "MET+SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Brain")
ggsave("brian_ratio_MN_MS.pdf", width = 4, height = 4)
ggsave("brian_ratio_MN_MS.png", width = 4, height = 4)


### heart
heart_n_o_y <- subset(heart.harmony, subset = origin %in% c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD"))

### bys
result_heart <- read.csv("prediction_heart_aging_bys.csv")
all(result_heart$cell_name == rownames(heart_n_o_y@meta.data))
heart_n_o_y@meta.data$log_label <- result_heart$y_pred
x <- sort((as.matrix(table(result_heart$y_pred, heart_n_o_y@meta.data$origin))[1, 3:8]) / colSums(as.matrix(table(result_heart$y_pred, heart_n_o_y@meta.data$origin))[, 3:8]))
x <- data.frame(origin = names(x), Ratio = x)
x$origin <- factor(x$origin, levels = c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD"))

celltype_pre_age <- heart_n_o_y@meta.data %>% group_by(celltype) %>% reframe(Ratio = (as.matrix(table(log_label, origin))[1, 3:8]) / colSums(as.matrix(table(log_label, origin))[, 3:8])) %>% as.data.frame()
celltype_pre_age$origin <- c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
celltype_pre_age <- celltype_pre_age %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_pre_age) <- celltype_pre_age[, 1]
celltype_pre_age <- celltype_pre_age[, -1]

pheatmap(celltype_pre_age,
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         scale = "none",
         breaks = seq(from = 0, to = 1, length.out = 100),
         filename = "heart_ratio_by_cell_type_heatmap_all.pdf",
         width = 5, height = 3
)
pheatmap(celltype_pre_age[c(1, 2, 5), ],
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         filename = "heart_ratio_by_cell_type_heatmap_MN.pdf",
         width = 5, height = 2
)
pheatmap(celltype_pre_age[c(1, 4, 6), ],
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         filename = "heart_ratio_by_cell_type_heatmap_MS.pdf",
         width = 5, height = 2
)
pheatmap(celltype_pre_age[-3, ],
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         filename = "heart_ratio_by_cell_type_heatmap_MN_MS.pdf",
         width = 5, height = 2
)

### barplot
ggplot(data = x %>% filter(origin %in% c("MET", "NR", "D+Q", "SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Heart")
ggsave("heart_ratio.pdf", width = 4, height = 4)
ggsave("heart_ratio.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "NR", "MET+NR")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Heart")
ggsave("heart_ratio_MN.pdf", width = 4, height = 4)
ggsave("heart_ratio_MN.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "SPD", "MET+SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Heart")
ggsave("heart_ratio_MS.pdf", width = 4, height = 4)
ggsave("heart_ratio_MS.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "NR", "SPD", "MET+NR", "MET+SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Heart")
ggsave("heart_ratio_MN_MS.pdf", width = 4, height = 4)
ggsave("heart_ratio_MN_MS.png", width = 4, height = 4)



### kidney
kidney_n_o_y <- subset(kidney.harmony, subset = origin %in% c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD"))

### BYS
result_kidney <- read.csv("prediction_kidney_aging_bys.csv")
all(result_kidney$cell_name == rownames(kidney_n_o_y@meta.data))
kidney_n_o_y@meta.data$log_label <- result_kidney$y_pred
x <- sort((as.matrix(table(result_kidney$y_pred, kidney_n_o_y@meta.data$origin))[1, 3:8]) / colSums(as.matrix(table(result_kidney$y_pred, kidney_n_o_y@meta.data$origin))[, 3:8]))
x <- data.frame(origin = names(x), Ratio = x)
x$origin <- factor(x$origin, levels = c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD"))

celltype_pre_age <- kidney_n_o_y@meta.data %>% group_by(celltype) %>% reframe(Ratio = (as.matrix(table(log_label, origin))[1, 3:8]) / colSums(as.matrix(table(log_label, origin))[, 3:8])) %>% as.data.frame()
celltype_pre_age$origin <- c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
celltype_pre_age <- celltype_pre_age %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_pre_age) <- celltype_pre_age[, 1]
celltype_pre_age <- celltype_pre_age[, -1]

pheatmap(celltype_pre_age,
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         breaks = seq(from = 0, to = 1, length.out = 100),
         names = "Ratio",
         filename = "kidney_ratio_by_cell_type_heatmap_all.pdf",
         width = 5, height = 3
)
pheatmap(celltype_pre_age[c(1, 2, 5), ],
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         filename = "kidney_ratio_by_cell_type_heatmap_MN.pdf",
         width = 5, height = 2
)
pheatmap(celltype_pre_age[c(1, 4, 6), ],
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         filename = "kidney_ratio_by_cell_type_heatmap_MS.pdf",
         width = 5, height = 2
)
pheatmap(celltype_pre_age[-3, ],
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         filename = "kidney_ratio_by_cell_type_heatmap_MN_MS.pdf",
         width = 5, height = 2
)

### barplot
ggplot(data = x %>% filter(origin %in% c("MET", "NR", "D+Q", "SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Kidney")
ggsave("kidney_ratio.pdf", width = 4, height = 4)
ggsave("kidney_ratio.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "NR", "MET+NR")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Kidney")
ggsave("kidney_ratio_MN.pdf", width = 4, height = 4)
ggsave("kidney_ratio_MN.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "SPD", "MET+SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Kidney")
ggsave("kidney_ratio_MS.pdf", width = 4, height = 4)
ggsave("kidney_ratio_MS.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "NR", "SPD", "MET+NR", "MET+SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Kidney")
ggsave("kidney_ratio_MN_MS.pdf", width = 4, height = 4)
ggsave("kidney_ratio_MN_MS.png", width = 4, height = 4)


### liver
liver_n_o_y <- subset(liver.harmony, subset = origin %in% c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD"))

### BYS
result_liver <- read.csv("prediction_liver_aging_bys.csv")
all(result_liver$cell_name == rownames(liver_n_o_y@meta.data))
liver_n_o_y@meta.data$log_label <- result_liver$y_pred
x <- sort((as.matrix(table(result_liver$y_pred, liver_n_o_y@meta.data$origin))[1, 3:8]) / colSums(as.matrix(table(result_liver$y_pred, liver_n_o_y@meta.data$origin))[, 3:8]))
x <- data.frame(origin = names(x), Ratio = x)
x$origin <- factor(x$origin, levels = c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD"))

celltype_pre_age <- liver_n_o_y@meta.data %>% group_by(celltype) %>% reframe(Ratio = (as.matrix(table(log_label, origin))[1, 3:8]) / colSums(as.matrix(table(log_label, origin))[, 3:8])) %>% as.data.frame()
celltype_pre_age$origin <- c("MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
celltype_pre_age <- celltype_pre_age %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_pre_age) <- celltype_pre_age[, 1]
celltype_pre_age <- celltype_pre_age[, -1]

pheatmap(celltype_pre_age,
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         breaks = seq(from = 0, to = 1, length.out = 100),
         names = "Ratio",
         filename = "liver_ratio_by_cell_type_heatmap_all.pdf",
         width = 3, height = 3
)
pheatmap(celltype_pre_age[c(1, 2, 5), ],
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         filename = "liver_ratio_by_cell_type_heatmap_MN.pdf",
         width = 3, height = 2
)
pheatmap(celltype_pre_age[c(1, 4, 6), ],
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         filename = "liver_ratio_by_cell_type_heatmap_MS.pdf",
         width = 3, height = 2
)
pheatmap(celltype_pre_age[-3, ],
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         names = "Ratio",
         filename = "liver_ratio_by_cell_type_heatmap_MN_MS.pdf",
         width = 3, height = 2
)
### barplot
ggplot(data = x %>% filter(origin %in% c("MET", "NR", "D+Q", "SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Liver")
ggsave("liver_ratio.pdf", width = 4, height = 4)
ggsave("liver_ratio.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "NR", "MET+NR")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Liver")
ggsave("liver_ratio_MN.pdf", width = 4, height = 4)
ggsave("liver_ratio_MN.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "SPD", "MET+SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Liver")
ggsave("liver_ratio_MS.pdf", width = 4, height = 4)
ggsave("liver_ratio_MS.png", width = 4, height = 4)

ggplot(data = x %>% filter(origin %in% c("MET", "NR", "SPD", "MET+NR", "MET+SPD")),
       aes(x = origin, y = Ratio, fill = origin, color = origin)) +
  geom_point(size = 5) +
  scale_fill_manual(values = color_nmn) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  guides(color = "none", fill = "none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "", y = "Ratio", title = "Liver")
ggsave("liver_ratio_MN_MS.pdf", width = 4, height = 4)
ggsave("liver_ratio_MN_MS.png", width = 4, height = 4)



### per celltype

celltype_x <- brain_n_o_y@meta.data %>% group_by(celltype, origin) %>% reframe(Ratio = table(svm_label)[1] / (table(svm_label)[1] + table(svm_label)[2])) %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_x) <- celltype_x$origin
celltype_x <- celltype_x[, -1]

pheatmap::pheatmap(as.matrix(celltype_x), cluster_rows = F, border_color = NA)

celltype_x <- brain_n_o_y@meta.data %>% group_by(celltype, origin) %>% reframe(Ratio = table(svm_label)[1] / (table(svm_label)[1] + table(svm_label_aging)[2])) %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_x) <- celltype_x$origin
celltype_x <- celltype_x[, -1]

pheatmap::pheatmap(as.matrix(celltype_x), cluster_rows = F, border_color = NA)



celltype_x <- heart_n_o_y@meta.data %>% group_by(celltype, origin) %>% reframe(Ratio = table(svm_label)[1] / (table(svm_label)[1] + table(svm_label)[2])) %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_x) <- celltype_x$origin
celltype_x <- celltype_x[, -1]

pheatmap::pheatmap(as.matrix(celltype_x), cluster_rows = F, border_color = NA)

celltype_x <- heart_n_o_y@meta.data %>% group_by(celltype, origin) %>% reframe(Ratio = table(svm_label)[1] / (table(svm_label)[1] + table(svm_label_aging)[2])) %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_x) <- celltype_x$origin
celltype_x <- celltype_x[, -1]

pheatmap::pheatmap(as.matrix(celltype_x), cluster_rows = F, border_color = NA)


celltype_x <- kidney_n_o_y@meta.data %>% group_by(celltype, origin) %>% reframe(Ratio = table(svm_label)[1] / (table(svm_label)[1] + table(svm_label)[2])) %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_x) <- celltype_x$origin
celltype_x <- celltype_x[, -1]

pheatmap::pheatmap(as.matrix(celltype_x), cluster_rows = F, border_color = NA)

celltype_x <- kidney_n_o_y@meta.data %>% group_by(celltype, origin) %>% reframe(Ratio = table(svm_label)[1] / (table(svm_label)[1] + table(svm_label_aging)[2])) %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_x) <- celltype_x$origin
celltype_x <- celltype_x[, -1]

pheatmap::pheatmap(as.matrix(celltype_x), cluster_rows = F, border_color = NA)



celltype_x <- liver_n_o_y@meta.data %>% group_by(celltype, origin) %>% reframe(Ratio = table(svm_label)[1] / (table(svm_label)[1] + table(svm_label)[2])) %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_x) <- celltype_x$origin
celltype_x <- celltype_x[, -1]

pheatmap::pheatmap(as.matrix(celltype_x), cluster_rows = F, border_color = NA)

celltype_x <- liver_n_o_y@meta.data %>% group_by(celltype, origin) %>% reframe(Ratio = table(svm_label)[1] / (table(svm_label)[1] + table(svm_label_aging)[2])) %>% pivot_wider(names_from = celltype, values_from = Ratio) %>% as.data.frame()
rownames(celltype_x) <- celltype_x$origin
celltype_x <- celltype_x[, -1]

pheatmap::pheatmap(as.matrix(celltype_x), cluster_rows = F, border_color = NA)