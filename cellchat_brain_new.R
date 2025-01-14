library(Seurat)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(reshape2)
library(wesanderson)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)


object.list <- list(Old = cellchat_old,
                    Young = cellchat_young,
                    MET = cellchat_met,
                    NR = cellchat_nr,
                    `D+Q` = cellchat_dq,
                    SPD = cellchat_spd,
                    MN = cellchat_mn,
                    MS = cellchat_ms
)

color_nmn <- c("#46ACC8", "#E2D200", "#A6CEE3", "#B40F20", "#FB9A99",
               "#984EA3", "#F1B412", "#0066B2")
names(color_nmn) <- c("Old", "Young", "MET", "NR", "D+Q", "SPD", "MN", "MS")


p1 <- compareInteractions(cellchat,
                    show.legend = F,
                    group = c("Young", "Old", "MET", "NR", "D+Q", "SPD", "MN", "MS"),
                    color.use = color_nmn,
                    size.text = 12
)

p2 <- compareInteractions(cellchat,
                    show.legend = F,
                    group = c("Young", "Old", "MET", "NR", "D+Q", "SPD", "MN", "MS"),
                    color.use = color_nmn,
                    size.text = 12,
                    measure = "weight"
)
p1 + p2
ggsave("Brain_number_strength_barplot.png", width = 8, height = 3)
ggsave("Brain_number_strength_barplot.pdf", width = 8, height = 3)

p3 <- netVisual_heatmap_nmn(cellchat,
                  cluster.rows = T,
                  cluster.cols = T,
                  # signaling = pathway_select,
                  comparison = c("Old", "MS"),
                  remove.isolate = T,
                  measure = "count",
                  font.size = 10
)

p4 <- netVisual_heatmap_nmn(cellchat,
                  cluster.rows = T,
                  cluster.cols = T,
                  # signaling = pathway_select,
                  comparison = c("Old", "MS"),
                  remove.isolate = T,
                  measure = "weight",
                  font.size = 10
)


pdf("Brain_O_MS_heatmap.pdf", width = 10, height = 5)
p3 + p4
dev.off()


gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", stacked = T, do.stat = TRUE, comparison = c(5, 6))
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", stacked = F, do.stat = TRUE, comparison = c(5, 6))

gg1 + gg2


pathway.union <- union(union(union(object.list[[1]]@netP$pathways,
                       object.list[[2]]@netP$pathways),
                       object.list[[5]]@netP$pathways),
                       object.list[[6]]@netP$pathways)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[1], width = 8, height = 13)
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[2], width = 8, height = 13)
ht3 <- netAnalysis_signalingRole_heatmap(object.list[[5]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[5], width = 8, height = 13)
ht4 <- netAnalysis_signalingRole_heatmap(object.list[[6]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[6], width = 8, height = 13)

draw(ht1 + ht2 + ht3 + ht4, ht_gap = unit(0.5, "cm"))


data <- netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45, max.dataset = 1, remove.isolate = T)
data <- data$data

netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45, max.dataset = 2, remove.isolate = T)
netVisual_bubble(cellchat_old, angle.x = 45)
netVisual_bubble(cellchat_young, angle.x = 45)


### comparison_2 - comparison_1
heatmap_tow_group <- function(comparison_1, comparison_2){
  comparison_1_name <- names(color_nmn)[comparison_1]
  comparison_2_name <- names(color_nmn)[comparison_2]
  o_y <- netVisual_bubble(cellchat, comparison = c(comparison_1, comparison_2), angle.x = 45, remove.isolate = F)
  o_y <- o_y$data
  o_y$prob[which(is.na(o_y$prob))] <- 0
  o_y_lr <- o_y %>% group_by(interaction_name_2, dataset) %>% reframe(sum(prob))
  all_lr <- levels(o_y_lr$interaction_name_2)
  o_y_lr_diff <- c()
  for(i in all_lr){
    o_y_lr_i <- o_y_lr %>% filter(interaction_name_2 == i)
       if(!comparison_2_name %in% o_y_lr_i$dataset){
         o_y_lr_diff <- c(o_y_lr_diff, - o_y_lr_i$`sum(prob)`[which(o_y_lr_i$dataset == comparison_1_name)])
       } else {
          if(!comparison_1_name %in% o_y_lr_i$dataset){
            o_y_lr_diff <- c(o_y_lr_diff,
                             o_y_lr_i$`sum(prob)`[which(o_y_lr_i$dataset == comparison_2_name)])
          } else {
            o_y_lr_diff <- c(o_y_lr_diff,
                             o_y_lr_i$`sum(prob)`[which(o_y_lr_i$dataset == comparison_2_name)] - o_y_lr_i$`sum(prob)`[which(o_y_lr_i$dataset == comparison_1_name)])
          }
       }
  }
o_y_lr_diff <- data.frame(o_y_lr_diff, names = all_lr)
o_y_lr_diff <- o_y_lr_diff %>% arrange(-o_y_lr_diff)

all_cc <- unique(o_y$group.names)
all_cc_matrix <- data.frame(all_cc, row.names = all_cc)
for(i in all_lr){
  o_y_cc_i <- o_y %>% filter(interaction_name_2 == i)

  o_y_cc_diff <- c()
      for(j in all_cc){
         if(j %in% o_y_cc_i$group.names){
           o_y_cc_i_j <- o_y_cc_i %>% filter(group.names == j)
              if(!comparison_2_name %in% o_y_cc_i_j$dataset){
                o_y_cc_diff <- c(o_y_cc_diff, - o_y_cc_i_j$`prob`[which(o_y_cc_i_j$dataset == comparison_1_name)])
              } else {
                 if(!comparison_1_name %in% o_y_cc_i_j$dataset){
                   o_y_cc_diff <- c(o_y_cc_diff,
                                     o_y_cc_i_j$`prob`[which(o_y_cc_i_j$dataset == comparison_2_name)])
                 } else {
                   o_y_cc_diff <- c(o_y_cc_diff,
                                     o_y_cc_i_j$`prob`[which(o_y_cc_i_j$dataset == comparison_2_name)] - o_y_cc_i_j$`prob`[which(o_y_cc_i_j$dataset == comparison_1_name)])
                 }
              }
         } else {
           o_y_cc_diff <- c(o_y_cc_diff, 0)
         }
      }

  o_y_cc_diff <- data.frame(o_y_cc_diff)
  all_cc_matrix <- cbind(all_cc_matrix, o_y_cc_diff)
}
all_cc_matrix <- all_cc_matrix[, -1]
colnames(all_cc_matrix) <- all_lr
rownames(all_cc_matrix) <- all_cc

all_cc_matrix <- all_cc_matrix[, o_y_lr_diff$names]

top_vector <- o_y_lr_diff$o_y_lr_diff
names(top_vector) <- o_y_lr_diff$names
top_anno <- HeatmapAnnotation(score = top_vector,
                              col = list(score = colorRamp2(c(min(o_y_lr_diff$o_y_lr_diff), 0, max(o_y_lr_diff$o_y_lr_diff)), c("blue", "white", "red"))))
hc_o_y <- Heatmap(all_cc_matrix, top_annotation = top_anno)
pdf(paste0("hc_", comparison_1_name, "_", comparison_2_name, ".pdf"), width = 12, height = 60)
draw(hc_o_y)
dev.off()

}

### 1:Old 2:Young
### comparison_2 - comparison_1
### Young vs Old
heatmap_tow_group(2, 1)
### Old vs D+Q
heatmap_tow_group(1, 5)
### Old vs SPD
heatmap_tow_group(1, 6)
### Young vs SPD
heatmap_tow_group(2, 6)



### Find lr-pair(-)
pathway_heatmap_lr <- function(comparison_1, comparison_2, comparison_3, which_lr){
  comparison_1_name <- names(color_nmn)[comparison_1]
  comparison_2_name <- names(color_nmn)[comparison_2]
  comparison_3_name <- names(color_nmn)[comparison_3]

  pathway.show <- o_y$pathway_name[which(o_y$interaction_name_2 == which_lr)[1]]
  ht_list <- list()
  ht_list[[1]] <- netVisual_heatmap(object.list[[comparison_1]],
                  signaling = pathway.show,
                  color.heatmap = "Reds",
                  title.name = paste(pathway.show, "signaling ", comparison_1_name))

  ht_list[[2]] <- netVisual_heatmap(object.list[[comparison_2]],
                  signaling = pathway.show,
                  color.heatmap = "Reds",
                  title.name = paste(pathway.show, "signaling ", comparison_2_name))

  ht_list[[3]] <- netVisual_heatmap(object.list[[comparison_3]],
                  signaling = pathway.show,
                  color.heatmap = "Reds",
                  title.name = paste(pathway.show, "signaling ", comparison_3_name))

  pdf(paste(pathway.show, "signaling_Y_O_SPD.pdf"), width = 12, height = 4)
  ComplexHeatmap::draw(ht_list[[3]] + ht_list[[1]] + ht_list[[2]], ht_gap = unit(0.5, "cm"))
  dev.off()
}


which_lr <- "Cadm3  - Cadm3"
pathway_heatmap_lr(1, 6, 2, which_lr)


netVisual_aggregate(object.list[[comparison_3]],
                      signaling = pathway.show,
                      layout = "chord",
                      signaling.name = paste(pathway.show, "signaling ", comparison_3_name))


### GO function
library(clusterProfiler)
library(org.Mm.eg.db)
comparison_1_name <- names(color_nmn)[comparison_1]
comparison_2_name <- names(color_nmn)[comparison_2]

### comparison_1_name
o_y_1 <- netVisual_bubble(cellchat, comparison = c(comparison_1, comparison_2), angle.x = 45, remove.isolate = T, max.dataset = comparison_1)
o_y_1 <- o_y_1$data %>% filter(dataset == comparison_1_name)

entrez_ids <- bitr(unique(o_y_1$ligand), fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # 选择生物过程（BP）、细胞组分（CC）或分子功能（MF）
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
simplified_go_results1_l <- simplify(go_results)
simplified_go_results1_l <- simplified_go_results1_l@result

entrez_ids <- bitr(unique(o_y_1$receptor), fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # 选择生物过程（BP）、细胞组分（CC）或分子功能（MF）
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
simplified_go_results1_r <- simplify(go_results)
simplified_go_results1_r <- simplified_go_results1_r@result


### comparison_2_name
o_y_2 <- netVisual_bubble(cellchat, comparison = c(comparison_1, comparison_2), angle.x = 45, remove.isolate = T, max.dataset = comparison_2)
o_y_2 <- o_y_2$data %>% filter(dataset == comparison_2_name)

entrez_ids <- bitr(unique(o_y_2$ligand), fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # 选择生物过程（BP）、细胞组分（CC）或分子功能（MF）
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
simplified_go_results2_l <- simplify(go_results)
simplified_go_results2_l <- simplified_go_results2_l@result

entrez_ids <- bitr(unique(o_y_2$receptor), fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # 选择生物过程（BP）、细胞组分（CC）或分子功能（MF）
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
simplified_go_results2_r <- simplify(go_results)
simplified_go_results2_r <- simplified_go_results2_r@result


all_diff_lr <- c()
for(comparison_2 in 2:8){
  o_y_1 <- netVisual_bubble(cellchat, comparison = c(1, comparison_2), angle.x = 45, remove.isolate = F)
  o_y_1 <- o_y_1$data
  diff_lr <- unique(o_y_1$interaction_name_2)
  all_diff_lr <- union(all_diff_lr, diff_lr)
}

all_diff_lr_matrix <- matrix(NA, length(all_diff_lr), 7)
rownames(all_diff_lr_matrix) <- all_diff_lr
colnames(all_diff_lr_matrix) <- names(color_nmn)[2:8]

for(comparison_2 in 2:8){
  comparison_2_name <- names(color_nmn)[comparison_2]
  o_y_1 <- netVisual_bubble(cellchat, comparison = c(1, comparison_2), angle.x = 45, remove.isolate = F)
  o_y_1 <- o_y_1$data
  for(i in all_diff_lr){
    if(i %in% as.character(o_y_1$interaction_name_2)){
      o_y_1_i <- o_y_1 %>% filter(interaction_name_2 == i) %>% group_by(dataset) %>% reframe(sum(prob))
         if(!"Old" %in% o_y_1_i$dataset){
           all_diff_lr_matrix[i, comparison_2_name] <- as.numeric(o_y_1_i[which(o_y_1_i$dataset == comparison_2_name), 2])
         } else {
           if(!comparison_2_name %in% o_y_1_i$dataset){
             all_diff_lr_matrix[i, comparison_2_name] <- -1 * as.numeric(o_y_1_i[which(o_y_1_i$dataset == "Old"), 2])
           } else {
             all_diff_lr_matrix[i, comparison_2_name] <- as.numeric(o_y_1_i[which(o_y_1_i$dataset == comparison_2_name), 2]) - as.numeric(o_y_1_i[which(o_y_1_i$dataset == "Old"), 2])
           }
         }
    }

  }
}
all_diff_lr_matrix[which(is.na(all_diff_lr_matrix), arr.ind = T)] <- 0
pheatmap(all_diff_lr_matrix, breaks = seq(from = -10, to = 10, length.out = 50))

all_diff_lr_df <- as.data.frame(all_diff_lr_matrix)
all_diff_lr_df$score <- apply(all_diff_lr_df, 1, sum)
all_diff_lr_df$freq_all <- apply(all_diff_lr_df, 1, function(x){
  length(which(x[2:7] != 0))
})
all_diff_lr_df$freq_up <- apply(all_diff_lr_df, 1, function(x){
  length(which(x[2:7] > 0))
})
all_diff_lr_df$freq_down <- apply(all_diff_lr_df, 1, function(x){
  length(which(x[2:7] < 0))
})



library(reshape2)
all_diff_lr_df_dan <- all_diff_lr_df[, 1:5]
all_diff_lr_df_dan$score <- apply(all_diff_lr_df_dan, 1, sum)
all_diff_lr_df_dan$freq_all <- apply(all_diff_lr_df_dan, 1, function(x){
  length(which(x[2:5] != 0))
})
all_diff_lr_df_dan$freq_up <- apply(all_diff_lr_df_dan, 1, function(x){
  length(which(x[2:5] > 0))
})
all_diff_lr_df_dan$freq_down <- apply(all_diff_lr_df_dan, 1, function(x){
  length(which(x[2:5] < 0))
})

all_diff_lr_df_shared <- all_diff_lr_df_dan %>% filter(freq_all >= 2) %>% arrange(-freq_all)
pheatmap(t(scale(t(all_diff_lr_df_shared[, 1:5]))), scale = "none", cluster_rows = T, border_color = NA)

all_diff_lr_df_shared <- cbind(t(scale(t(all_diff_lr_df_shared[, 1:5]))), all_diff_lr_df_shared[, 6:9])
all_diff_lr_df_shared <- melt(all_diff_lr_df_shared[, 1:5], variable.name = "origin", value.name = "value")
all_diff_lr_df_shared$LR_pair <- rownames(all_diff_lr_df_dan %>% filter(freq_all >= 2) %>% arrange(-freq_all))
all_diff_lr_df_shared$LR_pair <- factor(all_diff_lr_df_shared$LR_pair, levels = rownames(all_diff_lr_df_dan %>% filter(freq_all >= 2) %>% arrange(-freq_all)))
all_diff_lr_df_shared$origin <- factor(all_diff_lr_df_shared$origin, levels = c("Young", "MET", "NR", "D+Q", "SPD"))

ggplot(all_diff_lr_df_shared, aes(x = LR_pair, y = value, fill = origin)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_nmn) +
  coord_flip() +
  theme_classic()


### Specific
all_diff_lr_df_specific <- all_diff_lr_df %>% filter(freq_all == 1) %>% arrange(-score)
all_diff_lr_df_specific <- melt(all_diff_lr_df_specific[, 2:7], variable.name = "origin", value.name = "value")
all_diff_lr_df_specific$LR_pair <- rownames(all_diff_lr_df %>% filter(freq_all == 1) %>% arrange(-score))
all_diff_lr_df_specific$LR_pair <- factor(all_diff_lr_df_specific$LR_pair, levels = rownames(all_diff_lr_df %>% filter(freq_all == 1) %>% arrange(-score)))
all_diff_lr_df_specific$origin <- factor(all_diff_lr_df_specific$origin, levels = c("MET", "NR", "D+Q", "SPD"))

ggplot(all_diff_lr_df_specific, aes(x = LR_pair, y = value, fill = origin)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_nmn) +
  coord_flip() +
  theme_classic()


all_diff_lr <- c()
for(comparison_2 in 2:8){
  o_y_1 <- netVisual_bubble(cellchat, comparison = c(1, comparison_2), angle.x = 45, remove.isolate = F)
  o_y_1 <- o_y_1$data
  diff_lr <- unique(o_y_1$interaction_name_2)
  all_diff_lr <- union(all_diff_lr, diff_lr)
}

all_diff_lr_matrix <- matrix(NA, length(all_diff_lr), 8)
rownames(all_diff_lr_matrix) <- all_diff_lr
colnames(all_diff_lr_matrix) <- names(color_nmn)

for(comparison_2 in 2:8){
  comparison_2_name <- names(color_nmn)[comparison_2]
  o_y_1 <- netVisual_bubble(cellchat, comparison = c(1, comparison_2), angle.x = 45, remove.isolate = F)
  o_y_1 <- o_y_1$data
  for(i in all_diff_lr){
    if(i %in% as.character(o_y_1$interaction_name_2)){
      o_y_1_i <- o_y_1 %>% filter(interaction_name_2 == i) %>% group_by(dataset) %>% reframe(sum(prob))
         if(!"Old" %in% o_y_1_i$dataset){
           all_diff_lr_matrix[i, comparison_2_name] <- as.numeric(o_y_1_i[which(o_y_1_i$dataset == comparison_2_name), 2])
         } else {
           if(!comparison_2_name %in% o_y_1_i$dataset){
             all_diff_lr_matrix[i, 1] <- as.numeric(o_y_1_i[which(o_y_1_i$dataset == "Old"), 2])
           } else {
             all_diff_lr_matrix[i, comparison_2_name] <- as.numeric(o_y_1_i[which(o_y_1_i$dataset == comparison_2_name), 2])
             all_diff_lr_matrix[i, 1] <- as.numeric(o_y_1_i[which(o_y_1_i$dataset == "Old"), 2])
           }
         }
    }

  }
}


all_diff_lr_matrix[which(is.na(all_diff_lr_matrix), arr.ind = T)] <- 0

all_diff_lr_df_all_dan <- as.data.frame(all_diff_lr_matrix)
all_diff_lr_df_all_dan <- all_diff_lr_df_all_dan[, 1:6]
all_diff_lr_df_all_dan$freq_all <- apply(all_diff_lr_df_all_dan, 1, function(x){
  length(which(x != 0))
})


all_diff_lr_df_all_dan <- all_diff_lr_df_all_dan %>% filter(freq_all > 0)
all_diff_lr_df_specific_dan <- all_diff_lr_df_all_dan %>% filter(freq_all == 1)
all_diff_lr_df_shared_dan <- all_diff_lr_df_all_dan %>% filter(freq_all > 1)

### compare SPD and MET, NR
pdf("Brain_SPD_vs_MET_NR_heatmap.pdf", width = 6, height = 10)
pheatmap(t(scale(t(all_diff_lr_df_all_dan[!(all_diff_lr_df_all_dan$SPD == 0 & all_diff_lr_df_all_dan$NR == 0 & all_diff_lr_df_all_dan$MET == 0), c(3, 4, 6)]))),
         scale = "none",
         cluster_rows = T,
         border_color = NA,
         color = colorRampPalette(c(brewer.pal(6, 'GnBu')))(100)
)
dev.off()

hc <- hclust(dist(t(scale(t(all_diff_lr_df_all_dan[, c(3, 4, 6)]))), method = "euclidean"), method = "complete")
lr_pair <- cutree(hc, k = 2)

lr_pair_1 <- names(lr_pair)[which(lr_pair %in% 1:2)]
lr_pair_1 <- unique(unlist(strsplit(lr_pair_1, split = "  - ")))
entrez_ids <- bitr(lr_pair_1, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
simplified_go_results_1 <- simplify(go_results)
simplified_go_results_1 <- simplified_go_results_1@result %>% arrange(pvalue)


lr_pair_2 <- names(lr_pair)[which(lr_pair == 2)]
lr_pair_2 <- unique(unlist(strsplit(lr_pair_2, split = "  - ")))
entrez_ids <- bitr(lr_pair_2, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
simplified_go_results_2 <- simplify(go_results)
simplified_go_results_2 <- simplified_go_results_2@result %>% arrange(pvalue)




library(clusterProfiler)
library(org.Mm.eg.db)
d_q_compare_count <- apply(all_diff_lr_df_all_dan, 1, function(x){
  sum(x[5] > x[3:6])
})
d_q_compare <- unique(unlist(strsplit(names(which(d_q_compare_count >= 3)), split = "  - ")))

entrez_ids <- bitr(d_q_compare, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
simplified_go_results_d_q <- simplify(go_results)
simplified_go_results_d_q <- simplified_go_results_d_q@result %>% arrange(pvalue)



met_compare_count <- apply(all_diff_lr_df_all_dan, 1, function(x){
  sum(x[3] > x[3:6])
})
met_compare <- unique(unlist(strsplit(names(which(met_compare_count >= 3)), split = "  - ")))

entrez_ids <- bitr(met_compare, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
simplified_go_results_met <- simplify(go_results)
simplified_go_results_met <- simplified_go_results_met@result %>% arrange(pvalue)




nr_compare_count <- apply(all_diff_lr_df_all_dan, 1, function(x){
  sum(x[4] > x[3:6])
})
nr_compare <- unique(unlist(strsplit(names(which(nr_compare_count >= 3)), split = "  - ")))

entrez_ids <- bitr(nr_compare, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # 选择生物过程（BP）、细胞组分（CC）或分子功能（MF）
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
simplified_go_results_nr <- simplify(go_results)
simplified_go_results_nr <- simplified_go_results_nr@result %>% arrange(pvalue)


spd_compare_count <- apply(all_diff_lr_df_all_dan, 1, function(x){
  sum(x[6] > x[1])
})
spd_compare <- unique(unlist(strsplit(names(which(spd_compare_count == 1)), split = "  - ")))

entrez_ids <- bitr(spd_compare, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # 选择生物过程（BP）、细胞组分（CC）或分子功能（MF）
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
simplified_go_results_spd <- simplify(go_results)
simplified_go_results_spd <- simplified_go_results_spd@result %>% arrange(pvalue)



data_count_score <- data.frame(Origin = names(color_nmn[3:6]),
                               Count = c(sum(all_diff_lr_df_all_dan$MET > all_diff_lr_df_all_dan$Old, 3),
                                         sum(all_diff_lr_df_all_dan$NR > all_diff_lr_df_all_dan$Old, 4),
                                         sum(all_diff_lr_df_all_dan$`D+Q` > all_diff_lr_df_all_dan$Old, 5),
                                         sum(all_diff_lr_df_all_dan$SPD > all_diff_lr_df_all_dan$Old, 6),
                                         -sum(all_diff_lr_df_all_dan$MET < all_diff_lr_df_all_dan$Old, 3),
                                         -sum(all_diff_lr_df_all_dan$NR < all_diff_lr_df_all_dan$Old, 4),
                                         -sum(all_diff_lr_df_all_dan$`D+Q` < all_diff_lr_df_all_dan$Old, 5),
                                         -sum(all_diff_lr_df_all_dan$SPD < all_diff_lr_df_all_dan$Old, 6)
                               ),
                               Score = c(sum(all_diff_lr_df_all_dan[all_diff_lr_df_all_dan$MET > all_diff_lr_df_all_dan$Old, 3]),
                                         sum(all_diff_lr_df_all_dan[all_diff_lr_df_all_dan$NR > all_diff_lr_df_all_dan$Old, 4]),
                                         sum(all_diff_lr_df_all_dan[all_diff_lr_df_all_dan$`D+Q` > all_diff_lr_df_all_dan$Old, 5]),
                                         sum(all_diff_lr_df_all_dan[all_diff_lr_df_all_dan$SPD > all_diff_lr_df_all_dan$Old, 6]),
                                         -sum(all_diff_lr_df_all_dan[all_diff_lr_df_all_dan$MET < all_diff_lr_df_all_dan$Old, 3]),
                                         -sum(all_diff_lr_df_all_dan[all_diff_lr_df_all_dan$NR < all_diff_lr_df_all_dan$Old, 4]),
                                         -sum(all_diff_lr_df_all_dan[all_diff_lr_df_all_dan$`D+Q` < all_diff_lr_df_all_dan$Old, 5]),
                                         -sum(all_diff_lr_df_all_dan[all_diff_lr_df_all_dan$SPD < all_diff_lr_df_all_dan$Old, 6])
                               )
)

data_count_score$Origin <- factor(data_count_score$Origin, levels = names(color_nmn[3:6]))
data_count_score <- data.frame(data_count_score)
ggplot(data_count_score, aes(x = Origin, y = Score, fill = Origin)) +
  geom_bar(stat = "identity", show.legend = F, width = 0.75) +
  #geom_line(aes(y = Count * 4), size = 0.9, group = 1, show.legend = F, linewidth = 0.5) +
  geom_point(aes(y = Count * 4), size = 2, show.legend = F) +
  scale_fill_manual(values = color_nmn) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(-200, 400),
    sec.axis = sec_axis(~ . / 4, name = "Number fot L-R pairs") # 反向缩放以匹配第二个y轴
  ) +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11)
  ) +
  theme_classic() +
  labs(x = "", y = "Score") +
  coord_flip()

ggsave("Brain_4_drugs_compare_vs_old_bar_line.png", width = 5, height = 3)
ggsave("Brain_4_drugs_compare_vs_old_bar_line.pdf", width = 5, height = 3)

### pathway_visualization network
load("brain_up_pathway.RData")
### filter pathway(n = 15)
simplified_go_results_d_q <- simplified_go_results_d_q %>% filter(Count > 2)
simplified_go_results_d_q <- simplified_go_results_d_q[1:15, ]
simplified_go_results_d_q$geneID <- gsub("\\/", ",", simplified_go_results_d_q$geneID)
simplified_go_results_d_q <- data.frame(simplified_go_results_d_q[, c(1, 2, 5, 6)],
                                        Phenotype = 1,
                                        simplified_go_results_d_q[, 8])

simplified_go_results_met <- simplified_go_results_met %>% filter(Count > 2)
simplified_go_results_met <- simplified_go_results_met %>%
   filter(Description %in% c(simplified_go_results_met$Description[c(1:13, 37)], "positive regulation of leukocyte migration"))
simplified_go_results_met$geneID <- gsub("\\/", ",", simplified_go_results_met$geneID)
simplified_go_results_met <- data.frame(simplified_go_results_met[, c(1, 2, 5, 6)],
                                        Phenotype = 1,
                                        simplified_go_results_met[, 8])

simplified_go_results_nr <- simplified_go_results_nr %>% filter(Count > 2)
 simplified_go_results_nr <- simplified_go_results_nr[1:15, ]
simplified_go_results_nr$geneID <- gsub("\\/", ",", simplified_go_results_nr$geneID)
simplified_go_results_nr <- data.frame(simplified_go_results_nr[, c(1, 2, 5, 6)],
                                        Phenotype = 1,
                                        simplified_go_results_nr[, 8])


simplified_go_results_spd <- simplified_go_results_spd %>% filter(Count > 2)
simplified_go_results_spd <- simplified_go_results_spd[c(1:14, 31), ]
simplified_go_results_spd$geneID <- gsub("\\/", ",", simplified_go_results_spd$geneID)
simplified_go_results_spd <- data.frame(simplified_go_results_spd[, c(1, 2, 5, 6)],
                                        Phenotype = 1,
                                        simplified_go_results_spd[, 8])

write.table(simplified_go_results_d_q, file = "brain_enrich/enrich_d_q.txt", sep = "\t", quote = F, row.names = F)
write.table(simplified_go_results_met, file = "brain_enrich/enrich_met.txt", sep = "\t", quote = F, row.names = F)
write.table(simplified_go_results_nr, file = "brain_enrich/enrich_nr.txt", sep = "\t", quote = F, row.names = F)
write.table(simplified_go_results_spd, file = "brain_enrich/enrich_spd.txt", sep = "\t", quote = F, row.names = F)


all_diff_lr_df_all_dan <- as.data.frame(all_diff_lr_matrix)
all_diff_lr_df_all_dan <- all_diff_lr_df_all_dan[, c(1:3, 6, 8)]
all_diff_lr_df_all_dan$freq_all <- apply(all_diff_lr_df_all_dan, 1, function(x){
  length(which(x != 0))
})

all_diff_lr_df_all_dan <- all_diff_lr_df_all_dan %>% filter(freq_all > 0)
all_diff_lr_df_specific_dan <- all_diff_lr_df_all_dan %>% filter(freq_all == 1)
all_diff_lr_df_shared_dan <- all_diff_lr_df_all_dan %>% filter(freq_all > 1)

### compare SPD and MET, NR
pdf("Brain_MS_heatmap.pdf", width = 3.5, height = 8)
pheatmap(t(scale(t(all_diff_lr_df_all_dan[!(all_diff_lr_df_all_dan$SPD == 0 & all_diff_lr_df_all_dan$MET == 0 & all_diff_lr_df_all_dan$MS == 0), c(1, 3:5)]))),
         scale = "none",
         cluster_rows = T,
         border_color = NA,
         color = colorRampPalette(c(brewer.pal(6, 'GnBu')))(100)
)
dev.off()


### MN
all_diff_lr_df_all_dan <- as.data.frame(all_diff_lr_matrix)
all_diff_lr_df_all_dan <- all_diff_lr_df_all_dan[, c(1:4, 7)]
all_diff_lr_df_all_dan$freq_all <- apply(all_diff_lr_df_all_dan, 1, function(x){
  length(which(x != 0))
})

all_diff_lr_df_all_dan <- all_diff_lr_df_all_dan %>% filter(freq_all > 0)
all_diff_lr_df_specific_dan <- all_diff_lr_df_all_dan %>% filter(freq_all == 1)
all_diff_lr_df_shared_dan <- all_diff_lr_df_all_dan %>% filter(freq_all > 1)

### compare SPD and MET, NR
pdf("Brain_MN_heatmap.pdf", width = 3.5, height = 8)
pheatmap(t(scale(t(all_diff_lr_df_all_dan[!(all_diff_lr_df_all_dan$NR == 0 & all_diff_lr_df_all_dan$MET == 0 & all_diff_lr_df_all_dan$MN == 0), c(1, 3:5)]))),
         scale = "none",
         cluster_rows = T,
         border_color = NA,
         color = colorRampPalette(c(brewer.pal(6, 'GnBu')))(100)
)
dev.off()

GO_data <- clusterProfiler:::get_GO_data("org.Mm.eg.db", "BP", "SYMBOL")
names(GO_data$PATHID2EXTID) <- GO_data$PATHID2NAME[match(names(GO_data$PATHID2EXTID), names(GO_data$PATHID2NAME))]

pathway_lr <- apply(pathway_type, 1, function(x){
  as.character(unlist(GO_data$PATHID2EXTID[names(GO_data$PATHID2EXTID) == x[1]]))
})
names(pathway_lr) <- pathway_type$Description

### comparsion <- c(1, 2), 1 > 2
pathway_circos <- function(comparsion){
  young_old_lr <- rownames(all_diff_lr_df_all_dan %>% filter(!!sym(comparsion[1]) > !!sym(comparsion[2])))
  young_old_score_df <- data.frame(pathway_type, origin = paste0(comparsion, collapse = '/'), score = 0)
  for(i in 1:length(pathway_type$Description)){
    score_all <- 0
    for(j in young_old_lr){
      l_r_i <- unlist(strsplit(j, split = "  - "))
      ### 如果配受体有一个在通路里，则计数
      if(l_r_i[1] %in% pathway_lr[[i]]){
        score_all <- score_all + all_diff_lr_df_all_dan[j, comparsion[1]]
      }
      if(l_r_i[2] %in% pathway_lr[[i]]){
        score_all <- score_all + all_diff_lr_df_all_dan[j, comparsion[1]]
      }
    }
    young_old_score_df[i, 4] <- score_all
  }
  return(young_old_score_df)
}


young_old_score_df <- pathway_circos(c("Young", "Old"))
young_old_score_df <- young_old_score_df[-c(38, 39), ]

met_old_score_df <- pathway_circos(c("MET", "Old"))
met_old_score_df <- met_old_score_df[-c(38, 39), ]

dq_old_score_df <- pathway_circos(c("D+Q", "Old"))
dq_old_score_df <- dq_old_score_df[-c(38, 39), ]

nr_old_score_df <- pathway_circos(c("NR", "Old"))
nr_old_score_df <- nr_old_score_df[-c(38, 39), ]

spd_old_score_df <- pathway_circos(c("SPD", "Old"))
spd_old_score_df <- spd_old_score_df[-c(38, 39), ]

dq_old_score_df_down <- pathway_circos(c("Old", "D+Q"))
dq_old_score_df_down <- dq_old_score_df_down[c(38, 39), ]

spd_old_score_df_down <- pathway_circos(c("Old", "SPD"))
spd_old_score_df_down <- spd_old_score_df_down[c(38, 39), ]

### combined
pathway_circos_combine <- rbind(young_old_score_df, rbind(met_old_score_df, rbind(dq_old_score_df, rbind(nr_old_score_df, rbind(spd_old_score_df, rbind(dq_old_score_df_down, spd_old_score_df_down))))))

node_info_1 <- pathway_circos_combine %>% group_by(Description) %>% reframe(Size = sum(score), type = type) %>% as.data.frame() %>% unique()
node_info <- pathway_circos_combine %>% group_by(origin) %>% reframe(Size = sum(score)) %>% as.data.frame()
node_info$type <- "Group"
colnames(node_info)[1] <- 'Node'
colnames(node_info_1)[1] <- 'Node'
node_info <- rbind(node_info, node_info_1)

write.table(node_info, file = "pathway_circos_combine_node.txt", sep = "\t", quote = F, row.names = F)
write.table(pathway_circos_combine, file = "pathway_circos_combine.txt", sep = "\t", quote = F, row.names = F)

pathway_circos_combine <- read.table("pathway_circos_combine.txt", sep = "\t", header = T)

library(circlize)
chordDiagram(pathway_circos_combine[, -1], scale = F)


### barplot
bar_info <- pathway_circos_combine %>% group_by(origin, type) %>% reframe(Size = sum(score)) %>% as.data.frame()

ggplot(data = bar_info %>% filter(type == "development"), aes(x = Size, y = origin, fill = Size)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_gradient2(high = '#85B94B', low = "white") +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11)) +
  labs(y = "", fill = "Score")
ggsave("brain_development_bar_plot.pdf", width = 5, height = 2)

ggplot(data = bar_info %>% filter(type == "synapse-related"), aes(x = Size, y = origin, fill = Size)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_gradient2(brewer.pal(8, "Purples")) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11)) +
  labs(y = "", fill = "Score")
ggsave("brain_synapse_bar_plot.pdf", width = 5, height = 2)

ggplot(data = bar_info %>% filter(type == "Immune related"), aes(x = Size, y = origin, fill = Size)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_gradient2(high = '#C568CC', low = "white") +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11)) +
  labs(y = "", fill = "Score")
ggsave("brain_Immune_bar_plot.pdf", width = 5, height = 2)

ggplot(data = bar_info %>% filter(type == "brain function"), aes(x = Size, y = origin, fill = Size)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_gradient2(high = '#E77CAF', low = "white") +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11)) +
  labs(y = "", fill = "Score")
ggsave("brain_function_bar_plot.pdf", width = 5, height = 2)