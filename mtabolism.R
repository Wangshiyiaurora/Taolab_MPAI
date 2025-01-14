library(dplyr)
library(pheatmap)


expression_metabolism <- read.table("metabolites_full_table.txt", sep = "\t", quote = "", comment.char = "", header = T, row.names = 1)

### PLS-DA
library(mixOmics)
expr <- as.data.frame(t(expression_metabolism[, 2:37]))
pos_expr <- as.data.frame(t(expression_metabolism[2080:3795, 2:37]))
neg_expr <- as.data.frame(t(expression_metabolism[1:2079, 2:37]))
class <- rep(c("MET",  "NR", "Old", "D+Q", "SPD", "Young"), each = 6)

plsda.srbct <- plsda(expr, class, ncomp = 10)
plsda.srbct <- plsda(pos_expr, class, ncomp = 10)
plsda.srbct <- plsda(neg_expr, class, ncomp = 10)

set.seed(123)
perf.plsda.srbct <- mixOmics::perf(plsda.srbct, validation = "Mfold", folds = 3, progressBar = F, nrepeat = 10)
plot(perf.plsda.srbct, sd = T, legend.posotion = "horizontal")

final.plsda.srbct <- plsda(expr, class, ncomp = 3)
color_nmn <- c("#46ACC8", "#E2D200", "#A6CEE3", "#B40F20", "#FB9A99",
               "#984EA3", "#F1B412", "#0066B2")
names(color_nmn) <- c("Old", "Young", "MET", "NR", "D+Q", "SPD", "MET+NR", "MET+SPD")
plotIndiv(final.plsda.srbct, ind.names = FALSE, legend = TRUE,
          comp = c(1, 2), ellipse = TRUE,
          title = 'PLS-DA on SRBCT comp 1-2',
          col.per.group = color_nmn[c("D+Q", "MET", "NR", "Old", "SPD", "Young")]
)
ggsave("PLS-DA.pdf", width = 8, height = 6)

library(ropls)
sacurine.pca <- opls(expr)
plot(sacurine.pca, typeVc = "x-score", parAsColFcVn = class)

sacurine.plsda <- opls(expr, class)
plot(sacurine.plsda, typeVc = "x-score", parAsColFcVn = class)


expression_metabolism_o_y_up <- expression_metabolism %>% filter(C_vs_F_VIP > 1 & C_vs_F_log2FC > 0.58 & C_vs_F_Pvalue < 0.05)
expression_metabolism_o_y_down <- expression_metabolism %>% filter(C_vs_F_VIP > 1 & C_vs_F_log2FC < -0.58 & C_vs_F_Pvalue < 0.05)
UP_KEGG <- unlist(apply(expression_metabolism_o_y_up[rownames(fc_up), ], 1, function(x){
  strsplit(x[103], split = ";;")[[1]]
}))
sort(table(as.character(UP_KEGG)))

DOWN_KEGG <- unlist(apply(expression_metabolism_o_y_down[rownames(fc_down), ], 1, function(x){
  strsplit(x[103], split = ";;")[[1]]
}))
sort(table(as.character(DOWN_KEGG)))

fc_all <- expression_metabolism[, c(42, 70, 54, 46, 38)]
fc_all[, 4] <- 1 / fc_all[, 4]
fc_all[, 5] <- 1 / fc_all[, 5]
colnames(fc_all) <- c("Old/Young", "MET/Old", "NR/Old", "D+Q/Old", "SPD/Old")


fc_up <- fc_all[rownames(expression_metabolism %>% filter(`C_vs_F_.FC_1_Pvalue_0.05_VIP_1._regulated` == "up" & C_vs_F_log2FC > 1)), ]
fc_up$name <- expression_metabolism[rownames(fc_up), "name"]
pheatmap::pheatmap(as.matrix(fc_up),
                   scale = "none",
                   cluster_cols = F,
                   border_color = NA,
                   breaks = seq(-4, 4, length.out = 100),
                   filename = "UP_FC_1_heatmap.pdf",
                   width = 6, height = 10
)
fc_down <- fc_all[rownames(expression_metabolism_o_y_down %>% filter(`C_vs_F_.FC_1_Pvalue_0.05_VIP_1._regulated` == "down" & C_vs_F_log2FC < -1)), ]
fc_down$name <- expression_metabolism[rownames(fc_down), "name"]
pheatmap::pheatmap(as.matrix(fc_down),
                   scale = "none",
                   cluster_cols = F,
                   border_color = NA,
                   breaks = seq(-4, 4, length.out = 100),
                   filename = "DOWN_FC_1_heatmap.pdf",
                   width = 6, height = 10
)


up_metabolites_name_pos <- expression_metabolism_o_y_up$name[grep("pos_", rownames(expression_metabolism_o_y_up))]
down_metabolites_name_pos <- expression_metabolism_o_y_down$name[grep("pos_", rownames(expression_metabolism_o_y_down))]
up_metabolites_name_neg <- expression_metabolism_o_y_up$name[grep("neg_", rownames(expression_metabolism_o_y_up))]
down_metabolites_name_neg <- expression_metabolism_o_y_down$name[grep("neg_", rownames(expression_metabolism_o_y_down))]

write.table(up_metabolites_name_pos, file = "up_pos_id.txt", sep = "\n", quote = F, row.names = F, col.names = F)
write.table(down_metabolites_name_pos, file = "down_pos_id.txt", sep = "\n", quote = F, row.names = F, col.names = F)
write.table(up_metabolites_name_neg, file = "up_neg_id.txt", sep = "\n", quote = F, row.names = F, col.names = F)
write.table(down_metabolites_name_neg, file = "down_neg_id.txt", sep = "\n", quote = F, row.names = F, col.names = F)

write.table(expression_metabolism_o_y_up$name, file = "up_id.txt", sep = "\n", quote = F, row.names = F, col.names = F)
write.table(expression_metabolism_o_y_down$name, file = "down_id.txt", sep = "\n", quote = F, row.names = F, col.names = F)


library(qs)
hmdb <- qs::qread("compound_db.qs")
smpdb_pathway <- qs::qread("smpdb_pathway.qs")
syn.db <- qs::qread("syn_nms.qs")
kegg_pathway <- qs::qread("kegg_pathway.qs")

up_yes <- expression_metabolism_o_y_up[which(expression_metabolism_o_y_up$name %in% hmdb$name), "name"]
up_no <- expression_metabolism_o_y_up[which(!expression_metabolism_o_y_up$name %in% hmdb$name), c("KEGG_annotation", "HMDB_ID")]
up_no$remap <- apply(up_no, 1, function(x){
  ifelse(x[2] %in% hmdb$hmdb_id, hmdb$name[which(hmdb$hmdb_id == x[2])],
         ifelse(x[1] %in% hmdb$kegg_id, hmdb$name[which(hmdb$kegg_id == x[1])], NA))
})
up_no <- up_no %>% filter(!is.na(remap))
up_all <- c(up_yes, up_no$remap)
write.table(up_all, file = "up_id_map.txt", sep = "\n", quote = F, row.names = F, col.names = F)


down_yes <- expression_metabolism_o_y_down[which(expression_metabolism_o_y_down$name %in% hmdb$name), "name"]
down_no <- expression_metabolism_o_y_down[which(!expression_metabolism_o_y_down$name %in% hmdb$name), c("KEGG_annotation", "HMDB_ID")]
down_no$remap <- apply(down_no, 1, function(x){
  ifelse(x[2] %in% hmdb$hmdb_id, hmdb$name[which(hmdb$hmdb_id == x[2])],
         ifelse(x[1] %in% hmdb$kegg_id, hmdb$name[which(hmdb$kegg_id == x[1])], NA))
})
down_no <- down_no %>% filter(!is.na(remap))
down_all <- c(down_yes, down_no$remap)
write.table(down_all, file = "down_id_map.txt", sep = "\n", quote = F, row.names = F, col.names = F)

expression_metabolism_o_y_up[which(!(expression_metabolism_o_y_up$KEGG_annotation %in% hmdb$kegg_id | expression_metabolism_o_y_up$HMDB_ID %in% hmdb$hmdb_id)), c("KEGG_annotation", "HMDB_ID")]


library(ggplot2)
library(RColorBrewer)

### UP
up_enrich <- read.csv("pathway_results_up_enrich.csv")
up_enrich <- up_enrich %>% filter(Hits > 1 & Impact > 0)
up_enrich$rich_factor <- up_enrich$Hits / up_enrich$Total

up_msea <- read.csv("msea_ora_result_up_msea.csv")
up_msea <- up_msea %>% filter(hits > 1)

### DOWN
down_enrich <- read.csv("pathway_results_down_enrich.csv")
down_enrich <- down_enrich %>% filter(Hits > 1 & Impact > 0)
down_enrich$rich_factor <- down_enrich$Hits / down_enrich$Total

down_msea <- read.csv("msea_ora_result_down_msea.csv")
down_msea <- down_msea %>% filter(hits > 1)

up_enrich_draw <- up_enrich[1:4, ] %>% arrange(rich_factor)
up_enrich_draw$X <- factor(up_enrich_draw$X, levels = up_enrich_draw$X)
ggplot(data = up_enrich_draw, aes(x = rich_factor, y = X, fill = Impact)) +
  geom_bar(stat = 'identity', width = 0.65) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13)
  ) +
  labs(x = "Rich factor", y = "")

ggsave("O_Y_up_pathway_barplot.pdf", width = 8, height = 3)

down_enrich_draw <- down_enrich[1:4, ] %>% arrange(rich_factor)
down_enrich_draw$X <- factor(down_enrich_draw$X, levels = down_enrich_draw$X)
ggplot(data = down_enrich_draw, aes(x = rich_factor, y = X, fill = Impact)) +
  geom_bar(stat = 'identity', width = 0.65) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13)
  ) +
  labs(x = "Rich factor", y = "")

ggsave("O_Y_down_pathway_barplot.pdf", width = 8, height = 3)

median_intensity <- apply(expression_metabolism, 1, function(x){
  in_Old <- median(as.numeric(x[14:19]))
  in_Y <- median(as.numeric(x[32:37]))

  ### MET
  in_i <- median(as.numeric(x[2:7]))
  met_label <- ifelse(in_i > in_Old & in_Y > in_Old & x[72] > 1 & x[44] > 1, in_i - in_Old,
         ifelse(in_i < in_Old & in_Y < in_Old & x[72] > 1 & x[44] > 1, in_i - in_Old, 0))

  ### NR
  in_i <- median(as.numeric(x[8:13]))
  nr_label <- ifelse(in_i > in_Old & in_Y > in_Old & x[56] > 1 & x[44] > 1, in_i - in_Old,
         ifelse(in_i < in_Old & in_Y < in_Old & x[56] > 1 & x[44] > 1, in_i - in_Old, 0))

  ### D+Q
  in_i <- median(as.numeric(x[20:25]))
  dq_label <- ifelse(in_i > in_Old & in_Y > in_Old & x[48] > 1 & x[44] > 1, in_i - in_Old,
         ifelse(in_i < in_Old & in_Y < in_Old & x[48] > 1 & x[44] > 1, in_i - in_Old, 0))

  ### SPD
  in_i <- median(as.numeric(x[26:31]))
  spd_label <- ifelse(in_i > in_Old & in_Y > in_Old & x[40] > 1 & x[44] > 1, in_i - in_Old,
         ifelse(in_i < in_Old & in_Y < in_Old & x[40] > 1 & x[44] > 1, in_i - in_Old, 0))

  result <- c(met_label, nr_label, dq_label, spd_label)
  return(result)
})

median_intensity <- as.data.frame(t(median_intensity))
colnames(median_intensity) <- c("MET", "NR", "D+Q", "SPD")

apply(median_intensity, 2, function(x){
  length(which(x > 0))
})

apply(median_intensity, 2, function(x){
  length(which(x < 0))
})


median_intensity_fan <- apply(expression_metabolism, 1, function(x){
  in_Old <- median(as.numeric(x[14:19]))
  in_Y <- median(as.numeric(x[32:37]))

  ### MET
  in_i <- median(as.numeric(x[2:7]))
  met_label <- ifelse(in_i > in_Old & in_Y < in_Old & x[44] > 1, "up",
         ifelse(in_i < in_Old & in_Y > in_Old & x[44] > 1, "down", "none"))

  ### NR
  in_i <- median(as.numeric(x[8:13]))
  nr_label <- ifelse(in_i > in_Old & in_Y < in_Old & x[44] > 1, "up",
         ifelse(in_i < in_Old & in_Y > in_Old & x[44] > 1, "down", "none"))

  ### D+Q
  in_i <- median(as.numeric(x[20:25]))
  dq_label <- ifelse(in_i > in_Old & in_Y < in_Old & x[44] > 1, "up",
         ifelse(in_i < in_Old & in_Y > in_Old & x[44] > 1, "down", "none"))


  ### SPD
  in_i <- median(as.numeric(x[26:31]))
  spd_label <- ifelse(in_i > in_Old & in_Y < in_Old & x[44] > 1, "up",
         ifelse(in_i < in_Old & in_Y > in_Old & x[44] > 1, "down", "none"))


  result <- c(met_label, nr_label, dq_label, spd_label)
  return(result)
})

median_intensity_fan <- as.data.frame(t(median_intensity_fan))
colnames(median_intensity_fan) <- c("MET", "NR", "D+Q", "SPD")

apply(median_intensity_fan, 2, function(x){
  length(which(x == "up"))
}) +
apply(median_intensity_fan, 2, function(x){
  length(which(x == "down"))
})


metabolite_par_anti <- data.frame(group = c("MET", "NR", "D+Q", "SPD", "MET", "NR", "D+Q", "SPD"),
                                  count = c(-191, -298, -475, -495, 263+184, 194+200, 190+123, 179+140))
metabolite_par_anti$group <- factor(metabolite_par_anti$group, levels = c("SPD", "D+Q", "NR", "MET"))
ggplot(data = metabolite_par_anti, aes(x = count, y = group, color = group)) +
  geom_segment(aes(x = 0, xend = count, y = group, yend = group),
               linetype = "solid",
               linewidth = 2) +
  geom_point(aes(size = abs(count))) +
  scale_size_continuous(range = c(5, 10)) +
  scale_color_manual(values = color_nmn) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13)
  ) +
  labs(x = "Count", y = "")

ggsave("metabolite_size_bangbangtang.pdf", width = 6, height = 3)

up_met <- expression_metabolism[which(median_intensity$`D+Q` < 0), ]
up_met_yes <- intersect(up_met$name, hmdb$name)
up_met_no <- up_met[which(!up_met$name %in% hmdb$name), c("KEGG_annotation", "HMDB_ID")]
up_met_no$remap <- apply(up_met_no, 1, function(x){
  ifelse(x[2] %in% hmdb$hmdb_id, hmdb$name[which(hmdb$hmdb_id == x[2])],
         ifelse(x[1] %in% hmdb$kegg_id, hmdb$name[which(hmdb$kegg_id == x[1])], NA))
})
up_met_no <- up_met_no %>% filter(!is.na(remap))
up_met_all <- c(up_met_yes, up_met_no$remap)
length(up_met_all)

up_yes_fc <- up_met[which(up_met$name %in% hmdb$name), ]
up_not_fc <- up_met[rownames(up_met_no), ]
up_not_fc$name <- up_met_no$remap

up_all_fc <- rbind(up_yes_fc[, c("name", "C_vs_D_log2FC")], up_not_fc[, c("name", "C_vs_D_log2FC")])
up_all_fc[, 2] <- -1 * up_all_fc[, 2]

write.table(up_met_all, file = "up_id_map_dq.txt", sep = "\n", quote = F, row.names = F, col.names = F)
write.table(up_all_fc, file = "up_all_fc_dq.txt", sep = "\t", quote = F, row.names = F, col.names = F)

write.table(up_met_all, file = "down_id_map_dq.txt", sep = "\n", quote = F, row.names = F, col.names = F)
write.table(up_all_fc, file = "down_all_fc_dq.txt", sep = "\t", quote = F, row.names = F, col.names = F)



met_up_pathway <- read.csv("met_up_enrich.csv")
met_up_pathway <- met_up_pathway %>% filter(Hits > 1 & Impact > 0)
met_up_pathway$rich_factor <- met_up_pathway$Hits / met_up_pathway$Total

met_down_pathway <- read.csv("met_down_enrich.csv")
met_down_pathway <- met_down_pathway %>% filter(Hits > 1 & Impact > 0)
met_down_pathway$rich_factor <- met_down_pathway$Hits / met_down_pathway$Total


nr_up_pathway <- read.csv("nr_up_enrich.csv")
nr_up_pathway <- nr_up_pathway %>% filter(Hits > 1 & Impact > 0)
nr_up_pathway$rich_factor <- nr_up_pathway$Hits / nr_up_pathway$Total

nr_down_pathway <- read.csv("nr_down_enrich.csv")
nr_down_pathway <- nr_down_pathway %>% filter(Hits > 1 & Impact > 0)
nr_down_pathway$rich_factor <- nr_down_pathway$Hits / nr_down_pathway$Total

dq_up_pathway <- read.csv("dq_up_enrich.csv")
dq_up_pathway <- dq_up_pathway %>% filter(Hits > 1 & Impact > 0)
dq_up_pathway$rich_factor <- dq_up_pathway$Hits / dq_up_pathway$Total

dq_down_pathway <- read.csv("dq_down_enrich.csv")
dq_down_pathway <- dq_down_pathway %>% filter(Hits > 1 & Impact > 0)
dq_down_pathway$rich_factor <- dq_down_pathway$Hits / dq_down_pathway$Total


spd_up_pathway <- read.csv("spd_up_enrich.csv")
spd_up_pathway <- spd_up_pathway %>% filter(Hits > 1 & Impact > 0)
spd_up_pathway$rich_factor <- spd_up_pathway$Hits / spd_up_pathway$Total

spd_down_pathway <- read.csv("spd_down_enrich.csv")
spd_down_pathway <- spd_down_pathway %>% filter(Hits > 1 & Impact > 0)
spd_down_pathway$rich_factor <- spd_down_pathway$Hits / spd_down_pathway$Total


pathway_rev_up <- rbind(data.frame(met_up_pathway, group = "MET"),
                        data.frame(nr_up_pathway, group = "NR"),
                        data.frame(dq_up_pathway, group = "D+Q"),
                        data.frame(spd_up_pathway, group = "SPD")
)
pathway_rev_up <- pathway_rev_up %>% filter(!X %in% all_fgsea_result$pathway[which(all_fgsea_result$NES > 0)])

pathway_rev_up$group <- factor(pathway_rev_up$group, levels = c("MET", "NR", "D+Q", "SPD"))
pathway_rev_up <- pathway_rev_up %>% arrange(group, desc(Impact))
pathway_rev_up$X <- factor(pathway_rev_up$X, levels = unique(pathway_rev_up$X))
save(pathway_rev_up, file = "pathway_rev_up.RData")

ggplot(pathway_rev_up, aes(x = group, y = X, color = Impact)) +
  geom_point(aes(size = Hits)) +
  scale_size_continuous(range = c(5, 8)) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  # scale_color_gradient(low = brewer.pal(9, 'Spectral')[4], high = brewer.pal(9, 'Spectral')[1]) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 11, color = "black")) +
  labs(x = "", y = "")
ggsave("rev_up_dot_plot.pdf", width = 6, height = 5.5)

### Down
pathway_rev_down <- rbind(data.frame(met_down_pathway, group = "MET"),
                        data.frame(nr_down_pathway, group = "NR"),
                        data.frame(dq_down_pathway, group = "D+Q"),
                        data.frame(spd_down_pathway, group = "SPD")
)
pathway_rev_down <- pathway_rev_down %>% filter(!X %in% all_fgsea_result$pathway[which(all_fgsea_result$NES < 0)])
pathway_rev_down$group <- factor(pathway_rev_down$group, levels = c("MET", "NR", "D+Q", "SPD"))
pathway_rev_down <- pathway_rev_down %>% arrange(group, desc(Impact))
pathway_rev_down$X <- factor(pathway_rev_down$X, levels = unique(pathway_rev_down$X))
save(pathway_rev_down, file = "pathway_rev_down.RData")

ggplot(pathway_rev_down, aes(x = group, y = X, color = Impact)) +
  geom_point(aes(size = Hits)) +
  scale_size_continuous(range = c(5, 8)) +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  # scale_color_gradient(low = brewer.pal(9, 'YlGnBu')[2], high = brewer.pal(9, 'YlGnBu')[9]) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 11, color = "black")) +
  labs(x = "", y = "")
ggsave("rev_down_dot_plot.pdf", width = 6.8, height = 3)