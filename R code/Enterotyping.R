install.packages("here")
install.packages("readxl")
install.packages("dplyr")
install.packages("tibble")
install.packages("writexl")
install.packages("Matrix")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("caret")
install.packages("factoextra")
install.packages("cluster")
install.packages("philentropy")
install.packages("ade4")
install.packages("clusterSim")
install.packages("fpc")
install.packages("tidyverse")
install.packages("forcats")

library(here)
library(readxl)
library(dplyr)
library(tibble)
library(writexl)
library(Matrix)
library(tidyr)
library(ggplot2)
library(caret)
library(factoextra)
library(cluster)
library(philentropy)
library(ade4)
library(clusterSim)
library(fpc)
library(tidyverse)
library(forcats)

packageVersion("cluster")
packageVersion("philentropy")
packageVersion("ade4")
packageVersion("clusterSim")
packageVersion("fpc")



renv::snapshot()

# Load data
Genus_raw <- read_xlsx(here::here("Raw data/02_MIC_Genus_Counts_corrected696.xlsx")) %>%
  as.data.frame()

meta_cols <- c("Pig_ID", "Sample_ID", "Study", "segment", "Batch")
genus_cols <- setdiff(names(Genus_raw), meta_cols)
unclassified_cols <- grep("^Unclassified_", genus_cols, value = TRUE)

print(unclassified_cols)

# Merge unclassified columns
Genus_counts2 <- Genus_raw
if (length(unclassified_cols) > 0) {
  Genus_counts2$Unclassified <- rowSums(Genus_counts2[, unclassified_cols, drop = FALSE], na.rm = TRUE)
  Genus_counts2 <- Genus_counts2[, !names(Genus_counts2) %in% unclassified_cols, drop = FALSE]
}

# Relative abundance
genus_cols2 <- setdiff(names(Genus_counts2), meta_cols)
row_totals <- rowSums(Genus_counts2[, genus_cols2, drop = FALSE], na.rm = TRUE)
print(sum(row_totals == 0))

Genus_rel <- Genus_counts2
Genus_rel[, genus_cols2] <- Genus_rel[, genus_cols2, drop = FALSE] / row_totals

genus_rel_mat <- as.matrix(Genus_rel[, genus_cols2, drop = FALSE])

# Filtering
keep_prev  <- colSums(genus_rel_mat > 0) >= ceiling(0.10 * nrow(genus_rel_mat))
keep_abund <- apply(genus_rel_mat, 2, max, na.rm = TRUE) >= 0.001

genus_rel_mat <- genus_rel_mat[, keep_prev & keep_abund, drop = FALSE]
genus_rel_mat <- genus_rel_mat[, !grepl("(?i)^unclassified", colnames(genus_rel_mat)), drop = FALSE]

# Re-standardize
genus_rel_mat <- sweep(genus_rel_mat, 1, rowSums(genus_rel_mat), "/")

# Jensen-Shannon distance matrix
jsd_dist <- philentropy::distance(genus_rel_mat, method = "jensen-shannon")
jsd_dist <- as.dist(jsd_dist)

# Determine optimal number of clusters using Calinski-Harabasz index and Silhouette width
k_range <- 2:6

cluster_eval <- lapply(k_range, function(k) {
  pam_fit <- pam(jsd_dist, k = k)
  
  ch <- clusterSim::index.G1(
    x = as.matrix(genus_rel_mat),
    cl = pam_fit$clustering,
    d = as.matrix(jsd_dist),
    centrotypes = "medoids"
  )
  
  sil <- pam_fit$silinfo$avg.width
  
  data.frame(
    k = k,
    CH_index = ch,
    Silhouette = sil
  )
}) %>%
  dplyr::bind_rows()

print(cluster_eval)

# Choose optimal k 

best_k_ch  <- cluster_eval$k[which.max(cluster_eval$CH_index)]
best_k_sil <- cluster_eval$k[which.max(cluster_eval$Silhouette)]

best_k_ch
best_k_sil

best_k <- cluster_eval$k[which.max(cluster_eval$Silhouette)]
best_k

# Visualize cluster evaluation metrics
png(
  filename = file.path("03_cluster_evaluation_metrics.png"),
  width = 7,
  height = 3.5,
  units = "in",
  res = 600
)

par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

# CH index
plot(
  cluster_eval$k, cluster_eval$CH_index,
  type = "b", pch = 19,
  xlab = "Number of clusters (k)",
  ylab = "Calinski-Harabasz index",
  main = "CH index"
)
abline(v = best_k_ch, lty = 2, col = "red")
points(best_k_ch, max(cluster_eval$CH_index), pch = 19, cex = 1.3, col = "red")

# Silhouette width
plot(
  cluster_eval$k, cluster_eval$Silhouette,
  type = "b", pch = 19,
  xlab = "Number of clusters (k)",
  ylab = "Average silhouette width",
  main = "Silhouette width"
)
abline(v = best_k_sil, lty = 2, col = "red")
points(best_k_sil, max(cluster_eval$Silhouette), pch = 19, cex = 1.3, col = "red")

par(mfrow = c(1, 1))
dev.off()


# Perform PAM clustering with the optimal number of clusters
pam_res <- pam(jsd_dist, k = best_k)

# jaccard stability by bootstrap
set.seed(1234)

cb <- clusterboot(
  genus_rel_mat,
  B = 100,
  bootmethod = "boot",
  clustermethod = pamkCBI,
  krange = 2,
  seed = 1234
)

cb$bootmean


min_jaccard <- min(cb$bootmean, na.rm = TRUE)
min_jaccard

# Add cluster assignments to the original data frame
cluster_result <- Genus_raw %>%
  mutate(cluster = paste0("Cluster_", pam_res$clustering))

table(cluster_result$cluster)
write_xlsx(cluster_result, here::here("02_Cluster assignments.xlsx"))

# Test association between cluster and Study
tab_study_cluster <- table(cluster_result$Study, cluster_result$cluster)
tab_study_cluster
chisq.test(tab_study_cluster)
round(prop.table(tab_study_cluster, margin = 1) * 100, 1)

chisq_res <- chisq.test(tab_study_cluster)
chisq_res$stdres

# Test association between cluster and Treatment
sample_meta <- read_xlsx(here::here("Raw data/01_MIC_Sample data for batch.xlsx")) %>%
  as.data.frame()

sample_meta2 <- sample_meta %>%
  dplyr::select(Sample_ID, Treatment)


cluster_treatment <- cluster_result %>%
  dplyr::select(Sample_ID, cluster) %>%
  dplyr::left_join(sample_meta2, by = "Sample_ID")

sum(is.na(cluster_treatment$Treatment))
cluster_treatment %>% dplyr::filter(is.na(Treatment))


tab_treat_cluster <- table(cluster_treatment$Treatment, cluster_treatment$cluster)
tab_treat_cluster
fisher.test(tab_treat_cluster, simulate.p.value = TRUE, B = 10000)

# Perform PCoA for visualization
pcoa_res <- ade4::dudi.pco(jsd_dist, scannf = FALSE, nf = 2)

pcoa_df <- data.frame(
  Sample_ID = Genus_raw$Sample_ID,
  Axis1 = pcoa_res$li[, 1],
  Axis2 = pcoa_res$li[, 2],
  cluster = cluster_result$cluster
)

var_explained <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)

# Compute cluster centers
centers_df <- pcoa_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    Axis1 = mean(Axis1, na.rm = TRUE),
    Axis2 = mean(Axis2, na.rm = TRUE),
    .groups = "drop"
  )

# Join cluster centers back to each sample
pcoa_df2 <- pcoa_df %>%
  dplyr::left_join(
    centers_df %>%
      dplyr::rename(center_Axis1 = Axis1, center_Axis2 = Axis2),
    by = "cluster"
  )

png(
  filename = file.path("02_PCoA_cluster_plot.png"),
  width = 6,
  height = 5,
  units = "in",
  res = 600
)

ggplot(pcoa_df2, aes(x = Axis1, y = Axis2, color = cluster)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.8) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.8) +
  geom_segment(
    aes(
      x = center_Axis1, y = center_Axis2,
      xend = Axis1, yend = Axis2
    ),
    linewidth = 0.7,
    alpha = 0.55
  ) +
  geom_point(size = 2.6, alpha = 0.9) +
  scale_color_manual(
    values = c("#72af70dc", "#b168a4"),
    name = "Cluster"
  ) +
  xlab(paste0("PCoA1 (", var_explained[1], "%)")) +
  ylab(paste0("PCoA2 (", var_explained[2], "%)")) +
  coord_equal() +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    axis.line = element_blank()
  ) 

dev.off()

# Identify top genera in each cluster
Genus_rel2 <- Genus_rel
Genus_rel2$cluster <- paste0("Cluster_", pam_res$clustering)
Genus_rel2 <- Genus_rel2[, !grepl("(?i)^unclassified", colnames(Genus_rel2)), drop = FALSE]

cluster_mean <- Genus_rel2 %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    across(-c(Pig_ID, Sample_ID, Study, segment, Batch), mean, na.rm = TRUE)
  )
top_genus_by_cluster_no_unclassified <- cluster_mean %>%
  tidyr::pivot_longer(-cluster, names_to = "Genus", values_to = "Mean_RA") %>%
  dplyr::mutate(Mean_RA = Mean_RA * 100) %>%  
  dplyr::filter(!grepl("(?i)^unclassified", Genus)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(dplyr::desc(Mean_RA), .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

print(top_genus_by_cluster_no_unclassified, n = Inf)

plot_comp <- top_genus_by_cluster_no_unclassified %>%
  mutate(
    cluster = factor(cluster, levels = c("Cluster_2", "Cluster_1")),
    signed_RA = if_else(cluster == "Cluster_2", -Mean_RA, Mean_RA),
    cluster_label = recode(cluster,
                           "Cluster_1" = "LacH",
                           "Cluster_2" = "LacL")
  ) %>%
  group_by(Genus) %>%
  mutate(max_abs = max(abs(signed_RA))) %>%
  ungroup() %>%
  arrange(signed_RA) %>%
  mutate(Genus = fct_reorder(Genus, signed_RA))

cluster_cols <- c("LacL" = "#b168a4", "LacH" = "#72af70")

p_comp <- ggplot(plot_comp, aes(x = signed_RA, y = Genus, fill = cluster_label)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  scale_fill_manual(values = cluster_cols) +
  scale_x_continuous(labels = \(x) abs(x)) +
  labs(
    x = "Mean relative abundance, %",
    y = NULL,
    fill = "Cluster"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black")
  )

p_comp


ggsave(
  filename = here::here("02_top_genus_by_cluster_mirrored_barplot.png"),
  plot = p_comp,
  width = 7.5,
  height = 6,
  dpi = 600
)
diff_df <- data.frame(
  Genus = colnames(cluster_mean[, -1]),
  Diff = as.numeric(cluster_mean[cluster_mean$cluster == "Cluster_2", -1] -
                      cluster_mean[cluster_mean$cluster == "Cluster_1", -1])
)

diff_df <- diff_df[order(abs(diff_df$Diff), decreasing = TRUE), ]

head(diff_df, 20)

# PCA + BCA for taxon weighting
pca_res <- ade4::dudi.pca(genus_rel_mat, center = TRUE, scale = FALSE, scannf = FALSE, nf = 5)

bca_res <- ade4::bca(pca_res, fac = as.factor(pam_res$clustering), scannf = FALSE, nf = 2)

# Calculate taxon weights based on the first BCA axis
taxon_weight <- data.frame(
  Genus = rownames(bca_res$co),
  Axis1 = bca_res$co[, 1],
  TaxonWeight = abs(bca_res$co[, 1])
)

taxon_weight <- taxon_weight[order(taxon_weight$TaxonWeight, decreasing = TRUE), ]

head(taxon_weight, 20)

taxon_weight_df <- data.frame(
  Genus = rownames(bca_res$co),
  Axis1 = bca_res$co[, 1]
) %>%
  filter(!grepl("(?i)^unclassified", Genus)) %>%
  mutate(
    TaxonWeight = abs(Axis1),
    Genus = str_replace_all(Genus, "\\.", " ")
  ) %>%
  slice_max(TaxonWeight, n = 20, with_ties = FALSE) %>%
  arrange(TaxonWeight) %>%
  mutate(Genus = fct_reorder(Genus, TaxonWeight))

p_comp1 <- ggplot(taxon_weight_df, aes(x = TaxonWeight, y = Genus)) +
  geom_col(width = 0.72, fill = "#4C6491") +
  labs(
    x = "Taxon weight",
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
  )

p_comp1

ggsave(
  filename = here::here("02_top_genus_by_cluster_taxon_weight_mirrored_barplot.png"),
  plot = p_comp1,
  width = 7.5,
  height = 6,
  dpi = 600
)

