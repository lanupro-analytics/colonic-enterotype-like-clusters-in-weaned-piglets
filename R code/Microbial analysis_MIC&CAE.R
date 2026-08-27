library(pheatmap)
library(RVAideMemoire)
library(PLSDAbatch)
library(ape)
library(MBECS)
library(stringi)
library(plyr)
library(tidyverse)
library(readxl)
library(tibble)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(reshape2)
library(devtools)
library(biomformat)
library(Biobase)
library(permute)
library(lattice)
library(plotly)
library(venn)
library(RCy3)
library(igraph)
library(philentropy)
library(vegan)
library(DESeq2)
library(phyloseq)
library(mixOmics)
library(ggpubr)
library(coin)
library(ggnewscale)
library(VennDiagram)
library(UpSetR)
library(ANCOMBC)
library(ARTool)
library(microbiome)
library(DT)
library(pacman)
library(ggplot2)
library(gplots)
library(dplyr)
library(vcd)
library(ggstatsplot)
library(gtools)
library(ggridges)
library(rstantools)
library(ggside)
library(ggcorrplot)
library(gginnards)
library(ggalluvial)
library(ggh4x)
library(magrittr)
library(corrr)
library(coin)
library(viridis)
library(RColorBrewer)
library(Hotelling)
library(knitr)
library(openxlsx)
library(rgl)
library(bladderbatch)
library(sva)
library(SRS)
library(doParallel)
library(patchwork)
library(apeglm)
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(here)

packageVersion("vegan")
packageVersion("DESeq2")
packageVersion("phyloseq")
packageVersion("mixOmics")
packageVersion("ANCOMBC")

renv::snapshot()

####### Load data and create phyloseq object for mid-colon #######
otu_mat<- read_xlsx(here::here("Raw data/03_MIC_Genus_Counts_corrected696.xlsx"))
tax_mat<- read_xlsx(here::here("Raw data/03_MIC_taxanomic table.xlsx"))
samples_df <- read_xlsx(here::here("Raw data/03_MIC_Sample data.xlsx"))
identical(rownames(tax_mat), rownames(otu_mat))

otu_matD <- data.frame(column_to_rownames(otu_mat, var = "otu"))
tax_matD <- data.frame(column_to_rownames(tax_mat, var = "otu"))
samples_dfD <- data.frame(column_to_rownames(samples_df, var = "sID"))
identical(rownames(tax_matD), rownames(otu_matD))

samples_dfD <- samples_dfD %>%
  mutate(cluster = as.factor(case_when(
    cluster == 1 ~ 'LacH',
    cluster == 2 ~ 'LacL',
    TRUE ~ as.character(cluster) 
  )))

otu_matM <- as.matrix(otu_matD)
tax_matM <- as.matrix(tax_matD)
samples_dfD$cluster <- factor(samples_dfD$cluster)
samples_dfD$Study <- factor(samples_dfD$Study)

OTU = phyloseq::otu_table(otu_matM, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax_matM)
samples = phyloseq::sample_data(samples_dfD)

OTU
TAX
samples

edi <- phyloseq(OTU, TAX, samples)

edi

####### ANCOMBC2 for mid-colon #######

sample_data(edi)$cluster <- factor(sample_data(edi)$cluster,
                                   levels = c("LacL", "LacH"))
sample_data(edi)$Study <- factor(sample_data(edi)$Study)

set.seed(1234)
out <- ancombc2(
  data = edi,
  tax_level = NULL,          
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 0,
  fix_formula = "Study + cluster",
  rand_formula = NULL,
  group = "cluster",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = FALSE,
  pairwise = FALSE,          
  dunnet = FALSE
)

res <- out$res

lfc_col    <- grep("^lfc_cluster", colnames(res), value = TRUE)
q_col      <- grep("^q_cluster", colnames(res), value = TRUE)
robust_col <- grep("^diff_robust_cluster", colnames(res), value = TRUE)

ancom_res <- res %>%
  transmute(
    Genus = taxon,
    lfc = .data[[lfc_col]],
    q_value = .data[[q_col]],
    diff_robust = .data[[robust_col]]
  ) %>%
  mutate(
    enriched_in = case_when(
      lfc > 0 ~ "LacH",
      lfc < 0 ~ "LacL",
      TRUE ~ NA_character_
    )
  )

#  Calculate mean relative abundance by cluster
ps_rel <- transform_sample_counts(edi, function(x) x / sum(x))

otu_rel <- as(otu_table(ps_rel), "matrix")
if (!taxa_are_rows(ps_rel)) otu_rel <- t(otu_rel)

mean_ra <- as.data.frame(t(otu_rel)) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(
    data.frame(sample_data(ps_rel)) %>%
      rownames_to_column("Sample_ID") %>%
      dplyr::select(Sample_ID, cluster),
    by = "Sample_ID"
  ) %>%
  group_by(cluster) %>%
  dplyr::summarise(across(-Sample_ID, mean), .groups = "drop") %>%
  pivot_longer(-cluster, names_to = "Genus", values_to = "Mean_RA") %>%
  pivot_wider(names_from = cluster, values_from = Mean_RA, names_prefix = "Mean_RA_")


#  Keep significant and robust genera

sig_genus <- ancom_res %>%
  left_join(mean_ra, by = "Genus") %>%
  filter(diff_robust == TRUE) %>%
  arrange(q_value, desc(abs(lfc)))

write.csv(
  sig_genus,
  here::here("03_ANCOMBC2_sig_genera_cluster_MIC.csv"),
  row.names = FALSE
)

plot_df <- sig_genus %>%
  slice_max(order_by = abs(lfc), n = 25) %>%
  mutate(
    Genus = forcats::fct_reorder(Genus, lfc),
    enriched_in = factor(enriched_in, levels = c("LacL", "LacH"))
  )

cluster_cols <- c("LacL" = "#b168a4", "LacH" = "#72af70")

p_top <- ggplot(plot_df, aes(x = lfc, y = Genus, fill = enriched_in)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  scale_fill_manual(values = cluster_cols) +
  labs(
    x = "ANCOM-BC2 log fold change (LacH vs LacL)",
    y = NULL,
    fill = "Enriched in"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black")
  )

p_top


ggsave(
  filename = here::here("03_ANCOMBC2_top_differential_genera_by_cluster_MIC.png"),
  plot = p_top,
  width = 7,
  height = 5.5,
  dpi = 600
)






####### Load data and create phyloseq object for cecum #######
otu_mat<- read_xlsx(here::here("Raw data/03_CAE_Genus_Counts_corrected664.xlsx"))
tax_mat<- read_xlsx(here::here("Raw data/03_CAE_taxonomic genera table.xlsx"))
samples_df <- read_xlsx(here::here("Raw data/03_CAE_Sample data.xlsx"))
identical(rownames(tax_mat), rownames(otu_mat))

otu_matD <- data.frame(column_to_rownames(otu_mat, var = "otu"))
tax_matD <- data.frame(column_to_rownames(tax_mat, var = "otu"))
samples_dfD <- data.frame(column_to_rownames(samples_df, var = "sID"))
identical(rownames(tax_matD), rownames(otu_matD))

samples_dfD <- samples_dfD %>%
  mutate(cluster = as.factor(case_when(
    cluster == 1 ~ 'LacH',
    cluster == 2 ~ 'LacL',
    TRUE ~ as.character(cluster) 
  )))

otu_matM <- as.matrix(otu_matD)
tax_matM <- as.matrix(tax_matD)
samples_dfD$cluster <- factor(samples_dfD$cluster)
samples_dfD$Study <- factor(samples_dfD$Study)

OTU = phyloseq::otu_table(otu_matM, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax_matM)
samples = phyloseq::sample_data(samples_dfD)

OTU
TAX
samples

edi <- phyloseq(OTU, TAX, samples)

edi

#Clear
#PRACTICE: clear
#Remove ASV's with total abundance lower than 0.01% (0.0001) of the total depth. Store the result in a "ediclear" object
#Note: in many papers 0.005% is taken, but then too many non relevant taxa are kept in tax table

head(taxa_sums(edi))
total.depth <-sum(taxa_sums(edi))
total.depth
threshold <- 0.0001*total.depth
taxa.filter <-taxa_sums(edi) > threshold
ediclear<- prune_taxa(taxa.filter, edi)
ediclear
edi

# Check the number of taxa before and after filtering
ntaxa(edi)
ntaxa(ediclear)
ntaxa(ediclear) / ntaxa(edi)

sum(taxa_sums(edi))
sum(taxa_sums(ediclear))
sum(taxa_sums(ediclear)) / sum(taxa_sums(edi))

removed_taxa <- taxa_names(edi)[!taxa.filter]
summary(taxa_sums(edi)[removed_taxa])

# Rarefy to the minimum sample depth after filtering for alpha diversity analyses
tab <- otu_table(ediclear)
class(tab) <- "matrix"
tab <- t(tab)# transpose observations to rows

rarecurve(tab, step = 50, xlab="Sample size", 
          ylab="ASV's", label=FALSE, abline(v=min(rowSums(tab))))
ediclear.rare <- rarefy_even_depth(ediclear, rngseed = 20170215)
sample_sums(ediclear.rare)
sample_sums(ediclear)
edi
ediclear
ediclear.rare

alpha_div <- estimate_richness(ediclear.rare, measures = c("Chao1", "Shannon", "Simpson"))
sample_data(ediclear.rare)$SampleID <- rownames(sample_data(ediclear.rare))

sample_data(ediclear.rare) <- cbind(sample_data(ediclear.rare), alpha_div)

alpha_div
write.csv(alpha_div, "03_CAE_Alpha diversity_Genus.csv")

alpha_melt <- melt(as.data.frame(sample_data(ediclear.rare)),
                   measure.vars = c("Chao1", "Shannon", "Simpson"))

p <- ggplot(alpha_melt, aes(x = cluster, y = value, fill = cluster)) + 
  geom_violin(trim = FALSE, alpha = 0.4, color = NA) +  
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7, linewidth = 0.7) +  
  geom_jitter(aes(color = cluster), width = 0.15, size = 2.2, alpha = 0.7) +  
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  labs(title = "Genus Level", x = "", y = "Alpha Index Value") +
  scale_fill_manual(values = c("LacH" = "#72af70dc", "LacL" = "#b168a4")) +
  scale_color_manual(values = c("LacH" = "#72af70dc", "LacL" = "#b168a4")) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

p


ggsave(
  filename = here("Figures", "05_MIC_Alpha_diversity_Genus.png"),
  plot = p,
  width = 10,
  height = 7,
  dpi = 600
)

dat20 <- read_xlsx(here("Raw data/03_CAE_Alpha diversity_Genus.xlsx"))
dat20 <- dat20 %>%
  mutate(cluster = as.factor(case_when(
    cluster == 1 ~ 'LACH',
    cluster == 2 ~ 'LACL',
    TRUE ~ as.character(cluster) 
  )))

dat20$cluster <- factor(dat20$cluster)
dat20$Study <- factor(dat20$Study)

dat21 <- dat20[!is.na(dat20$Chao1), ]
NonChao1 <- art(Chao1 ~ cluster + (1|Study), data = dat21)
anova(NonChao1)

dat22 <- dat20[!is.na(dat20$Shannon), ]
NonShannon <- art(Shannon ~ cluster + (1|Study), data = dat22)
anova(NonShannon)

dat23 <- dat20[!is.na(dat20$Simpson), ]
NonSimpson <- art(Simpson ~ cluster + (1|Study), data = dat23)
anova(NonSimpson)

# Beta diversity with edi.clear
ps_rel <- transform_sample_counts(ediclear, function(x) x / sum(x))
otu_matB <- otu_table(ps_rel)
sample_dataB <- sample_data(ps_rel)
bc_dist <- phyloseq::distance(ps_rel, method = "bray")
pcoa <- ordinate(ps_rel, method = "PCoA", distance = "bray")

pcoa_df <- as.data.frame(pcoa$vectors[,1:2])  # Assuming you only need the first two principal coordinates
colnames(pcoa_df) <- c("PCo1", "PCo2")
pcoa_df$cluster <- sample_data(ps_rel)$cluster
#posthoc_results
adonis_result <- adonis2(bc_dist ~ cluster, data = as(sample_dataB, "data.frame"), permutations = 9999)
adonis_result
sample_dataB$cluster <- factor(sample_dataB$cluster)
posthoc_results <-pairwise.perm.manova(bc_dist, sample_dataB$cluster,nperm = 9999, 
                                       test = c(wilks),p.method = "fdr", F = TRUE, R2 = TRUE)
print(posthoc_results)

#Create the PCoA plot with ellipses
var_explained <- pcoa$values$Relative_eig 

pcoa_plot <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = cluster)) +
  geom_point(size = 4) +
  stat_ellipse(
    type = "t",
    linetype = 1,
    linewidth = 1,
    alpha = 0.5,
    aes(fill = cluster),
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("LacH" = "#72af70dc", "LacL" = "#b168a4")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCoA of Bray-Curtis Distances",
    x = paste("PCo1 (", round(100 * var_explained[1], 2), "%)", sep = ""),
    y = paste("PCo2 (", round(100 * var_explained[2], 2), "%)", sep = "")
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(size = 13, face = "bold", color = "black"),
    legend.text = element_text(size = 12, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.1),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )

pcoa_plot
print(pcoa_plot)
# Calculate centroids
centroids <- aggregate(cbind(PCo1, PCo2) ~ cluster, data = pcoa_df, FUN = mean)
pcoa_plot <- pcoa_plot +
  geom_point(data = centroids, aes(x = PCo1, y = PCo2), 
             shape = 4, size = 5, stroke = 2, color = "black") +
  geom_text(data = centroids, aes(x = PCo1, y = PCo2, label = cluster), 
            vjust = -1, fontface = "bold", color = "black")
pcoa_plot
ggsave(
  filename = here::here("Figures", "05_MIC_PCoA_Bray_Curtis_Genus.png"),
  plot = pcoa_plot,
  width = 8,
  height = 6,
  dpi = 600
)


# Add to plot
pcoa_plot <- pcoa_plot +
  geom_point(data = centroids, aes(x = PCo1, y = PCo2), 
             shape = 4, size = 5, stroke = 2, color = "black") +
  geom_text(data = centroids, aes(x = PCo1, y = PCo2, label = cluster), 
            vjust = -1, fontface = "bold", color = "black")

sample_dataB_df <- as(sample_dataB, "data.frame")
bd <- betadisper(bc_dist, sample_dataB_df$cluster)
anova(bd)


####### ANCOMBC2 for cecum #######

sample_data(edi)$cluster <- factor(sample_data(edi)$cluster,
                                   levels = c("LacL", "LacH"))
sample_data(edi)$Study <- factor(sample_data(edi)$Study)

set.seed(1234)
out <- ancombc2(
  data = edi,
  tax_level = NULL,          
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 0,
  fix_formula = "Study + cluster",
  rand_formula = NULL,
  group = "cluster",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = FALSE,
  pairwise = FALSE,          
  dunnet = FALSE
)

res <- out$res

lfc_col    <- grep("^lfc_cluster", colnames(res), value = TRUE)
q_col      <- grep("^q_cluster", colnames(res), value = TRUE)
robust_col <- grep("^diff_robust_cluster", colnames(res), value = TRUE)

ancom_res <- res %>%
  transmute(
    Genus = taxon,
    lfc = .data[[lfc_col]],
    q_value = .data[[q_col]],
    diff_robust = .data[[robust_col]]
  ) %>%
  mutate(
    enriched_in = case_when(
      lfc > 0 ~ "LacH",
      lfc < 0 ~ "LacL",
      TRUE ~ NA_character_
    )
  )

#  Calculate mean relative abundance by cluster
ps_rel <- transform_sample_counts(edi, function(x) x / sum(x))

otu_rel <- as(otu_table(ps_rel), "matrix")
if (!taxa_are_rows(ps_rel)) otu_rel <- t(otu_rel)

mean_ra <- as.data.frame(t(otu_rel)) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(
    data.frame(sample_data(ps_rel)) %>%
      rownames_to_column("Sample_ID") %>%
      dplyr::select(Sample_ID, cluster),
    by = "Sample_ID"
  ) %>%
  group_by(cluster) %>%
  dplyr::summarise(across(-Sample_ID, mean), .groups = "drop") %>%
  pivot_longer(-cluster, names_to = "Genus", values_to = "Mean_RA") %>%
  pivot_wider(names_from = cluster, values_from = Mean_RA, names_prefix = "Mean_RA_")


#  Keep significant and robust genera

sig_genus <- ancom_res %>%
  left_join(mean_ra, by = "Genus") %>%
  filter(diff_robust == TRUE) %>%
  arrange(q_value, desc(abs(lfc)))



# Ancombc for caecum 
fill_cols <- c("LacH" = "#72af70dc", "LacL" = "#b168a4")
line_cols <- c("LacH" = "#4E8C4A", "LacL" = "#8D4F83")

top2 <- sig_genus$Genus

dat_plot <- phyloseq::psmelt(ps_rel) %>%
  dplyr::filter(Genus %in% top2) %>%
  dplyr::mutate(
    Genus = factor(Genus, levels = c("Oscillospiraceae UCG-005", "Dialister")),
    cluster = factor(cluster, levels = c("LacH", "LacL"))
  )

p_box <- ggplot(dat_plot, aes(x = cluster, y = Abundance, fill = cluster)) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85,
    color = "black",
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(color = cluster),
    width = 0.12,
    size = 2,
    alpha = 0.75
  ) +
  facet_wrap(~Genus, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = fill_cols) +
  scale_color_manual(values = line_cols) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(x = NULL, y = "Relative abundance, %") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(color = "black")
  )

p_box

# Visualization
library(patchwork)
library(ggplot2)

p1 <- p +
  labs(title = NULL) +
  theme(
    legend.position = "none",
    plot.title = element_blank()
  )

p2 <- pcoa_plot +
  labs(title = NULL) +
  theme(
    legend.position = "right",
    plot.title = element_blank()
  )

p3 <- p_box +
  labs(title = NULL) +
  theme(
    legend.position = "none",
    plot.title = element_blank()
  )

p_combined <- p1 / (p3 | p2) +
  plot_layout(
    heights = c(1.05, 1.00),
    widths = c(1.05, 0.95)
  ) +
  plot_annotation(
  tag_levels = "a",
  tag_prefix = "(",
  tag_suffix = ")"
)

p_combined


ggsave(
  filename = here::here("05_CAE_Microbial_analysis_summary_CAE1.png"),
  plot = p_combined,
  width = 20,
  height = 15,
  dpi = 600,
  bg = "white"
)
