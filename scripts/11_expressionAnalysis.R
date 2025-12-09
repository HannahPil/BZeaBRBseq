library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

setwd("C:/Users/Hannah Pil/Documents/gemmalab/BZea/BZea RNA-seq/BZeaBRBseq")

# ============================
# load data
# ============================
counts <- read.table("Zea_mays_counts.txt",
                     header = TRUE,
                     row.names = 1,
                     sep = "\t")

meta_all <- read.csv("metadata.csv")

# keep only complete rows (plate 1–4)
meta <- subset(meta_all, plate %in% c(1, 2, 3, 4))

# restrict metadata to samples present in counts
meta <- meta[meta$sample_id %in% colnames(counts), ]

# reorder counts to match metadata order
counts <- counts[, meta$sample_id]

# ============================
# filter + normalize
# ============================
keep <- rowSums(counts >= 10) >= 2
counts_filt <- counts[keep, ]

dge <- DGEList(counts = counts_filt)
dge <- calcNormFactors(dge)

# log2 cpm with prior count
mat <- cpm(dge, log = TRUE, prior.count = 1)

# ============================
# load allelic series (genes to include)
# ============================
allelic <- read.csv("Allelic_series_for_expression.csv",
                    check.names = FALSE)

genes_from_allelic <- colnames(allelic)

# turn wide → long for teosinte carriers
carriers <- allelic %>%
  pivot_longer(
    cols = everything(),
    names_to  = "gene_id",
    values_to = "genotype"
  ) %>%
  filter(!is.na(genotype), genotype != "")

# ============================
# candidate genes + categories
# ============================
cand <- read.csv("candidate_genes.csv")

keep_cats <- c("FT", "targ", "GWAS_GBS_landraces_N", "Fst_landraces_N")

cand_sub_raw <- cand %>%
  filter(
    gene_id %in% genes_from_allelic,   # only analyze genes in allelic-series CSV
    category %in% keep_cats,
    gene_id %in% rownames(mat)
  )

# collapse multiple categories per gene (rare)
cand_sub <- cand_sub_raw %>%
  group_by(gene_id) %>%
  summarise(
    category = paste(unique(category), collapse = ";"),
    .groups = "drop"
  )

# ============================
# subset normalized matrix to selected genes
# ============================
mat_sub <- mat[rownames(mat) %in% cand_sub$gene_id, ]

# ============================
# long format expression table
# ============================
df_long <- as.data.frame(mat_sub) %>%
  mutate(gene_id = rownames(mat_sub)) %>%
  pivot_longer(
    cols = -gene_id,
    names_to  = "sample_id",
    values_to = "logCPM"
  ) %>%
  left_join(cand_sub, by = "gene_id") %>%
  left_join(meta %>% select(sample_id, genotype, taxa),
            by = "sample_id") %>%
  mutate(Taxa = taxa)

# ============================
# average per gene × genotype × taxa
# ============================
df_gene_taxa <- df_long %>%
  group_by(gene_id, category, Taxa, genotype) %>%
  summarise(mean_logCPM = mean(logCPM, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(gene_id) %>%
  mutate(centered_logCPM = mean_logCPM - mean(mean_logCPM)) %>%
  ungroup()

# ============================
# add teosinte-carrier info
# + force B73 and Purple Check to always show expression
# ============================
df_gene_taxa <- df_gene_taxa %>%
  left_join(
    carriers %>% mutate(has_teo = TRUE),
    by = c("gene_id", "genotype")
  ) %>%
  mutate(
    has_teo = ifelse(is.na(has_teo), FALSE, TRUE),
    
    # B73 + Purple Check always displayed
    has_teo = ifelse(genotype %in% c("B73", "Purple Check"), TRUE, has_teo),
    
    centered_logCPM_masked = ifelse(has_teo, centered_logCPM, NA)
  )

# ============================
# ordering for axes
# ============================
gene_order <- df_gene_taxa %>%
  arrange(factor(category, levels = keep_cats), gene_id) %>%
  pull(gene_id) %>%
  unique()

geno_order <- df_gene_taxa %>%
  arrange(Taxa, genotype) %>%
  pull(genotype) %>%
  unique()

df_gene_taxa <- df_gene_taxa %>%
  mutate(
    gene_id  = factor(gene_id, levels = rev(gene_order)),
    genotype = factor(genotype, levels = geno_order),
    category = factor(category, levels = keep_cats),
    Taxa     = factor(Taxa)
  )

# ============================
# big heatmap!
# – colored only where teosinte exists
# – gray where maize-only
# – B73 + Purple Check always shown
# ============================
ggplot(df_gene_taxa, aes(x = genotype, y = gene_id, fill = centered_logCPM_masked)) +
  geom_tile() +
  # Use both scales and space set to "free_x" to vary panel widths
  facet_grid(~ Taxa, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low     = "blue",
    mid     = "white",
    high    = "red",
    midpoint = 0,
    na.value = "grey80"
  ) +
  labs(
    x = "bzea line (genotype)",
    y = "gene (grouped by category)",
    fill = "centered logCPM",
    title = "expression of candidate genes across bzea lines and taxa\n(red/blue only where teosinte introgression is present;\nB73 & Purple Check always shown)"
  ) +
  theme_bw() +
  theme(
    # This removes the white spacing (gutters) between the facets
    panel.spacing = unit(0, "pt"),
    axis.text.x = element_text(angle = 90,
                               vjust = 0.5,
                               hjust = 1,
                               size = 6),
    axis.text.y = element_text(size = 6)
  )
