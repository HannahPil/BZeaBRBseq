library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(grid)  # for unit()

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
# filter + normalize (original keep rule)
# keep genes with >=10 counts in >=2 samples
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

# wide -> long: gene_id x genotype that carries teosinte
carriers <- allelic %>%
  pivot_longer(
    cols = everything(),
    names_to  = "gene_id",
    values_to = "genotype"
  ) %>%
  filter(!is.na(genotype), genotype != "")

# ============================
# candidate genes + categories (annotation with clear priority)
# ============================
cand <- read.csv("candidate_genes.csv")

keep_cats <- c("FT", "targ", "GWAS_GBS_landraces_N", "Fst_landraces_N")

cand_sub_raw <- cand %>%
  filter(
    gene_id %in% genes_from_allelic,
    gene_id %in% rownames(mat)
  )

cand_sub <- cand_sub_raw %>%
  group_by(gene_id) %>%
  summarise(
    category = {
      cats <- unique(category)
      # keep only categories in your main nitrogen-related set
      cats_keep <- keep_cats[keep_cats %in% cats]
      
      if (length(cats_keep) == 0) {
        NA_character_
      } else {
        # choose the first match based on the order in keep_cats
        cats_keep[1]
      }
    },
    .groups = "drop"
  )

# ============================
# load gene names (for aesthetics)
# expects columns: v5_gene_id, gene_name
# ============================
gene_names <- read.csv("gene_names.csv") %>%
  rename(gene_id = v5_gene_id)

# ============================
# subset normalized matrix to allelic genes that passed filter
# ============================
genes_to_use <- intersect(genes_from_allelic, rownames(mat))
mat_sub <- mat[rownames(mat) %in% genes_to_use, ]

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
# category factor for plotting
# ============================
df_gene_taxa <- df_gene_taxa %>%
  mutate(
    category_plot = ifelse(!is.na(category) & category %in% keep_cats,
                           category,
                           "other/NA"),
    category_plot = factor(category_plot, levels = c(keep_cats, "other/NA"))
  )

# ============================
# add gene names (aesthetic only)
# if no gene_name, fall back to gene_id
# ============================
df_gene_taxa <- df_gene_taxa %>%
  left_join(gene_names, by = "gene_id") %>%
  mutate(
    gene_label = ifelse(!is.na(gene_name) & gene_name != "",
                        gene_name,
                        gene_id)
  )

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
# NEW: drop BZea lines that never carry teosinte for any gene
# (but always keep B73 and Purple Check)
# ============================
genos_with_teo <- unique(carriers$genotype)

df_gene_taxa <- df_gene_taxa %>%
  filter(genotype %in% c(genos_with_teo, "B73", "Purple Check"))

# ============================
# ordering for axes
# ============================
# order genes within category by gene_id, but display gene_label
gene_order <- df_gene_taxa %>%
  arrange(category_plot, gene_id) %>%
  pull(gene_label) %>%
  unique()

geno_order <- df_gene_taxa %>%
  arrange(Taxa, genotype) %>%
  pull(genotype) %>%
  unique()

df_gene_taxa <- df_gene_taxa %>%
  mutate(
    gene_label = factor(gene_label, levels = rev(gene_order)),
    genotype   = factor(genotype, levels = geno_order),
    Taxa       = factor(Taxa)
  )

# ============================
# big heatmap!
# – colored only where teosinte exists
# – gray where maize-only
# – B73 + Purple Check always shown
# – only genotypes that have teosinte for at least one gene (plus checks)
# – rows faceted by category_plot
# – y-axis shows gene names where available
# ============================
ggplot(df_gene_taxa, aes(x = genotype, y = gene_label, fill = centered_logCPM_masked)) +
  geom_tile() +
  facet_grid(
    category_plot ~ Taxa,
    scales = "free",
    space  = "free"
  ) +
  scale_fill_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    na.value = "grey80"
  ) +
  labs(
    x = "BZea genotype",
    y = "gene",
    fill = "centered logCPM",
    title = "expression of candidate genes across bzea lines"
  ) +
  theme_bw() +
  theme(
    panel.spacing    = unit(0, "pt"),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x      = element_text(angle = 90,
                                    vjust = 0.5,
                                    hjust = 1,
                                    size  = 6),
    axis.text.y      = element_text(size = 6)
  )
