library(edgeR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(limma)
library(purrr)

# files
counts_file <- "Zea_mays_counts.txt"
meta_file   <- "metadata.csv"

# reference taxa for all contrasts
ref_taxa <- "B73"

# make output folder for plots
out_dir <- file.path("output", "genomeWide")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- read counts correctly: gene ids are rownames ----
counts <- read.table(
  counts_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

count_mat <- as.matrix(counts)

# ---- read metadata (csv) ----
meta <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)

# keep only plates 1-4 (the sequenced plates)
meta <- subset(meta, plate %in% c(1, 2, 3, 4))

# restrict metadata to samples present in counts, then reorder to match counts
meta <- meta[meta$sample_id %in% colnames(count_mat), ]
count_mat <- count_mat[, meta$sample_id]

# coerce columns
meta <- meta %>%
  mutate(
    taxa  = factor(taxa),
    Row   = as.numeric(Row),
    Range = as.numeric(Range)
  )

# set B73 as the reference level for taxa (so logFC = taxa vs B73)
meta$taxa <- relevel(meta$taxa, ref = ref_taxa)

# drop any unused taxa levels after subsetting
meta <- meta %>% mutate(taxa = droplevels(taxa))

# ---- edgeR on raw counts ----
y <- DGEList(counts = count_mat)
y <- calcNormFactors(y)

# filter low expression genes (standard)
keep <- filterByExpr(y, group = meta$taxa)
y <- y[keep, , keep.lib.sizes = FALSE]

# design: taxa + Row + Range
design <- model.matrix(
  object = reformulate(termlabels = c("taxa", "Row", "Range")),
  data   = meta
)

# fit and test (omnibus taxa effect controlling for space)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

taxa_coefs <- grep("^taxa", colnames(design))
qlf_taxa <- glmQLFTest(fit, coef = taxa_coefs)

# results: one row per gene
res_taxa <- topTags(qlf_taxa, n = Inf)$table %>%
  tibble::rownames_to_column("gene_id")

write.csv(res_taxa, file.path(out_dir, "edgeR_results_taxa_plus_space.csv"), row.names = FALSE)

# log2 CPM for plotting
logcpm <- cpm(y, log = TRUE, prior.count = 1)
write.csv(
  data.frame(gene_id = rownames(logcpm), logcpm, check.names = FALSE),
  file.path(out_dir, "edgeR_log2cpm_TMM_filtered.csv"),
  row.names = FALSE
)

#---------------------------------------------------------------------------------------
#-------------------- volcanoes + summaries --------------------------------------------
#---------------------------------------------------------------------------------------

keep_cats <- c("FT", "targ", "GWAS_GBS_landraces_N", "Fst_landraces_N")

# =========================================
# 1) faceted volcano plots: each taxa vs B73
# =========================================
taxa_levels <- levels(meta$taxa)
other_taxa <- setdiff(taxa_levels, ref_taxa)

# build contrast names (treatment coding yields coefficients named like "taxa<LEVEL>")
contrast_names <- paste0("taxa", other_taxa)
contrast_names <- contrast_names[contrast_names %in% colnames(design)]

C <- makeContrasts(contrasts = contrast_names, levels = design)

volcano_df <- map_dfr(colnames(C), function(cname) {
  tt <- glmQLFTest(fit, contrast = C[, cname]) %>%
    topTags(n = Inf) %>%
    (\(x) x$table)() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    mutate(contrast = cname)
  tt
})

volcano_df <- volcano_df %>%
  mutate(
    taxa_comp = sub("^taxa", "", contrast),
    label = paste0(taxa_comp, " vs ", ref_taxa),
    neglog10FDR = -log10(FDR)
  )

# -----------------------------------------
# candidate categories: join + focus
# -----------------------------------------
cand <- read.csv("candidate_genes.csv", stringsAsFactors = FALSE)

cand_cat <- cand %>%
  group_by(gene_id) %>%
  summarise(category = first(category), .groups = "drop")

volcano_df <- volcano_df %>%
  left_join(cand_cat, by = "gene_id") %>%
  mutate(
    category_focus = ifelse(!is.na(category) & category %in% keep_cats, category, "other"),
    category_focus = factor(category_focus, levels = c("other", keep_cats))
  )

# -----------------------------------------
# plot: color only keep_cats, soften everything else (including other candidate categories)
# -----------------------------------------
p_volcano <- ggplot(volcano_df, aes(x = logFC, y = neglog10FDR)) +
  
  # background: everything else (including other candidate categories)
  geom_point(
    data = subset(volcano_df, category_focus == "other"),
    color = "grey70",
    alpha = 0.7,
    size = 0.9
  ) +
  
  # foreground: focus categories, drawn last so they sit on top
  geom_point(
    data = subset(volcano_df, category_focus != "other"),
    aes(color = category_focus),
    alpha = 0.9,
    size = 1.2
  ) +
  
  facet_wrap(~ label, scales = "free_y") +
  
  scale_color_manual(
    values = c(
      "FT" = "#e41a1c",
      "targ" = "#377eb8",
      "GWAS_GBS_landraces_N" = "#4daf4a",
      "Fst_landraces_N" = "#984ea3"
    ),
    drop = FALSE
  ) +
  
  labs(
    x = "log2 fold change (taxa vs B73)",
    y = "-log10(FDR)",
    color = "candidate category",
    title = "Volcano plots by taxa contrast (focus categories emphasized; controlling for Row + Range)"
  ) +
  theme_bw()

p_volcano

ggsave(file.path(out_dir, "volcano_taxa_vs_B73_by_candidate_category_colors.png"), p_volcano, width = 12, height = 8, dpi = 300)


#---------------gene-specific-modeling---------------------
#---------------gene-specific-modeling---------------------
#---------------gene-specific-modeling---------------------

library(tidyverse)

# ----------------------------
# inputs
# ----------------------------
gene_id <- "Zm00001eb121780"

expr_file     <- file.path(out_dir, "edgeR_log2cpm_TMM_filtered.csv")
meta_file     <- "metadata.csv"
allelic_file  <- "Allelic_series_for_expression.csv"

# ----------------------------
# load expression
# ----------------------------
expr <- read.csv(expr_file, check.names = FALSE)

expr_long <- expr %>%
  filter(gene_id == !!gene_id) %>%
  pivot_longer(
    cols = -gene_id,
    names_to  = "sample_id",
    values_to = "logCPM"
  )

# ----------------------------
# load metadata (plates 1-4 only)
# ----------------------------
meta <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  filter(plate %in% c(1, 2, 3, 4)) %>%
  mutate(
    taxa  = factor(taxa),
    Row   = as.numeric(Row),
    Range = as.numeric(Range)
  )

df <- expr_long %>%
  left_join(meta, by = "sample_id")

# ----------------------------
# load allelic series
# ----------------------------
allelic <- read.csv(allelic_file, check.names = FALSE)

# genotypes that carry teosinte for THIS gene
teo_carriers <- allelic[[gene_id]] %>%
  na.omit() %>%
  unique()

df <- df %>%
  mutate(
    has_teo = genotype %in% teo_carriers,
    
    # always keep checks as reference (same logic you used before)
    has_teo = ifelse(genotype %in% c("B73", "Purple Check"),
                     FALSE,
                     has_teo),
    
    has_teo = factor(has_teo, levels = c(FALSE, TRUE))
  )

# ----------------------------
# fit gene-level model
# ----------------------------
fit <- lm(
  logCPM ~ taxa + has_teo + Row + Range,
  data = df
)

summary(fit)

#-----------------------------------------gtf stuff------------------
library(vroom)
library(dplyr)
library(stringr)

gtf_path <- "Zea_mays.gtf"  # can be .gtf or .gtf.gz

gtf_tx <- vroom::vroom(
  file = gtf_path,
  delim = "\t",
  col_names = FALSE,
  comment = "#",
  col_select = c(X1, X3, X4, X5, X9),
  progress = TRUE
) %>%
  setNames(c("seqname", "type", "start", "end", "attr")) %>%
  filter(type == "transcript") %>%
  mutate(
    gene_id = str_match(attr, 'gene_id "([^"]+)"')[,2],
    chr = as.character(seqname)
  ) %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id, chr) %>%
  summarise(
    start = min(as.integer(start)),
    end   = max(as.integer(end)),
    .groups = "drop"
  )

cat("genes parsed:", nrow(gtf_tx), "\n")
print(head(gtf_tx, 10))

# optional: quick check if names match teogeno-style chromosomes
print(head(sort(unique(gtf_tx$chr)), 20))
# print(head(sort(unique(seg_df$chr)), 20))  # requires seg_df from 13a/13b



