# ==============================================================================
# BZEA RNA: GENOME-WIDE EQTL ARCHITECTURE (TOP 10% CIS + TOP 10% TRANS)
# ==============================================================================
# inputs:
#   cis  : output/cis_eQTL/BZea_Significant_eQTLs.csv
#   trans: output/trans_by_source_chr/BZea_trans_ALL_sources_combined_pv1e-04.tsv.gz
#   gtf  : Zea_mays.gtf
#
# plot:
#   x = SOURCE gene genomic position
#   y = TARGET gene genomic position
#   keep ONLY the top 10% (smallest p-values) within cis and within trans
#
# outputs:
#   output/eqtl_architecture_from_csv/
#     plot_eqtl_contactmap_TOP10pct_bin2d.png
#     plot_eqtl_contactmap_TOP10pct_joined_data.csv
#
# notes:
# - trans .tsv.gz is read directly by read.delim()
# - includes dplyr::select and dplyr::rename to avoid masking issues
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(scales)
})

# ==============================================================================
# 0. FILES + SETTINGS
# ==============================================================================
cis_csv  <- file.path("output", "cis_eQTL", "BZea_Significant_eQTLs.csv")
trans_gz <- file.path("output", "trans_by_source_chr", "BZea_trans_ALL_sources_combined_pv1e-04.tsv.gz")
data_dir <- "data"
gtf_file <- file.path(data_dir, "Zea_mays.gtf")

out_dir <- file.path("output", "eqtl_architecture_from_csv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

top_frac <- 0.10
bins_n   <- 350

out_png <- file.path(out_dir, paste0("plot_eqtl_contactmap_TOP", round(top_frac*100), "pct_bin2d.png"))
out_csv <- file.path(out_dir, paste0("plot_eqtl_contactmap_TOP", round(top_frac*100), "pct_joined_data.csv"))

stopifnot(file.exists(cis_csv), file.exists(trans_gz), file.exists(gtf_file))

cat("settings:\n")
cat("  top_frac = ", top_frac, "\n", sep = "")
cat("  bins_n   = ", bins_n, "\n", sep = "")
cat("  cis_csv  = ", cis_csv, "\n", sep = "")
cat("  trans_gz = ", trans_gz, "\n", sep = "")
cat("  gtf_file = ", gtf_file, "\n", sep = "")
cat("  out_dir  = ", out_dir, "\n", sep = "")

# ==============================================================================
# 1. LOAD RESULTS
# ==============================================================================
cat("\n1. Loading cis + trans...\n")

cis_raw <- read.csv(cis_csv, stringsAsFactors = FALSE)

# trans: tab-delimited gz
trans_raw <- read.delim(
  trans_gz,
  sep = "\t",
  stringsAsFactors = FALSE
)

cat("  cis rows:   ", nrow(cis_raw), "\n", sep = "")
cat("  trans rows: ", nrow(trans_raw), "\n", sep = "")

# standardize column names
normalize_cols <- function(df) {
  
  nms <- names(df)
  
  if(!("SNP" %in% nms) && ("snps" %in% nms)) df <- df |> dplyr::rename(SNP = snps)
  if(!("SNP" %in% nms) && ("snp"  %in% nms)) df <- df |> dplyr::rename(SNP = snp)
  
  if(!("gene" %in% nms) && ("geneid" %in% nms)) df <- df |> dplyr::rename(gene = geneid)
  
  nms <- names(df)
  if(!("p.value" %in% nms) && ("pvalue" %in% nms)) df <- df |> dplyr::rename(p.value = pvalue)
  if(!("p.value" %in% nms) && ("pval"   %in% nms)) df <- df |> dplyr::rename(p.value = pval)
  if(!("p.value" %in% nms) && ("p_value"%in% nms)) df <- df |> dplyr::rename(p.value = p_value)
  
  df
}

cis_df <- normalize_cols(cis_raw) |>
  mutate(
    type = "Cis",
    p.value = as.numeric(p.value)
  ) |>
  filter(!is.na(p.value))

trans_df <- normalize_cols(trans_raw) |>
  mutate(
    type = "Trans",
    p.value = as.numeric(p.value)
  ) |>
  filter(!is.na(p.value))

stopifnot("SNP" %in% names(cis_df), "gene" %in% names(cis_df))
stopifnot("SNP" %in% names(trans_df), "gene" %in% names(trans_df))

# keep TOP fraction within each set (smallest p-values)
keep_top_frac <- function(df, frac) {
  n_keep <- ceiling(nrow(df) * frac)
  df |>
    arrange(p.value) |>
    slice_head(n = n_keep)
}

cis_top   <- keep_top_frac(cis_df, top_frac)
trans_top <- keep_top_frac(trans_df, top_frac)

cat("  cis kept (top):   ", nrow(cis_top), "\n", sep = "")
cat("  trans kept (top): ", nrow(trans_top), "\n", sep = "")

eqtl_top <- bind_rows(cis_top, trans_top)
cat("  total kept: ", nrow(eqtl_top), "\n", sep = "")

# ==============================================================================
# 2. BUILD GENE COORDINATES
# ==============================================================================
cat("\n2. Building gene coordinates...\n")

gtf <- import(gtf_file)
gtf_df <- as.data.frame(gtf)

gene_coords <- gtf_df |>
  filter(!is.na(gene_id)) |>
  group_by(gene_id) |>
  summarise(
    chr   = as.character(seqnames[1]),
    start = min(start),
    end   = max(end),
    .groups = "drop"
  ) |>
  filter(grepl("^chr[0-9]+$", chr)) |>
  mutate(chr_num = as.numeric(gsub("chr", "", chr, ignore.case = TRUE))) |>
  arrange(chr_num, start)

chrom_info <- gene_coords |>
  group_by(chr, chr_num) |>
  summarise(length_bp = max(end), .groups = "drop") |>
  arrange(chr_num) |>
  mutate(
    tot_offset = lag(cumsum(length_bp), default = 0),
    center = tot_offset + length_bp / 2
  )

get_global_pos <- function(chr, pos, ref_df) {
  ref_df$tot_offset[match(chr, ref_df$chr)] + pos
}

# ==============================================================================
# 3. JOIN SOURCE + TARGET POSITIONS
# ==============================================================================
cat("\n3. Joining positions...\n")

src_pos <- gene_coords |>
  dplyr::select(gene_id, chr, start) |>
  dplyr::rename(SNP = gene_id, SNP_chr = chr, SNP_pos = start)

tgt_pos <- gene_coords |>
  dplyr::select(gene_id, chr, start) |>
  dplyr::rename(gene = gene_id, gene_chr = chr, gene_pos = start)

plot_df <- eqtl_top |>
  left_join(src_pos, by = "SNP") |>
  left_join(tgt_pos, by = "gene") |>
  filter(!is.na(SNP_chr), !is.na(gene_chr), !is.na(SNP_pos), !is.na(gene_pos)) |>
  mutate(
    x_global = get_global_pos(SNP_chr, SNP_pos, chrom_info) / 1e6,
    y_global = get_global_pos(gene_chr, gene_pos, chrom_info) / 1e6
  )

cat("  rows after coord join: ", nrow(plot_df), "\n", sep = "")
stopifnot(nrow(plot_df) > 0)

write.csv(plot_df, out_csv, row.names = FALSE)
cat("  saved joined table: ", out_csv, "\n", sep = "")

# ==============================================================================
# 4. PLOT CONTACT MAP (BIN2D)
# ==============================================================================
cat("\n4. Plotting...\n")

x_breaks <- chrom_info$center / 1e6
x_labels <- gsub("chr", "", chrom_info$chr)
y_breaks <- chrom_info$center / 1e6
y_labels <- gsub("chr", "", chrom_info$chr)

p <- ggplot(plot_df, aes(x = x_global, y = y_global)) +
  geom_vline(
    data = chrom_info,
    aes(xintercept = tot_offset / 1e6),
    color = "grey93",        # lighter
    linetype = "solid",      # solid instead of dashed
    linewidth = 0.5          # thinner
  ) +
  geom_hline(
    data = chrom_info,
    aes(yintercept = tot_offset / 1e6),
    color = "grey93",
    linetype = "solid",
    linewidth = 0.3
  ) +
  geom_bin2d(bins = bins_n) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0, 0)) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels, expand = c(0, 0)) +
  scale_fill_viridis_c(option = "magma", trans = "sqrt") +
  labs(
    title = paste0("Genome-wide eQTL Architecture (top ", round(top_frac * 100), "% cis + top ", round(top_frac * 100), "% trans)"),
    subtitle = "x = source gene; y = target gene",
    x = "Source genomic position (chromosome)",
    y = "Target genomic position (chromosome)",
    fill = "Density"
  ) +
  
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(face = "bold")
  )

p
ggsave(out_png, p, width = 9, height = 8)

cat("\nSaved:\n")
cat("  ", out_png, "\n", sep = "")
cat("  ", out_csv, "\n", sep = "")


# ==============================================================================
# 5. TRANS HOTSPOT CANDIDATES (TOP 10% TRANS SUBSET)
# ==============================================================================
# goal:
#   find source genes (SNP) that hit many target genes in trans_top
# outputs:
#   output/eqtl_architecture_from_csv/
#     trans_hotspot_sources_TOP10pct.csv
#     trans_hotspot_top_source_targets.csv
# ==============================================================================

cat("\n5. Ranking trans source genes by number of targets (top 10% trans)...\n")

# if you already have trans_top from the plotting script, use it.
# otherwise, rebuild trans_top here from the .tsv.gz using the same logic:
# trans_raw <- read.delim(trans_gz, sep = "\t", stringsAsFactors = FALSE)
# trans_df  <- normalize_cols(trans_raw) |> mutate(type = "Trans", p.value = as.numeric(p.value)) |> filter(!is.na(p.value))
# trans_top <- trans_df |> arrange(p.value) |> slice_head(n = ceiling(nrow(trans_df) * top_frac))

stopifnot(exists("trans_top"))

# count unique targets per source, plus optional sig counts
hotspot_rank <- trans_top |>
  dplyr::select(SNP, gene, p.value, FDR) |>
  group_by(SNP) |>
  summarise(
    n_targets = n_distinct(gene),
    n_tests   = dplyr::n(),
    n_sig_targets = n_distinct(gene[FDR < 0.05]),
    min_p     = min(p.value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_targets), desc(n_sig_targets), min_p)

# annotate hotspot sources with coordinates so you can see if they sit at chr6/7 edge
src_pos <- gene_coords |>
  dplyr::select(gene_id, chr, start, end) |>
  dplyr::rename(SNP = gene_id, SNP_chr = chr, SNP_start = start, SNP_end = end)

hotspot_rank_annot <- hotspot_rank |>
  left_join(src_pos, by = c("SNP")) |>
  mutate(
    SNP_chr_num = suppressWarnings(as.numeric(gsub("chr", "", SNP_chr, ignore.case = TRUE))),
    SNP_mid = (SNP_start + SNP_end) / 2
  ) |>
  arrange(desc(n_targets), desc(n_sig_targets), min_p)

out_hotspot <- file.path(out_dir, paste0("trans_hotspot_sources_TOP", round(top_frac*100), "pct.csv"))
write.csv(hotspot_rank_annot, out_hotspot, row.names = FALSE)

cat("  saved hotspot rank table:\n  ", out_hotspot, "\n", sep = "")
cat("\n  top 15 sources by n_targets:\n")
print(hotspot_rank_annot |> dplyr::select(SNP, SNP_chr, SNP_start, n_targets, n_sig_targets, min_p) |> head(15))

# pull targets for the top source gene so you can inspect what it's hitting
top_source <- hotspot_rank_annot$SNP[1]
cat("\n  top source SNP = ", top_source, "\n", sep = "")

top_source_targets <- trans_top |>
  filter(SNP == top_source) |>
  arrange(p.value) |>
  left_join(
    tgt_pos,  # from earlier script: gene -> gene_chr/gene_pos
    by = "gene"
  ) |>
  dplyr::select(SNP, gene, gene_chr, gene_pos, p.value, FDR, beta) |>
  distinct()

out_targets <- file.path(out_dir, "trans_hotspot_top_source_targets.csv")
write.csv(top_source_targets, out_targets, row.names = FALSE)

cat("  saved targets for top source:\n  ", out_targets, "\n", sep = "")
cat("  n targets (rows): ", nrow(top_source_targets), "\n", sep = "")

# ==============================================================================
# LOCUSZOOM-STYLE PLOT FOR SOURCE CHR6 (ALL TRANS, NO LD COLORING)
# - uses ALL trans results (not top 10%)
# - one point per source gene (SNP) in window
# - y = -log10(min p-value across its trans targets)
# ==============================================================================

# (tidyverse already loaded above)

window_mb <- 0.5
out_lz <- file.path(out_dir, "locuszoom_like_source_chr6_allTrans_noLD.png")

# ------------------------------ load ALL trans --------------------------------
# if trans_df already exists from earlier, this will reuse it; otherwise rebuild
if (!exists("trans_df")) {
  
  trans_raw <- read.delim(
    trans_gz,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  trans_df <- normalize_cols(trans_raw) |>
    dplyr::mutate(
      type = "Trans",
      p.value = as.numeric(p.value)
    ) |>
    dplyr::filter(!is.na(p.value))
  
  stopifnot("SNP" %in% names(trans_df), "gene" %in% names(trans_df))
}

# src_pos should already exist from your main script; if not, rebuild it here
if (!exists("src_pos")) {
  src_pos <- gene_coords |>
    dplyr::select(gene_id, chr, start) |>
    dplyr::rename(SNP = gene_id, SNP_chr = chr, SNP_pos = start)
}

# ------------------------------ choose lead source (forced) -------------------

lead_source <- "Zm00001eb297900"

lead_pos_bp <- src_pos |>
  dplyr::filter(SNP == lead_source) |>
  dplyr::slice(1) |>
  dplyr::pull(SNP_pos)

stopifnot(length(lead_pos_bp) == 1)

center_pos_mb <- lead_pos_bp / 1e6

cat("\nLead source on chr6 (forced):\n")
cat("  ", lead_source, " at ", round(center_pos_mb, 4), " Mb\n", sep = "")

# ------------------------------ build regional df ------------------------------
lz_df <- trans_df |>
  dplyr::left_join(src_pos, by = "SNP") |>
  dplyr::filter(SNP_chr == "chr6") |>
  dplyr::mutate(pos_mb = SNP_pos / 1e6) |>
  dplyr::filter(
    pos_mb >= (center_pos_mb - window_mb),
    pos_mb <= (center_pos_mb + window_mb)
  ) |>
  dplyr::group_by(SNP) |>
  dplyr::summarise(
    pos_mb = dplyr::first(pos_mb),
    p_min  = min(p.value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(log_p = -log10(p_min))

stopifnot(nrow(lz_df) > 0)

# ------------------------------ plot ------------------------------------------
p_lz <- ggplot(lz_df, aes(x = pos_mb, y = log_p)) +
  
  geom_point(
    color = "black",
    size = 3,
    alpha = 0.8
  ) +
  
  geom_point(
    data = dplyr::filter(lz_df, SNP == lead_source),
    shape = 23,
    size = 5,
    fill = "#8E44AD",
    color = "black"
  ) +
  
  geom_text(
    data = dplyr::filter(lz_df, SNP == lead_source),
    aes(label = paste0("Lead source\n", round(center_pos_mb, 2), " Mb")),
    vjust = -0.7,
    fontface = "bold",
    size = 4
  ) +
  
  labs(
    title = "Regional Association Plot",
    subtitle = paste0("All trans-eQTL sources on chr6 around ", round(center_pos_mb, 2), " Mb"),
    x = "Genomic Position on Chromosome 6 (Mb)",
    y = expression(-log[10](P-value))
  ) +
  
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  scale_y_continuous(
    limits = c(0, 35),
    expand = c(0, 0)
  )

p_lz
ggsave(out_lz, p_lz, width = 8, height = 6)

cat("\nSaved:\n  ", out_lz, "\n", sep = "")
