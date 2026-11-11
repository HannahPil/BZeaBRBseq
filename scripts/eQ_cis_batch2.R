# ==============================================================================
# eQ_cis_batch2.R
#
# cis-eQTL scan on one or more sequencing batches, with Fausto's calls
# (ZeaL_QTL data/zeal/zeal_rtiger_mosaic.rds). Model, expression transform and
# segment construction are the rtiger run of eQ_cis_genotype_compare.R:
#   - log2(CPM+1) ~ introgression + plate, MatrixEQTL, cisDist = 10
#   - samples matched on pedigree == metadata genotype; B73 checks, blanks and
#     Purple Checks drop out because they are not in the RTIGER set
#   - genotype strings used by >1 sample in the analysed set are excluded (one call per
#     pedigree cannot be assigned to a specific plant)
#   - state 1 (het) and 2 (teosinte) both = introgressed (1)
# Batch comes from metadata$batch. BATCHES env var picks them (default "2";
# "1 2" = combined). Plate is the covariate; it nests batch, so batch needs no term.
#
# Input:   data/processed/Zea_mays_counts.txt (8-plate merge, 2026-11-11)
# Outputs: output/cis_eQTL_<tag>/<tag>_cis_results.tsv, manhattan_<tag>.png
#          tag = batch2 or batch1_2
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MatrixEQTL)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggrepel)
})

rtiger_path <- Sys.getenv("RTIGER_PATH",
  unset = "C:/Users/Hannah Pil/Documents/gemmalab/BZea/ZeaL_QTL/data/zeal/zeal_rtiger_mosaic.rds")
batches <- strsplit(trimws(Sys.getenv("BATCHES", unset = "2")), "[ ,]+")[[1]]
tag     <- paste0("batch", paste(batches, collapse = "_"))
out_dir <- file.path("output", paste0("cis_eQTL_", tag))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
PPL2 <- "Zm00001eb012750"
N_LABEL <- 5

# --- inputs -------------------------------------------------------------------
counts   <- read.delim(file.path("data", "processed", "Zea_mays_counts.txt"),
                       check.names = FALSE, row.names = 1)
metadata <- read.csv(file.path("data", "metadata.csv"), stringsAsFactors = FALSE)
fai      <- read.delim(file.path("data", "external", "Zea_mays.fasta.fai"), header = FALSE)[, 1:2]
chr_len  <- setNames(fai$V2, fai$V1)

gtf_df <- as.data.frame(import(file.path("data", "external", "Zea_mays.gtf")))
gene_coords <- gtf_df |>
  filter(!is.na(gene_id)) |>
  group_by(gene_id) |>
  summarise(chr = as.character(seqnames[1]), start = min(start), end = max(end), .groups = "drop") |>
  filter(gene_id %in% rownames(counts)) |>
  arrange(gene_id)
genes_gr <- GRanges(gene_coords$chr, IRanges(gene_coords$start, gene_coords$end),
                    gene_id = gene_coords$gene_id)

rt <- readRDS(rtiger_path)
mk <- as.data.frame(rt$markers)

b2 <- metadata |> filter(batch %in% batches, sample_id %in% colnames(counts))
reused <- names(which(table(b2$genotype) > 1))
samples <- b2 |>
  filter(genotype %in% colnames(rt$state), !(genotype %in% reused)) |>
  transmute(sample_id, genotype, plate = factor(plate)) |>
  arrange(sample_id)
cat(sprintf("Batch %s samples in counts: %d; in RTIGER and unique pedigree: %d\n",
            paste(batches, collapse = "+"), nrow(b2), nrow(samples)))
print(table(samples$plate))

# --- introgression segments ----------------------------------------------------
rt_segs_for <- function(ped) {
  s <- rt$state[, ped]; out <- list()
  for (ch in unique(mk$chr)) {
    i <- which(mk$chr == ch); pos <- mk$pos[i]; v <- !is.na(s[i]) & s[i] > 0
    if (!any(v)) next
    r <- rle(v); e <- cumsum(r$lengths); b <- e - r$lengths + 1; k <- which(r$values)
    cname <- paste0("chr", ch)
    st <- ifelse(b[k] == 1, 1L, as.integer(floor((pos[pmax(b[k] - 1, 1)] + pos[b[k]]) / 2) + 1))
    en <- ifelse(e[k] == length(pos), as.integer(chr_len[cname]),
                 as.integer(floor((pos[e[k]] + pos[pmin(e[k] + 1, length(pos))]) / 2)))
    out[[length(out) + 1]] <- data.frame(chr = cname, start = st, end = en)
  }
  if (length(out)) do.call(rbind, out) else data.frame(chr = character(), start = integer(), end = integer())
}
segs <- map_dfr(seq_len(nrow(samples)), function(i) {
  s <- rt_segs_for(samples$genotype[i])
  if (nrow(s)) cbind(sample_id = samples$sample_id[i], s) else NULL
})

# --- MatrixEQTL ----------------------------------------------------------------
cnt <- counts[, samples$sample_id]
sg <- GRanges(segs$chr, IRanges(segs$start, segs$end), sample_id = segs$sample_id)
h <- findOverlaps(sg, genes_gr)
geno <- matrix(0, nrow(gene_coords), nrow(samples),
               dimnames = list(gene_coords$gene_id, samples$sample_id))
geno[cbind(genes_gr$gene_id[subjectHits(h)], sg$sample_id[queryHits(h)])] <- 1
geno <- geno[rowSums(geno) > 0 & rowSums(geno) < ncol(geno), , drop = FALSE]

gc  <- gene_coords[match(rownames(geno), gene_coords$gene_id), ]
lib <- colSums(cnt)
ex  <- log2(t(t(as.matrix(cnt[rownames(geno), ])) / lib) * 1e6 + 1)

gs <- SlicedData$new(); gs$CreateFromMatrix(ex);   gs$ResliceCombined(sliceSize = 2000)
ss <- SlicedData$new(); ss$CreateFromMatrix(geno); ss$ResliceCombined(sliceSize = 2000)
cv <- t(model.matrix(~ plate, data = samples)[, -1, drop = FALSE])
cs <- SlicedData$new(); cs$CreateFromMatrix(cv)

out_file <- file.path(out_dir, paste0(tag, "_cis_results.tsv"))
Matrix_eQTL_main(
  snps = ss, gene = gs, cvrt = cs, output_file_name = NULL, pvOutputThreshold = 0,
  useModel = modelLINEAR, errorCovariance = numeric(), verbose = FALSE,
  output_file_name.cis = out_file, pvOutputThreshold.cis = 1,
  snpspos = setNames(as.data.frame(gc[, c("gene_id", "chr", "start")]), c("snp", "chr", "pos")),
  genepos = setNames(as.data.frame(gc[, c("gene_id", "chr", "start", "end")]), c("geneid", "chr", "left", "right")),
  cisDist = 10, pvalue.hist = FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE)

res <- read.delim(out_file, stringsAsFactors = FALSE)
res$carriers <- rowSums(geno)[res$SNP]
write.table(res, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
self <- res[res$SNP == res$gene, ]
cat(sprintf("genes tested=%d  pairs=%d  FDR<0.05: %d pairs / %d self-gene\n",
            nrow(geno), nrow(res), sum(res$FDR < 0.05), sum(self$FDR < 0.05)))
x <- self[self$gene == PPL2, ]
if (nrow(x)) cat(sprintf("ppl2: carriers=%d  beta=%.2f  t=%.1f  FDR=%.2g  rank=%d of %d\n",
                         x$carriers, x$beta, x$t.stat, x$FDR,
                         rank(self$p.value)[self$gene == PPL2], nrow(self))) else
  cat("ppl2 not tested\n")

# --- Manhattan (as eQ_cis_manhattan_compare.R) ---------------------------------
fd <- read.delim(gzfile(file.path("data", "external", "Zm00001eb.1.fulldata.txt.gz")),
                 header = FALSE, quote = "", stringsAsFactors = FALSE)
symbols <- fd |>
  filter(V12 != "-", V12 != "", V11 != "-") |>
  group_by(gene_id = V2) |>
  summarise(symbol = paste(sort(unique(V11)), collapse = "/"), .groups = "drop")
symbols <- bind_rows(symbols |> filter(gene_id != PPL2), tibble(gene_id = PPL2, symbol = "ppl2"))

gene_pos <- gtf_df |>
  filter(!is.na(gene_id)) |>
  group_by(gene_id) |>
  summarise(chr = as.character(seqnames[1]), pos = min(start), .groups = "drop") |>
  mutate(chr_num = suppressWarnings(as.numeric(gsub("chr", "", chr, ignore.case = TRUE)))) |>
  filter(!is.na(chr_num))
offsets <- gene_pos |>
  group_by(chr_num) |> summarise(max_bp = as.numeric(max(pos)), .groups = "drop") |>
  arrange(chr_num) |> mutate(bp_add = lag(cumsum(max_bp), default = 0))

d <- res |>
  inner_join(gene_pos, by = c("SNP" = "gene_id")) |>
  inner_join(offsets, by = "chr_num") |>
  mutate(bp_cum = pos + bp_add, log_p = -log10(p.value),
         color_group = case_when(FDR < 0.05 & chr_num %% 2 == 1 ~ "Sig_Odd",
                                 FDR < 0.05                     ~ "Sig_Even",
                                 chr_num %% 2 == 1              ~ "Base_Odd",
                                 TRUE                           ~ "Base_Even"))
axis_set <- d |> group_by(chr_num) |>
  summarise(center = (max(bp_cum) + min(bp_cum)) / 2, .groups = "drop")
cutoff <- -log10(max(d$p.value[d$FDR < 0.05]))
lab <- d |> filter(SNP == gene) |> arrange(p.value) |> slice_head(n = N_LABEL)
lab <- bind_rows(lab, d |> filter(SNP == PPL2, gene == PPL2, !(gene %in% lab$gene))) |>
  left_join(symbols, by = c("gene" = "gene_id")) |>
  mutate(label = if_else(is.na(symbol), gene, paste0(symbol, " (", gene, ")")))

p <- ggplot(d, aes(bp_cum, log_p, colour = color_group)) +
  geom_point(alpha = 0.75, size = 1.1) +
  geom_hline(yintercept = cutoff, colour = "red", linetype = "dashed") +
  geom_text_repel(data = lab, aes(label = label), colour = "black", size = 4.6,
                  min.segment.length = 0, segment.colour = "grey40", segment.size = 0.3,
                  box.padding = 0.5, max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = c(Base_Odd = "grey45", Base_Even = "grey70",
                                 Sig_Odd = "firebrick", Sig_Even = "darkred")) +
  scale_x_continuous(breaks = axis_set$center, labels = axis_set$chr_num) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(d$log_p[is.finite(d$log_p)]) * 1.08)) +
  labs(title = sprintf("cis-eQTL Manhattan plot, %s (plates %s), RTIGER calls",
                       if (length(batches) > 1) paste("batches", paste(batches, collapse = " and ")) else paste("batch", batches),
                       paste(sort(unique(as.integer(as.character(samples$plate)))), collapse = ", ")),
       subtitle = sprintf(paste0("MatrixEQTL log2(CPM+1) ~ introgression + plate; cisDist = 10 bp; ",
                                 "n = %d samples; dashed line = BH FDR 0.05"), nrow(samples)),
       caption = sprintf("%d gene pairs significant at FDR < 0.05. Labelled: the %d lowest-p genes and ppl2 (Zm00001eb012750).",
                         sum(d$FDR < 0.05), N_LABEL),
       x = "Chromosome", y = "-log10(p)") +
  theme_bw(base_size = 17) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#ffffff", colour = NA),
        plot.title.position = "plot", plot.caption.position = "plot",
        plot.subtitle = element_text(size = 14),
        plot.caption = element_text(hjust = 0, size = 13))

ggsave(file.path(out_dir, paste0("manhattan_", tag, ".png")), p, width = 14, height = 7, dpi = 200, bg = "#ffffff")
cat("Wrote", file.path(out_dir, paste0("manhattan_", tag, ".png")), "\n")
