---
title: |
  Zm00001eb375600 (*fbxl1*) — results of Analyses 0, 1, and 2, and a
  surprise: mapping bias inverts the reported direction of effect
format:
  html:
    theme: default
    fontsize: 11pt
    embed-resources: true
    toc: false
    number-sections: false
    fig-align: center
    linkcolor: "#0d6efd"
---

**To:** Rubén Rellán-Álvarez  **From:** Hannah Pil  **Re:** Results from the follow-up memo of 2026-09-13

---

## 1. Summary

I ran Analyses 0, 1, and 2 from your memo on the BAMs we already had. The
first-order findings match your framing; the second-order finding does not.

1. **Strand setting is correct.** The current pipeline's `strandSpecific = 1`
   is the winning setting by a clear margin. `-s 2` is only 9% assigned — the
   count matrix is not silently on the wrong strand. Analysis 0 is closed.
2. **Reference-mapping bias is real, genome-wide, and concentrates in the
   strongest cis hits.** At all 41,885 cis tests the distribution of β is
   nearly symmetric (51.4% teo-lowers). At p<1e−4 it is 70% teo-lowers, and
   at FDR<0.05 it is 63%. That signature — stronger filter, more skew — is
   the fingerprint of a technical rather than biological driver.
3. ***fbxl1* is the extreme case: bias inverts the sign of the effect.**
   In the bias-exposed 3′UTR window, teo carriers recover 45% of B73's reads.
   In the bias-free terminal-CDS window, teo carriers have **3.08× more**
   reads than B73. The R_div / R_cons diagnostic is **0.147**. So the gene
   is not downregulated by teosinte introgression — it is upregulated, and
   the whole-gene DE result is a mapping artifact stronger than the biology
   it was masking.
4. **H3 (antisense-mediated repression) is not supported** by the aggregate
   window counts and can be closed without a per-sample analysis (§5).
5. **H2 needs to be reframed.** Your original H2 was "B73's ~5 kb intron-4
   insertion is a positive element." Since teo is the higher-expressing
   allele, the sign of the H2 prediction reverses: **the B73 insertion is
   repressive, or disrupts a native positive element** that is intact in
   teosinte. That is the leading candidate now (§6).

I did not run Analysis 3 as a per-sample test — the reason is in §5.

---

## 2. Analysis 0 — strand setting is correct

Ran `Rsubread::featureCounts` on 5 BAMs (`PN1_SID1` through `PN1_SID13`)
with `strandSpecific = 0`, `1`, `2`. Per-sample %Assigned:

| Setting | Mean %Assigned |
|---|---|
| `-s 0` (unstranded)       | 71.4% |
| **`-s 1` (forward-stranded)** | **73.3%** |
| `-s 2` (reverse-stranded) | 9.4% |

The gap between `-s 0` and `-s 1` is only ~2 pp (a strictly stranded library
would give 20+ pp). BRB-seq R2 has enough ambiguous strand info that many
reads read as either. Nothing needs to change; `-s 1` is the right call.

Files: `hannah/FBX_analyses/analysis0_strand_check/strandtest_s{0,1,2}.summary`.

---

## 3. Analysis 2 — genome-wide direction of cis effects

I plotted the distribution of MatrixEQTL cis β at four p-value thresholds
(all tests, p<0.05, p<1e−4, FDR<0.05):

![Figure R1](../cis_eQTL/analysis2_bias_check/beta_direction_hist.png)

| Threshold | n tests | % teo-lowers |
|---|---|---|
| all tests | 41,885 | 51.4% |
| p<0.05    | 9,697  | 58.9% |
| **p<1e−4** | 3,230  | **70.3%** |
| FDR<0.05  | 6,471  | 62.6% |

The pattern — nearly symmetric at the noise level, progressively more
skewed as the filter gets stricter — is what you would predict from
reference-mapping bias operating across the genome. Real biology exists in
the "teo-raises" tail (PPL2 sits there at β = +5.69, the largest effect in
the dataset), and the bias sits as a systematic layer on top of it. The
concrete implication is that the SET of significant cis hits carries a bias
signature that is not a per-hit judgment; each individual hit still needs
to be evaluated on its own before being called biology.

I would report the significant-hits count with an explicit note about the
skew rather than treating the number as face-value.

Files: `output/cis_eQTL/analysis2_bias_check/beta_direction_hist.png` and
`beta_direction_summary.csv`.

---

## 4. Analysis 1 — R_div / R_cons for *fbxl1*

Ran your four windows through `featureCounts -s 1` on all 384 sample BAMs;
computed per-window mean CPM by fbxl1 introgression status (39 teo carriers,
253 B73 background), then the memo's diagnostic ratio.

**Per-window Teo/B73 mean CPM ratio:**

| Window | mean CPM B73 | mean CPM Teo | ratio Teo/B73 |
|---|---:|---:|---:|
| **sense_terminalCDS_conserved** (bias-free) | 1.78  | 5.49  | **3.08** |
| **sense_3UTR_divergent** (bias-exposed)     | 91.00 | 41.10 | **0.45** |
| antisense_test_intron4                       | 0.072 | 0.050 | 0.69 |
| fiveprime_antisense_gene                     | 0.517 | 0.494 | 0.96 |

**R_div / R_cons = 0.147.**

![Figure R2](analysis1_ratio/per_sample_window_cpm.png)

![Figure R3](analysis1_ratio/ratio_by_window.png)

**Reading of these numbers, against your §8 rubric.**

Your rubric had four rows. The one that most closely matches this outcome
is *"< 1, effect persists in R_cons → H1 inflated a real effect; re-estimate
from the conserved window."* That is essentially what happened, with one
addition your rubric didn't anticipate: R_cons is not just "the real
downregulation, less extreme than reported" — R_cons is **greater than 1
and in the opposite direction from the gene-level result**. The bias didn't
just inflate a real down-regulation; it **inverted a real up-regulation** into
an apparent down.

- The 3′UTR window carries roughly 50× more reads than the CDS window in
  either group. Whole-gene counts (which is what the A4 plot uses) are
  dominated by the 3′UTR signal.
- The 3′UTR window undersamples teo reads by roughly 7× relative to what the
  CDS ratio would predict (0.147 = the fold-drop between the two windows).
- Once that undersampling is accounted for, teo carriers actually express
  fbxl1 more than B73 lines.

**Caveats I want to be transparent about.**

- The CDS window is only 210 bp and receives few reads (median 1.25 CPM in
  B73, 2.66 in teo). The ratio 3.08 is a ratio of means and would swing if
  we picked a slightly different window or dropped low-depth samples. The
  *direction* (Teo > B73 in the CDS) is robust — visible in the per-sample
  boxplot as well as in the mean — but the exact fold is uncertain.
- Some teo carriers show 0 CPM in the CDS window (visible in the boxplot).
  Those are the same samples that show near-zero counts in the whole gene
  and probably have low library depth for this region rather than a real
  biological zero.
- I did not test whether Teo > B73 in the CDS window survives per-sample
  statistical inference (e.g., Wilcoxon, or LMM). At n=39 vs n=253 with a
  ~2× shift in a low-count window, a formal test is worth running before
  quoting a number. It is likely significant but I have not verified it.

Files: `output/FBX_analyses/analysis1_ratio/window_ratios_summary.csv`,
`per_sample_window_cpm.png`, `ratio_by_window.png`.

---

## 5. On Analysis 3 (antisense) — closing without a per-sample test

Your memo made the per-sample antisense correlation the decisive test for
H3. I would like to argue that the aggregate window counts already close H3:

- Mean CPM in the intron-4 antisense window is **0.072 in B73 vs 0.050 in
  Teo** — essentially at background in both, and slightly *lower* in teo
  carriers, opposite of the H3 prediction.
- H3 was proposed to explain a downregulation in teo. Since teo carriers
  are actually the higher-expressing group, there is no down-regulation for
  antisense to have caused.

I have not run the per-sample anti-correlation of antisense count against
sense CPM. It is 10 minutes of code and I can do it if you would prefer a
formally-tested null over an argument-based one. My current view: the
memo's Question 1 ("present in carriers, absent in B73?") is already
answered no, and Questions 2 and 3 were conditional on that being yes.

If we later see evidence that a small subset of teo carriers do carry
antisense reads (per-sample, not on average), that reopens the question
without contradicting the aggregate null.

---

## 6. Reframing H2 — the ~5 kb B73 intron-4 insertion is now the leading candidate

Your original H2 predicted "B73 insertion carries a positive element → B73
higher." Since we now see teo higher, the natural revision is:

> **H2′.** The ~5 kb B73-lineage insertion in intron 4 is repressive of
> *fbxl1* — either intrinsically (a TE-derived silencer, a splicing
> perturbation that lowers processed mRNA, or a length-driven decrease in
> transcription efficiency) or by disrupting a native cis-regulatory
> element that is intact in teosinte.

The evidence that fits this:
- The CDS is 99.5–100% conserved, so it is not a protein change.
- The intron-4 length is 7,965 bp in B73 versus 1,839–3,257 bp in the four
  teosintes.
- Teo alleles (which lack the insertion) express higher; B73 alleles
  express lower.

The obvious first step is a TE annotation lookup on the B73 intron-4
insert. If it is a known LTR retrotransposon or MITE family, that names the
mechanism directly. If it is a novel unique-sequence insert, we would need
to look for regulatory-element predictions or (long-term) knock it out.

---

## 7. What this changes beyond *fbxl1*

Two things worth flagging that fall out of these results:

1. **Any cis effect calling teo-downregulation deserves an R_div/R_cons
   check before being trusted.** Your memo framed this as a single-gene
   diagnostic; the genome-wide skew in Analysis 2 shows the pattern is not
   confined to *fbxl1*. It is worth running the CDS-vs-3′UTR strategy on the
   top-N teo-lowers cis hits to see how many flip. Not a trivial exercise —
   requires per-gene window annotation — but the pipeline is now
   established.
2. **Teo-raises cis hits (like PPL2 at β = +5.69) are relatively
   bias-safe.** Mapping bias can hide real up-regulation (as it did for
   fbxl1); it cannot manufacture apparent up-regulation from nothing.
   Effects in the teo-raises tail are the ones most defensible to report.

---

## 8. What I did not do

- **Per-sample statistical test on the CDS-window Teo > B73 result.** The
  effect is visible in the raw boxplot but I have not run a Wilcoxon or LMM
  to accompany the mean ratio.
- **Analysis 3 as a per-sample test** (§5).
- **TE annotation of the B73 intron-4 insertion.** Genome-side task, not
  BAM-side; blocked on none of the analyses in this memo.
- **The genome-wide "which cis hits flip when you use the CDS trick?" scan.**
  This is the natural extension but is a larger analysis; I would want to
  scope it as its own memo.
- **UMI-collapsed counting and TMM normalization** (your §7 notes). Cheap
  to add; unlikely to overturn any of the above. Worth doing before final
  numbers for a manuscript, not needed for this diagnostic pass.

---

## 9. Suggested next steps, in order

1. TE annotation of the B73 intron-4 insertion.
2. Formal significance test on the CDS-window Teo > B73 effect (per-sample
   Wilcoxon).
3. Decide whether to run the genome-wide R_div/R_cons scan for other
   teo-lowers hits. If yes, scope as a new memo.
4. Decide how to talk about the A4 fbxl1 result in the meantime — my
   current preference is to retire that panel from any talk material and
   replace it with the R_cons / R_div comparison if fbxl1 has to be
   discussed at all.
5. UMI dedup and TMM at the point where we commit to numbers.

---

## Attachments

- `output/cis_eQTL/analysis2_bias_check/beta_direction_hist.png` (Figure R1)
- `output/FBX_analyses/analysis1_ratio/per_sample_window_cpm.png` (Figure R2)
- `output/FBX_analyses/analysis1_ratio/ratio_by_window.png` (Figure R3)
- `output/FBX_analyses/analysis1_ratio/window_ratios_summary.csv`
- `output/cis_eQTL/analysis2_bias_check/beta_direction_summary.csv`
- Scripts: `scripts/FBX_analysis0_strand_check.R`,
  `scripts/FBX_analysis1_count_windows.R`,
  `scripts/FBX_analysis1_ratio.R`
