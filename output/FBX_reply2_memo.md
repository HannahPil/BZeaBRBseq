---
title: |
  Zm00001eb375600 (*fbxl1*) — inversion retracted; independent soft-clip
  audit in the CDS window confirms your geometry model
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

**To:** Rubén Rellán-Álvarez  **From:** Hannah Pil  **Re:** Results from your reply memo of 2026-09-14

---

## 1. Summary

Three points.

1. **Analysis 1 inversion claim retracted.** Your §4 geometry model
   reproduces the 3.08× with no expression difference; on medians the
   residual is 0.52–1.00×; my flagged medians-vs-means gap is
   load-bearing, not a caveat. R_div/R_cons is not a bias diagnostic in
   the way I used it. Accepted.
2. **A new independent line falls out cleanly on the geometry side.**
   The soft-clip audit you sketched in §6.2, run on every read overlapping
   the CDS window, shows Teo carriers have a **12× larger mean trailing
   soft-clip length** than B73 (3.82 bp vs 0.31 bp) with a 1.5× higher
   soft-clip *rate* (39% vs 26%). Reads at this locus in Teo carriers
   are systematically getting clipped at their 3′ end — exactly the
   displaced-UTR-read fingerprint. Independent of the R_div/R_cons or
   coverage-profile arguments, and pointing at the same conclusion.
3. **Per-taxon stratification (your §9) — mixed.** Hueh is the extreme
   case (63% clip rate, ~10 bp mean trailing clip) — its analysis-3
   ratio of 2.01× is an artifact of donor divergence, as you suspected.
   But Zlux and Zdip look surprisingly mild (Zdip's mean softR is
   *below* B73 background). The parviglumis and mexicana taxa you
   verified as CDS-conserved also show heavily elevated soft-clip
   rates (Bals 50%, Chal 54%), which restates your §4 point from a
   different angle: geometry hits every donor whose 3′UTR diverges
   from B73, not just the distal ones.

## 2. Analysis 1 — inversion retracted

Accepting your §4 in full. Three specific things I want to record so it
doesn't come back later as a "wait, does that hold" question:

- The 3.08× sits inside the 2.12–4.06× geometry-model range with zero
  expression difference.
- On medians (2.13×), the residual attributable to real expression is
  0.52–1.00× — no upregulation left to explain.
- The soft-clip audit below is independent of the geometry model. It
  says the CDS-window reads themselves are not clean CDS reads; they
  extend into the divergent UTR. So even if geometry weren't the
  explanation, the R_cons number I quoted was built on reads that
  shouldn't have been treated as within-CDS to begin with.

R_div / R_cons as a "bias diagnostic" is retired.

## 3. Analysis 2 — presentational fix, and one keep-open

Reordered by actual stringency:

| Threshold | *n* tests | % teo-lowers |
|---|---:|---:|
| all tests | 41,885 | 51.4% |
| p<0.05    | 9,697  | 58.9% |
| FDR<0.05  | 6,471  | 62.6% |
| p<1e−4    | 3,230  | 70.3% |

Monotone as it should be — you're right that FDR<0.05 ≈ p<0.008 which
sits looser than 1e−4, so ordering by nominal p is the right axis.

Your alternative interpretation (regulatory mismatch produces the same
progressive skew) stays on the table. Agreed on the test — 3′UTR
variant-density regression against β sign, once BC1 delivers per-donor
variants. Scoped as its own memo when that data is in.

## 4. Analysis 4 — soft-clip and start-position audit in the CDS window

Following your §6.2. For every primary mapped read overlapping the CDS
window (chr9:17933693-17933902) across all 384 BAMs, I extracted the
leftmost aligned position and the leading/trailing soft-clip lengths
from the CIGAR string. `samtools view -F 260` filters unmapped and
secondary alignments. 4,303 reads in total; per-sample mean ~12 reads.

**Sample-level summary (Teo carrier vs B73 background):**

| | Teo (*n*=39) | B73 (*n*=252) | ratio |
|---|---:|---:|---:|
| mean reads per sample | 12.3 | 11.6 | 1.06 |
| mean soft-clip rate | **39.3%** | 25.7% | 1.53 |
| mean trailing clip length (bp) | **3.82** | 0.31 | **12.3** |

![Figure R4](FBX_analysis4_softclip/softclip_rate_by_group.png){width=60%}

![Figure R5](FBX_analysis4_softclip/softR_length_density.png){width=80%}

**Reading.** The soft-clip *rate* difference (1.5×) is real but modest.
The soft-clip *length* difference (12×) is the load-bearing signal.
B73 trailing clips concentrate at 1–2 bp — typical alignment noise.
Teo trailing clips are flat over 0–30+ bp — many reads have their 3′
end clipped by tens of bases because it doesn't match the B73
reference. That is the mark of a read whose 5′ end sits in the
conserved CDS but whose 3′ end extends into a divergent 3′UTR.

**Start-position distribution** (Figure R6): both groups pile up in the
last ~50 bp before the 3′ edge, which is the usual BRB-seq 3′ bias.
Teo has proportionally more mass in the upstream 100–200 bp of the
window; not the sharp 3′-edge clustering that a pure "UTR-read
mismapped-into-CDS" story would predict, but consistent with reads
distributed across the CDS whose 3′ ends get clipped depending on how
far they reach into the divergent UTR.

![Figure R6](FBX_analysis4_softclip/start_position_density.png){width=80%}

Files: `output/FBX_analysis4_softclip/per_sample_summary.csv`,
`softclip_rate_by_group.png`, `softR_length_density.png`,
`start_position_density.png`.

## 5. Per-taxon stratification — the §9 check

Same soft-clip audit, stratified by donor taxon. Sample-level means
(so a chatty single sample can't drive a taxon):

| taxon | *n* samples | *n* reads | mean clip rate | mean softR (bp) | CDS conservation |
|---|---:|---:|---:|---:|---|
| **Hueh** | 2 | 17 | **63.3%** | **9.99** | NOT verified (Rubén §9) |
| Chal | 4 | 52 | 54.0% | 4.77 | verified mexicana |
| Bals | 7 | 126 | 50.0% | 5.55 | verified parviglumis |
| Mesa | 6 | 102 | 33.8% | 2.74 | verified mexicana |
| Dura | 10 | 95 | 33.0% | 3.85 | verified mexicana |
| **Zlux** | 4 | 26 | 32.1% | 1.22 | NOT verified |
| Nobo | 4 | 46 | 30.2% | 2.71 | verified mexicana |
| **Zdip** | 2 | 17 | 29.2% | 0.17 | NOT verified |
| B73 ref | 252 | 2922 | 25.7% | 0.31 | reference |

![Figure R7](FBX_analysis4_softclip/softclip_rate_by_taxon.png){width=80%}

![Figure R8](FBX_analysis4_softclip/softR_length_by_taxon.png){width=80%}

**Three things this says:**

1. **Hueh (huehuetenangensis) is what your §9 worried about.** 63%
   clip rate, ~10 bp mean trailing clip, half of clips >5 bp. Only 17
   reads across 2 samples — big uncertainty on the point estimate —
   but the effect size is stark. Its analysis-3 CDS-ratio of 2.01× is
   an artifact of donor divergence, not a real expression signal.
2. **Verified-conserved taxa also show heavy soft-clipping** — Bals
   50%, Chal 54%, mean softR 5–6 bp. Even in donors with 100% CDS
   sequence identity, the CDS window's proximity to a divergent 3′UTR
   means many reads extending across the boundary get clipped. This
   is the geometry point from your §4, arriving from a different
   direction.
3. **Zlux and Zdip look B73-like on the soft-clip metrics.** Zlux mean
   softR 1.22 bp, Zdip 0.17 bp (below B73 background). Their CDS
   sequences are probably not massively divergent from B73, whatever
   else is going on. Their analysis-3 ratios (Zlux 0.61, Zdip 1.26)
   may reflect real but small expression differences rather than
   mapping artifact. Pulling the PanAnd terminal-CDS is still worth
   doing to settle them formally, but no longer urgent.

Files: `output/FBX_analysis4_softclip/per_taxon_summary.csv`,
`softclip_rate_by_taxon.png`, `softR_length_by_taxon.png`.

## 6. Where the fbxl1 story sits now

- The 3.08× "Teo carriers upregulate *fbxl1*" claim is dead.
- The 0.45× 3′UTR-window depletion is real but confounded with
  divergent-region mapping loss, so its magnitude is not
  interpretable as an expression difference either.
- Direction and magnitude of any real *fbxl1* expression effect at
  this locus remain unresolved.
- H2 stays in its original form — no H2′ reframe.
- The A4 whole-gene panel stays, with the "magnitude confounded by 3′
  coverage bias" caveat you suggested rather than being withdrawn.

## 7. One thing worth flagging in my analysis 3 script

My analysis-3 taxon-stratified CDS-ratio plot has a broken lookup
between my lab-shorthand taxa labels (`Bals`, `Mesa`, `Chal`, `Dura`,
`Nobo`) and your donor 3′UTR-length table (`parviglumis`, `mexicana`).
The `donor_utr_bp` column comes back NA for those five taxa, so I
never actually ran the "does the ratio track 3′UTR length?"
falsification test you specified in §6.1. Fixable once I have (or
build) the lab-code → taxon → assembly mapping. The taxon-stratified
ratios themselves (0.61–2.01×) are correct — it's just the
correlation-with-UTR-length step that isn't executed yet.

## 8. Revised revised next steps

Following your §10 ordering:

1. Per-base coverage profile (§5) — **DONE** (analysis 3). Shape-only
   panel shows Teo has proportionally more reads in the upstream CDS
   region than B73, consistent with geometry.
2. Taxon-stratified CDS-window ratio (§6.1) — **PARTIALLY DONE**;
   ratios computed, correlation with donor UTR length not yet run
   because of §7 above.
3. TE annotation of B73 intron-4 insertion — still to do.
4. Terminal CDS from distal-taxon assemblies (§9) — soft-clip audit
   promotes Hueh to high priority; Zdip and Zlux now lower priority
   but worth including in one pass.
5. Formal Wilcoxon per taxon on medians (§6.3) — with the geometry
   story now well-supported the p-values are bounded regardless of
   significance, but happy to run this if the reviewer axis matters.
6. Genome-wide R_div/R_cons scan — TBD, its own memo.
7. UMI dedup + TMM — write-up time.

## 9. What I did not do

- Fix the taxon-to-UTR-length lookup for the §6.1 falsification test.
- Wilcoxon on medians (§6.3).
- PanAnd terminal-CDS pull.
- TE annotation of the intron-4 insertion.

## Attachments

- `output/FBX_analysis4_softclip/per_sample_summary.csv`
- `output/FBX_analysis4_softclip/per_taxon_summary.csv`
- Figures R4–R8 as above
- Scripts: `scripts/FBX_analysis4_softclip_hpc.sh`,
  `scripts/FBX_analysis4_softclip_plot.R`
