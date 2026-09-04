06.3-differential-methylation-DSS
================
Kathleen Durkin
2026-08-31

- [1 Setup](#1-setup)
- [2 Read Bismark coverage into a BSseq
  object](#2-read-bismark-coverage-into-a-bsseq-object)
- [3 Multi-factor model (design-level
  tests)](#3-multi-factor-model-design-level-tests)
  - [3.1 Parental treatment effect](#31-parental-treatment-effect)
  - [3.2 Stage-dependence of the treatment
    effect](#32-stage-dependence-of-the-treatment-effect)
- [4 Save outputs](#4-save-outputs)
- [5 Two-group smoothed test within each
  stage](#5-two-group-smoothed-test-within-each-stage)
- [6 Save outputs](#6-save-outputs)
  - [6.1 Directional concordance with parental gamete
    DMLs](#61-directional-concordance-with-parental-gamete-dmls)

Differential methylation analysis using the `DSS` package (Dispersion
Shrinkage for Sequencing data)

Why `DSS` here:

- The experimental design has **two crossed factors** (parental
  treatment: Control/Exposed; life stage: Zygote/Larvae). `DSS` allows
  treatment, stage, and their interaction can be modelled directly.
- `DSS`/`bsseq` retains the **union** of CpG positions and tolerate
  missing / low coverage,rather than requiring coverage in *every*
  sample the way `methylKit::unite()` does.

Caveats to keep in mind:

- `DMLtest.multiFactor()` does **not** spatially smooth (smoothing is
  only in the two-group `DMLtest()`). Below we use the multi-factor fit
  for the design-level tests, and a separate smoothed two-group
  `DMLtest()` within each stage to mirror the `06.2` split.
- `DSS` aggregates DMRs from per-site statistics and gives no
  region-level FDR
- Only fixed-effect covariates are supported (no family/cross random
  effects)

Inputs: Bismark coverage files `[].deduplicated.bismark.cov.gz`, the
same ones used in `06.x` (stored under
`../data/bismark-methyl-extraction/`).

# 1 Setup

I ended up needing to run this code from using a .sh script from the
terminal for two reasons. First, Klone doesn’t seem to have the DSS
package available, and I don’t seem to have the requisite permissions to
install it, so I need to run all of this code within a conda
environment. Second, some of these runs are so comuptationally intensive
that I need \>200G of memory. I therefore needed to request a node with
a really high amount of memory on which to run this code.

Below is a copy of the simple `.sh` script I’m using for the run. From a
Klone login node, I can just use the `sbatch` command to run this `.sh`
script, which will request a suffiently powered node, open a conda
environment, and run all R code in this `.Rmd` file.

``` bash
#!/bin/bash
#SBATCH --account=srlab
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --cpus-per-task=8
#SBATCH --mem=350G
#SBATCH --time=1-00:00:00
#SBATCH --output=dss_%j.out
#SBATCH --error=dss_%j.err
#SBATCH --chdir=/gscratch/srlab/kdurkin1/kathleen-coral/project/code

apptainer exec --bind /gscratch:/gscratch \
 /gscratch/srlab/kdurkin1/srlab-R4.4-bioinformatics-container-703094b.sif \
  bash -c 'source /srlab/programs/miniforge3-24.7.1-0/etc/profile.d/conda.sh && conda activate /gscratch/srlab/kdurkin1/.conda/envs/dss && Rscript -e "rmarkdown::render(\"10.Rmd\")"'
```

``` r
# Install if needed:
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DSS", "bsseq"))

library(DSS)
library(bsseq)
library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
```

# 2 Read Bismark coverage into a BSseq object

Bismark `.cov` columns (no header):
`chr  start  end  meth%  count_methylated  count_unmethylated`. For
`DSS` we need, per sample, a df with columns `chr`, `pos`, `N` (total
coverage), and `X` (methylated count).

(Positions are 1-based in both Bismark `.cov` and `bsseq`, so no offset
is needed)

``` r
file.list  <- list.files("../data/bismark-methyl-extraction", pattern = "\\.bismark\\.cov\\.gz$", full.names = TRUE)

# As in 06.2, drop the low-yield libraries (<100M Cs)
# BUT, could keep all samples and let DSS handle the coverage differences -- come back to later.
drop <- c("CF01-CM01-Zygote", "CF08-CM04-Larvae", "CF08-CM05-Larvae", "EF04-EM04-Zygote", "EF05-EM05-Zygote")

file.list <- file.list[!str_detect(basename(file.list), str_c(drop, collapse = "|"))]

# format sample IDs
sample.ids <- basename(file.list) |>
  str_replace("_R1_001\\.fastp-trim_bismark_bt2_pe\\.deduplicated\\.bismark\\.cov\\.gz$", "") |>
  str_replace("_R1_001\\.fastp-trim\\.REPAIRED_bismark_bt2_pe\\.deduplicated\\.bismark\\.cov\\.gz$", ".REP")

# Grab treatment and lifestage info from Sample IDs
#   first letter of the cross code = parental treatment (C = control, E = exposed)
#   stage encoded in the ID string (Zygote / Larvae)
meta <- tibble(
  sample    = sample.ids,
  file      = file.list,
  treatment = case_when(str_starts(sample.ids, "E") ~ "Exposed",
                        str_starts(sample.ids, "C") ~ "Control",
                        TRUE ~ NA_character_),
  stage     = case_when(str_detect(sample.ids, "Zygote") ~ "Zygote",
                        str_detect(sample.ids, "Larvae") ~ "Larvae",
                        TRUE ~ NA_character_)
) |>
  # reference levels: Control (so the coefficient = effect of exposure) and Zygote (earlier stage)
  mutate(treatment = factor(treatment, levels = c("Control", "Exposed")),
         stage     = factor(stage,     levels = c("Zygote",  "Larvae")))

knitr::kable(count(meta, treatment, stage))
```

| treatment | stage  |   n |
|:----------|:-------|----:|
| Control   | Zygote |   7 |
| Control   | Larvae |   5 |
| Exposed   | Zygote |   4 |
| Exposed   | Larvae |  11 |

``` r
colData <- DataFrame(treatment = meta$treatment,
                     stage     = meta$stage,
                     row.names = meta$sample)

BSobj <- read.bismark(files          = meta$file,
                      colData        = colData,
                      rmZeroCov      = TRUE,          # drops positions covered in 0 samples
                      strandCollapse = FALSE,         # .cov isn't a cytosine report; leave off
                      BPPARAM        = MulticoreParam(workers = 12, progressbar = TRUE),
                      verbose        = TRUE)
```

    ## [read.bismark] Parsing files and constructing valid loci ...

    ##   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |=====                                                                 |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  15%  |                                                                              |=============                                                         |  19%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  35%  |                                                                              |===========================                                           |  38%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  58%  |                                                                              |===========================================                           |  62%  |                                                                              |==============================================                        |  65%  |                                                                              |================================================                      |  69%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |=========================================================             |  81%  |                                                                              |===========================================================           |  85%  |                                                                              |==============================================================        |  88%  |                                                                              |=================================================================     |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%

    ## Done in 56.6 secs

    ## [read.bismark] Parsing files and constructing 'M' and 'Cov' matrices ...

    ##   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  15%  |                                                                              |=============                                                         |  19%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  26%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  33%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  41%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  59%  |                                                                              |============================================                          |  63%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  78%  |                                                                              |=========================================================             |  81%  |                                                                              |============================================================          |  85%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%

    ## Done in 27 secs

    ## [read.bismark] Constructing BSseq object ...

``` r
BSobj    # positions are the UNION across samples; missing positions get N = 0
```

    ## An object of type 'BSseq' with
    ##   31450950 methylation loci
    ##   27 samples
    ## has not been smoothed
    ## All assays are in-memory

``` r
# DSS tolerates missingness, but a light filter helps reduce noise
# For now, require non-zero coverage in at least half the samples of each treatment x stage cell.
# (MAY WANT TO ADJUST LATER)
cov_mat <- getCoverage(BSobj, type = "Cov")           # loci x samples
cells   <- interaction(meta$treatment, meta$stage, drop = TRUE)

keep <- Reduce(`&`, lapply(levels(cells), function(cl) {
  idx <- which(cells == cl)
  rowSums(cov_mat[, idx, drop = FALSE] > 0) >= ceiling(length(idx) / 2)
}))

BSobj <- BSobj[keep, ]
cat("Loci retained after coverage filter:", nrow(BSobj), "\n")
```

    ## Loci retained after coverage filter: 7594836

# 3 Multi-factor model (design-level tests)

Fit an interaction model for tests of *stage-dependence* (does effect of
parental exposure differ between the life stages)

``` r
design <- data.frame(treatment = meta$treatment, stage = meta$stage)

# Interaction model: does the parental-treatment effect differ by stage?
fit_int <- DMLfit.multiFactor(BSobj, design = design,
                              formula = ~ treatment + stage + treatment:stage)
```

    ## Fitting DML model for CpG site: 100000 , 200000 , 300000 , 400000 , 500000 , 600000 , 700000 , 800000 , 900000 , 1000000 , 1100000 , 1200000 , 1300000 , 1400000 , 1500000 , 1600000 , 1700000 , 1800000 , 1900000 , 2000000 , 2100000 , 2200000 , 2300000 , 2400000 , 2500000 , 2600000 , 2700000 , 2800000 , 2900000 , 3000000 , 3100000 , 3200000 , 3300000 , 3400000 , 3500000 , 3600000 , 3700000 , 3800000 , 3900000 , 4000000 , 4100000 , 4200000 , 4300000 , 4400000 , 4500000 , 4600000 , 4700000 , 4800000 , 4900000 , 5000000 , 5100000 , 5200000 , 5300000 , 5400000 , 5500000 , 5600000 , 5700000 , 5800000 , 5900000 , 6000000 , 6100000 , 6200000 , 6300000 , 6400000 , 6500000 , 6600000 , 6700000 , 6800000 , 6900000 , 7000000 , 7100000 , 7200000 , 7300000 , 7400000 , 7500000 ,

``` r
# Inspect coefficient names so we test the right columns
colnames(fit_int$X)
```

    ## [1] "(Intercept)"                  "treatmentExposed"            
    ## [3] "stageLarvae"                  "treatmentExposed:stageLarvae"

``` r
# Expected: "(Intercept)" "treatmentExposed" "stageLarvae" "treatmentExposed:stageLarvae"
```

## 3.1 Parental treatment effect

Fit an additive model (no interaction) to test for signal of parental
treatment effect while *controlling for stage*.

Note that, for a multi-factor analysis, DSS doesn’t seem to rely on a
threshold of % methylation difference? Just an FDR pval.

``` r
fit_add  <- DMLfit.multiFactor(BSobj, design = design, formula = ~ treatment + stage)
```

    ## Fitting DML model for CpG site: 100000 , 200000 , 300000 , 400000 , 500000 , 600000 , 700000 , 800000 , 900000 , 1000000 , 1100000 , 1200000 , 1300000 , 1400000 , 1500000 , 1600000 , 1700000 , 1800000 , 1900000 , 2000000 , 2100000 , 2200000 , 2300000 , 2400000 , 2500000 , 2600000 , 2700000 , 2800000 , 2900000 , 3000000 , 3100000 , 3200000 , 3300000 , 3400000 , 3500000 , 3600000 , 3700000 , 3800000 , 3900000 , 4000000 , 4100000 , 4200000 , 4300000 , 4400000 , 4500000 , 4600000 , 4700000 , 4800000 , 4900000 , 5000000 , 5100000 , 5200000 , 5300000 , 5400000 , 5500000 , 5600000 , 5700000 , 5800000 , 5900000 , 6000000 , 6100000 , 6200000 , 6300000 , 6400000 , 6500000 , 6600000 , 6700000 , 6800000 , 6900000 , 7000000 , 7100000 , 7200000 , 7300000 , 7400000 , 7500000 ,

``` r
test_trt <- DMLtest.multiFactor(fit_add, coef = "treatmentExposed")
head(test_trt[order(test_trt$pvals), ])
```

    ##                 chr      pos      stat        pvals         fdrs
    ## 3252087 NC_035783.1 47063138 -11.20139 4.014000e-29 3.048567e-22
    ## 1613416 NC_035781.1 54566109 -11.07856 1.594096e-28 6.053448e-22
    ## 1904940 NC_035782.1 23797484  10.69016 1.131773e-26 2.757798e-20
    ## 3252085 NC_035783.1 47063100 -10.66700 1.452460e-26 2.757798e-20
    ## 6106198 NC_035787.1 45662859  10.24411 1.257759e-24 1.910495e-18
    ## 514989  NC_035780.1 36481455 -10.21824 1.642982e-24 2.079697e-18

``` r
# DMLs = FDR threshold (multiFactor gives per-site fdrs, not a methylKit-style % difference)
dml_trt <- test_trt[which(test_trt$fdrs < 0.05), ]
cat("Treatment DMLs (FDR < 0.05):", nrow(dml_trt), "\n")
```

    ## Treatment DMLs (FDR < 0.05): 8404

``` r
# Regions
dmr_trt <- callDMR(test_trt, p.threshold = 0.01, minlen = 50, minCG = 3, dis.merge = 100)
cat("Treatment DMRs:", nrow(dmr_trt), "\n")
```

    ## Treatment DMRs: 1788

## 3.2 Stage-dependence of the treatment effect

A *null* interaction at treatment-associated loci = inherited signal
preserved across stages (H1.2 supported). A *significant* interaction =
the effect differs between zygote and larva, i.e. evidence of
developmental reprogramming.

``` r
test_intx <- DMLtest.multiFactor(fit_int, coef = "treatmentExposed:stageLarvae")
head(test_intx[order(test_intx$pvals), ])
```

    ##                 chr      pos      stat        pvals         fdrs
    ## 5454373 NC_035786.1 26848260 -7.041529 1.901412e-12 1.444091e-05
    ## 121720  NC_035780.1  8153271 -6.072932 1.255955e-09 4.769388e-03
    ## 5959197 NC_035787.1 27908837 -5.810520 6.227925e-09 1.576669e-02
    ## 5357534 NC_035786.1 11658552  5.681226 1.337328e-08 2.053254e-02
    ## 6957431 NC_035788.1 66167519  5.638656 1.713825e-08 2.053254e-02
    ## 4963529 NC_035785.1 16140415  5.636172 1.738718e-08 2.053254e-02

``` r
dml_intx <- test_intx[which(test_intx$fdrs < 0.05), ]
cat("Loci with stage-dependent treatment effect (FDR < 0.05):", nrow(dml_intx), "\n")
```

    ## Loci with stage-dependent treatment effect (FDR < 0.05): 14

# 4 Save outputs

``` r
save_bed <- function(df, path, score_col) {
  if (!all(c("chr","start","end") %in% names(df))) {
    stop("save_bed: df lacks chr/start/end for ", basename(path),
         " — columns are: ", paste(names(df), collapse=", "))
  }
  score <- if (score_col %in% names(df)) df[[score_col]] else NA_real_
  bed <- data.frame(chr = df$chr, start = df$start, end = df$end, score = score)
  write.table(bed, path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# Site-level tables
write_tsv(test_trt,  file.path("../output/06.3-differential-methylation-DSS", "multifactor_treatment_test.tsv"))
write_tsv(test_intx, file.path("../output/06.3-differential-methylation-DSS", "multifactor_interaction_test.tsv"))

# Region-level BEDs (diff.Methy is exposed - control at the region level)
save_bed(dmr_trt, file.path("../output/06.3-differential-methylation-DSS", "treatment_DMR.bed"), "diff.Methy")
```

# 5 Two-group smoothed test within each stage

Also want to try within-stage tests with DSS. For non-multifactor tets,
DSS makes use of smoothing and dispersion shrinkage to essentiallys
“borrow” info across neighboring CpG sites. This cvan be useful for
handling low-coverage libraries, which is an existing problem with the
zygot libraries.

It also produces effect estimates of the same type as `methylKit` (%
methylation difference), which could be useful for more direct
comparisons to the conventions used in Rondon et al. 2017 and
Venkataraman et al. 2024. It would also be useful for questions
involving directional concordance (hyper- v hypo-methylation)

WARNING: Both the Zygote and Larvae stage-specific DMLtest() runs are
**very** memory-intensive. I’ve needed to up the node request to ~300G
just to successfully complete them, and it takes a while to finish.

``` r
run_stage <- function(stage_label) {
  s_idx  <- which(meta$stage == stage_label)
  BS_s   <- BSobj[, s_idx]
  g_exp  <- meta$sample[s_idx][meta$treatment[s_idx] == "Exposed"]
  g_ctrl <- meta$sample[s_idx][meta$treatment[s_idx] == "Control"]
  dml <- DMLtest(BS_s, group1 = g_exp, group2 = g_ctrl, smoothing = TRUE)
  out <- list(
    dml = callDML(dml, delta = 0.25, p.threshold = 0.01),
    dmr = callDMR(dml, delta = 0.10, p.threshold = 0.01,
                  minlen = 50, minCG = 3, dis.merge = 100, pct.sig = 0.5)
  )
  # write per-stage results to disk immediately, so a later crash can't lose them
  saveRDS(dml, file.path("../output/06.3-differential-methylation-DSS", paste0(stage_label, "_DMLtest.rds")))
  rm(dml, BS_s); gc()
  out
}
res_zyg <- run_stage("Zygote"); gc()
```

    ## Smoothing ...
    ## Estimating dispersion for each CpG site, this will take a while ...
    ## Computing test statistics ...

    ##             used   (Mb) gc trigger    (Mb)   max used    (Mb)
    ## Ncells  10572881  564.7   25790962  1377.4   25790962  1377.4
    ## Vcells 987298622 7532.5 2997064826 22865.8 2997050377 22865.7

``` r
res_lar <- run_stage("Larvae"); gc()
```

    ## Smoothing ...
    ## Estimating dispersion for each CpG site, this will take a while ...
    ## Computing test statistics ...

    ##             used   (Mb) gc trigger    (Mb)   max used    (Mb)
    ## Ncells  10572945  564.7   25790962  1377.4   25790962  1377.4
    ## Vcells 987406984 7533.4 2877246233 21951.7 3596555354 27439.6

``` r
cat("Zygote: DMLs =", nrow(res_zyg$dml), " DMRs =", nrow(res_zyg$dmr), "\n")
```

    ## Zygote: DMLs = 7287  DMRs = 5208

``` r
cat("Larvae: DMLs =", nrow(res_lar$dml), " DMRs =", nrow(res_lar$dmr), "\n")
```

    ## Larvae: DMLs = 5330  DMRs = 4661

Use existing `.rds` saves if available

``` r
# summarise_stage <- function(stage) {
#   rds <- file.path("../output/06.3-differential-methylation-DSS", paste0(stage, "_DMLtest.rds"))
#   if (!file.exists(rds)) {
#     message(stage, ": ", basename(rds), " not found — skipping")
#     return(NULL)
#   }
# 
#   dml_test <- readRDS(rds)   # the full DMLtest() output (per-CpG stats)
# 
#   # Called DMLs and DMRs (same thresholds as the 06.3 run)
#   dmls <- callDML(dml_test, delta = 0.25, p.threshold = 0.01)
#   dmrs <- callDMR(dml_test, delta = 0.10, p.threshold = 0.01,
#                   minlen = 50, minCG = 3, dis.merge = 100, pct.sig = 0.5)
# 
#   # callDMR returns NULL if no regions pass
#   n_dmr <- if (is.null(dmrs)) 0L else nrow(dmrs)
# 
#   cat("\n===", stage, "===\n")
#   cat("CpG sites tested:      ", nrow(dml_test), "\n")
#   cat("DMLs (|diff|>0.25, p<0.01):", nrow(dmls), "\n")
#   cat("  hyper (exposed > control):", sum(dmls$diff > 0), "\n")
#   cat("  hypo  (exposed < control):", sum(dmls$diff < 0), "\n")
#   cat("DMRs:                  ", n_dmr, "\n")
#   if (n_dmr > 0) {
#     cat("  mean DMR length (bp): ", round(mean(dmrs$length), 1), "\n")
#     cat("  mean CpGs per DMR:    ", round(mean(dmrs$nCG), 1), "\n")
#   }
# 
#   invisible(list(test = dml_test, dml = dmls, dmr = dmrs))
# }
# 
# zyg <- summarise_stage("Zygote")
# lar <- summarise_stage("Larvae")
```

# 6 Save outputs

``` r
# Per-stage DML BEDs (diff = mu1 - mu2 = exposed - control)
save_bed(res_zyg$dml, file.path("../output/06.3-differential-methylation-DSS", "zygote_DML.bed"), "diff")
save_bed(res_lar$dml, file.path("../output/06.3-differential-methylation-DSS", "larvae_DML.bed"), "diff")

# Region-level BEDs (diff.Methy is exposed - control at the region level)
save_bed(res_zyg$dmr, file.path("../output/06.3-differential-methylation-DSS", "zygote_DMR.bed"), "diff.Methy")
save_bed(res_lar$dmr, file.path("../output/06.3-differential-methylation-DSS", "larvae_DMR.bed"), "diff.Methy")
```

## 6.1 Directional concordance with parental gamete DMLs

`07.2` only tested *positional* overlap, but we’re also interested in
directional agreement (hyper/hypo).

As a priliminary check for shared signal, read in the parental beds
(sperm/egg DMLs, column 5 = parental methylation diff), join to the
offspring DMLs on `chr` + `pos`, and compare the sign of the difference.

``` r
offspring <- res_zyg$dml |> transmute(chr, pos, off_diff = diff)

sperm <- read_tsv("../data/adult_male_dml.bed",
                  col_names = c("chr", "start", "end", "strand", "par_diff"),
                  col_types = "ciicd") |>
  transmute(chr, pos = start + 1L, par_diff)       # BED is 0-based; +1 to match 1-based pos

shared <- inner_join(offspring, sperm, by = c("chr", "pos")) |>
  mutate(concordant = sign(off_diff) == sign(par_diff))

summarise(shared, n_shared = n(), n_concordant = sum(concordant),
          pct_concordant = mean(concordant) * 100)
```

    ##   n_shared n_concordant pct_concordant
    ## 1       86           86            100

``` r
# Enrichment vs. chance: is n_shared more than expected given the tested background?
# tested   <- nrow(BSobj)                          # universe of tested CpGs
# n_off    <- nrow(offspring); n_par <- nrow(sperm); k <- nrow(shared)
# phyper(k - 1, n_par, tested - n_par, n_off, lower.tail = FALSE)
```

``` r
offspring <- res_lar$dml |> transmute(chr, pos, off_diff = diff)

shared <- inner_join(offspring, sperm, by = c("chr", "pos")) |>
  mutate(concordant = sign(off_diff) == sign(par_diff))

summarise(shared, n_shared = n(), n_concordant = sum(concordant),
          pct_concordant = mean(concordant) * 100)
```

    ##   n_shared n_concordant pct_concordant
    ## 1       74           73       98.64865

``` r
# Enrichment vs. chance: is n_shared more than expected given the tested background?
# tested   <- nrow(BSobj)                          # universe of tested CpGs
# n_off    <- nrow(offspring); n_par <- nrow(sperm); k <- nrow(shared)
# phyper(k - 1, n_par, tested - n_par, n_off, lower.tail = FALSE)
```

For a real answer, however, we’d want to (a) re-evaluate parental
differential methylation signal using DSS as well (instead of relying on
the DMLs that were previously identified by Venkataraman et al using
`methylKit`), and (b) test whether parent-offspring DML overlap exceeds
chance (“enriched”)
