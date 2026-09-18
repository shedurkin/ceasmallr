10-diff-methyl-DSS
================
Kathleen Durkin
2026-09-10

- [1 Setup](#1-setup)
- [2 Download methylation calls](#2-download-methylation-calls)
- [3 Read Bismark coverage into a BSseq
  object](#3-read-bismark-coverage-into-a-bsseq-object)
- [4 Multi-factor model (design-level
  tests)](#4-multi-factor-model-design-level-tests)
  - [4.1 Parental treatment effect](#41-parental-treatment-effect)
  - [4.2 Stage-dependence of the treatment
    effect](#42-stage-dependence-of-the-treatment-effect)
- [5 Save outputs](#5-save-outputs)
- [6 Two-group smoothed test within each
  stage](#6-two-group-smoothed-test-within-each-stage)
- [7 Save outputs](#7-save-outputs)
  - [7.1 Directional concordance with parental gamete
    DMLs](#71-directional-concordance-with-parental-gamete-dmls)

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

# 2 Download methylation calls

(if necessary)

``` bash
# Run wget to retrieve FastQs and MD5 files
# Note: the --no-clobber command will skip re-downloading any files that are already present in the output directory
wget \
--directory-prefix ../data/bismark-methyl-extraction \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 4 \
--no-host-directories \
--no-parent \
--quiet \
--no-clobber \
--accept="*.deduplicated.bismark.cov.gz,checksums.md5" https://gannet.fish.washington.edu/gitrepos/ceasmallr/output/02.20-bismark-methylation-extraction/
```

``` bash
cd ../data/bismark-methyl-extraction/

echo "How many methylation calling files were downloaded?"
ls *.deduplicated.bismark.cov.gz | wc -l

echo ""

echo "Check file checksums:"
grep '.deduplicated.bismark.cov.gz' checksums.md5 | md5sum -c -
```

We should have N=32 (n=15 ControlxControl crosses, n=17 ExposedxExposed
crosses)

# 3 Read Bismark coverage into a BSseq object

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
# drop <- c("CF01-CM01-Zygote", "CF08-CM04-Larvae", "CF08-CM05-Larvae", "EF04-EM04-Zygote", "EF05-EM05-Zygote")
# file.list <- file.list[!str_detect(basename(file.list), str_c(drop, collapse = "|"))]

# format sample IDs
sample.ids <- basename(file.list) %>%
  str_replace("_R1_001\\.fastp-trim_bismark_bt2_pe\\.deduplicated\\.bismark\\.cov\\.gz$", "") %>%
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
) %>%
  # reference levels: Control (so the coefficient = effect of exposure) and Zygote (earlier stage)
  mutate(treatment = factor(treatment, levels = c("Control", "Exposed")),
         stage     = factor(stage,     levels = c("Zygote",  "Larvae")))

knitr::kable(count(meta, treatment, stage))
```

| treatment | stage  |   n |
|:----------|:-------|----:|
| Control   | Zygote |   8 |
| Control   | Larvae |   7 |
| Exposed   | Zygote |   6 |
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

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |=====                                                                 |   6%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  13%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  32%  |                                                                              |=========================                                             |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |================================                                      |  45%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  65%  |                                                                              |===============================================                       |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  87%  |                                                                              |===============================================================       |  90%  |                                                                              |=================================================================     |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Done in 69 secs

    ## [read.bismark] Parsing files and constructing 'M' and 'Cov' matrices ...

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |=======                                                               |   9%  |                                                                              |=========                                                             |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  19%  |                                                                              |===============                                                       |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  34%  |                                                                              |==========================                                            |  38%  |                                                                              |============================                                          |  41%  |                                                                              |===============================                                       |  44%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  59%  |                                                                              |============================================                          |  62%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  69%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  78%  |                                                                              |=========================================================             |  81%  |                                                                              |===========================================================           |  84%  |                                                                              |=============================================================         |  88%  |                                                                              |===============================================================       |  91%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Done in 26.5 secs

    ## [read.bismark] Constructing BSseq object ...

``` r
BSobj    # positions are the UNION across samples; missing positions get N = 0
```

    ## An object of type 'BSseq' with
    ##   31596324 methylation loci
    ##   32 samples
    ## has not been smoothed
    ## All assays are in-memory

``` r
# DSS tolerates missingness
# For now, require non-zero coverage in at least half the samples of each treatment x stage cell to reduce noise
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

    ## Loci retained after coverage filter: 6337469

EDIT: All of the downstream model fitting calls are *super* intensive
and take forever to run. Since I have to run this from Rscript, and thus
cannot easily retain R objects between runs, this is making script
troubleshooting super time-intensive and annoying. To facilitate
troubleshooting, I’m going to add a test run option that will randomly
subsample a tiny fraction of the data to proceed with. ONLY use this for
test-runs

``` r
# Keeps all samples but randomly thins loci genome-wide
# set test_run to FALSE for full-data runs
test_run  <- FALSE 
test_frac <- 0.02      # keep 2% of loci

if (test_run) {
  set.seed(1)
  idx   <- sort(sample(nrow(BSobj), size = floor(nrow(BSobj) * test_frac)))
  BSobj <- BSobj[idx, ]
  cat("TEST RUN — loci subsampled to:", nrow(BSobj), "\n")
}
```

Also save the set of loci being tested in the following analyses (for
reference when evaluating the effects of including/excluding my
low-coverage offspring samples)

``` r
gr <- granges(BSobj)
tested <- data.frame(chr   = as.character(seqnames(gr)),
                     start = start(gr) - 1L,
                     end   = start(gr))
write_tsv(tested, file.path("../output/10-diff-methyl-DSS", "tested_loci.bed.gz"), col_names = FALSE)
cat("Tested-locus set size:", nrow(tested), "loci\n")
```

    ## Tested-locus set size: 6337469 loci

``` r
# compact summary
knitr::kable(as.data.frame(table(tested$chr)),
             col.names = c("chr", "n_loci_tested"))
```

| chr         | n_loci_tested |
|:------------|--------------:|
| NC_007175.2 |           869 |
| NC_035780.1 |        740251 |
| NC_035781.1 |        688154 |
| NC_035782.1 |        777576 |
| NC_035783.1 |        693889 |
| NC_035784.1 |       1217134 |
| NC_035785.1 |        350838 |
| NC_035786.1 |        386124 |
| NC_035787.1 |        545166 |
| NC_035788.1 |        709957 |
| NC_035789.1 |        227511 |

# 4 Multi-factor model (design-level tests)

Fit an interaction model for tests of *stage-dependence* (does effect of
parental exposure differ between the life stages)

``` r
design <- data.frame(treatment = meta$treatment, stage = meta$stage)

# Interaction model: does the parental-treatment effect differ by stage?
fit_int <- DMLfit.multiFactor(BSobj, design = design,
                              formula = ~ treatment + stage + treatment:stage)
```

    ## Fitting DML model for CpG site: 100000 , 200000 , 300000 , 400000 , 500000 , 600000 , 700000 , 800000 , 900000 , 1000000 , 1100000 , 1200000 , 1300000 , 1400000 , 1500000 , 1600000 , 1700000 , 1800000 , 1900000 , 2000000 , 2100000 , 2200000 , 2300000 , 2400000 , 2500000 , 2600000 , 2700000 , 2800000 , 2900000 , 3000000 , 3100000 , 3200000 , 3300000 , 3400000 , 3500000 , 3600000 , 3700000 , 3800000 , 3900000 , 4000000 , 4100000 , 4200000 , 4300000 , 4400000 , 4500000 , 4600000 , 4700000 , 4800000 , 4900000 , 5000000 , 5100000 , 5200000 , 5300000 , 5400000 , 5500000 , 5600000 , 5700000 , 5800000 , 5900000 , 6000000 , 6100000 , 6200000 , 6300000 ,

``` r
# Inspect coefficient names so we test the right columns
colnames(fit_int$X)
```

    ## [1] "(Intercept)"                  "treatmentExposed"            
    ## [3] "stageLarvae"                  "treatmentExposed:stageLarvae"

``` r
# Expected: "(Intercept)" "treatmentExposed" "stageLarvae" "treatmentExposed:stageLarvae"
```

## 4.1 Parental treatment effect

Fit an additive model (no interaction) to test for signal of parental
treatment effect while *controlling for stage*.

Note that, for a multi-factor analysis, DSS doesn’t seem to rely on a
threshold of % methylation difference? Just an FDR pval.

``` r
fit_add  <- DMLfit.multiFactor(BSobj, design = design, formula = ~ treatment + stage)
```

    ## Fitting DML model for CpG site: 100000 , 200000 , 300000 , 400000 , 500000 , 600000 , 700000 , 800000 , 900000 , 1000000 , 1100000 , 1200000 , 1300000 , 1400000 , 1500000 , 1600000 , 1700000 , 1800000 , 1900000 , 2000000 , 2100000 , 2200000 , 2300000 , 2400000 , 2500000 , 2600000 , 2700000 , 2800000 , 2900000 , 3000000 , 3100000 , 3200000 , 3300000 , 3400000 , 3500000 , 3600000 , 3700000 , 3800000 , 3900000 , 4000000 , 4100000 , 4200000 , 4300000 , 4400000 , 4500000 , 4600000 , 4700000 , 4800000 , 4900000 , 5000000 , 5100000 , 5200000 , 5300000 , 5400000 , 5500000 , 5600000 , 5700000 , 5800000 , 5900000 , 6000000 , 6100000 , 6200000 , 6300000 ,

``` r
test_trt <- DMLtest.multiFactor(fit_add, coef = "treatmentExposed")
head(test_trt[order(test_trt$pvals), ])
```

    ##                 chr      pos      stat        pvals         fdrs
    ## 2758491 NC_035783.1 47063095 -11.36984 5.909464e-30 3.745105e-23
    ## 2758495 NC_035783.1 47063138 -11.19611 4.260543e-29 1.350053e-22
    ## 1609189 NC_035782.1 23797484  10.85410 1.906811e-27 4.028118e-21
    ## 2758493 NC_035783.1 47063100 -10.73816 6.737196e-27 1.067419e-20
    ## 1371039 NC_035781.1 54566055 -10.41559 2.104792e-25 2.667810e-19
    ## 5133966 NC_035787.1 45662859  10.23712 1.351964e-24 1.428005e-18

``` r
# DMLs = FDR threshold (multiFactor gives per-site fdrs, not a methylKit-style % difference)
dml_trt <- test_trt[which(test_trt$fdrs < 0.05), ]
cat("Treatment DMLs (FDR < 0.05):", nrow(dml_trt), "\n")
```

    ## Treatment DMLs (FDR < 0.05): 11356

``` r
# Regions
dmr_trt <- callDMR(test_trt, p.threshold = 0.01, minlen = 50, minCG = 3, dis.merge = 100)
cat("Treatment DMRs:", nrow(dmr_trt), "\n")
```

    ## Treatment DMRs: 1893

## 4.2 Stage-dependence of the treatment effect

A *null* interaction at treatment-associated loci = inherited signal
preserved across stages (H1.2 supported). A *significant* interaction =
the effect differs between zygote and larva, i.e. evidence of
developmental reprogramming.

``` r
test_intx <- DMLtest.multiFactor(fit_int, coef = "treatmentExposed:stageLarvae")
head(test_intx[order(test_intx$pvals), ])
```

    ##                 chr      pos      stat        pvals        fdrs
    ## 103885  NC_035780.1  8153271 -6.341555 2.274570e-10 0.001441502
    ## 5475321 NC_035788.1 13004293  6.101735 1.049235e-09 0.002371815
    ## 3392768 NC_035784.1 41997159 -6.090903 1.122758e-09 0.002371815
    ## 2876716 NC_035783.1 57225695  5.937695 2.890565e-09 0.004579716
    ## 5015998 NC_035787.1 27908837 -5.876962 4.178635e-09 0.005257702
    ## 4215337 NC_035785.1 16140415  5.847915 4.977731e-09 0.005257702

``` r
dml_intx <- test_intx[which(test_intx$fdrs < 0.05), ]
cat("Loci with stage-dependent treatment effect (FDR < 0.05):", nrow(dml_intx), "\n")
```

    ## Loci with stage-dependent treatment effect (FDR < 0.05): 20

# 5 Save outputs

``` r
save_bed <- function(df, path, score_col) {
  # callDMR() returns NULL when no regions pass; callDML() can return 0 rows.
  # Write an empty .bed in that case rather than erroring, so a knit completes.
  if (is.null(df) || nrow(df) == 0) {
    file.create(path)
    message("save_bed: no rows for ", basename(path), " — wrote empty file.")
    return(invisible(NULL))
  }
  # Region tables (callDMR) carry chr/start/end; single-base DML tables
  # (DMLtest/callDML) carry chr/pos — derive a 0-based, half-open interval.
  if (all(c("chr","start","end") %in% names(df))) {
    start <- df$start; end <- df$end
  } else if (all(c("chr","pos") %in% names(df))) {
    start <- df$pos - 1L; end <- df$pos
  } else {
    stop("save_bed: df lacks chr/start/end or chr/pos for ", basename(path),
         " - columns are: ", paste(names(df), collapse=", "))
  }
  score <- if (score_col %in% names(df)) df[[score_col]] else NA_real_
  bed <- data.frame(chr = df$chr, start = start, end = end, score = score)
  write.table(bed, path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# Site-level tables
write_tsv(test_trt,  file.path("../output/10-diff-methyl-DSS", "multifactor_treatment_test.tsv"))
write_tsv(test_intx, file.path("../output/10-diff-methyl-DSS", "multifactor_interaction_test.tsv"))

# Region-level BEDs (multifactor DMR uses `areaStat` as region-level test statistic)
save_bed(dmr_trt, file.path("../output/10-diff-methyl-DSS", "treatment_DMR.bed"), "areaStat")
```

# 6 Two-group smoothed test within each stage

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
  saveRDS(dml, file.path("../output/10-diff-methyl-DSS", paste0(stage_label, "_DMLtest.rds")))
  write_tsv(as.data.frame(dml), file.path("../output/10-diff-methyl-DSS", paste0(stage_label, "_DMLtest.tsv.gz")))
  
  rm(dml, BS_s); gc()
  out
}
res_zyg <- run_stage("Zygote"); gc()
```

    ## Smoothing ...
    ## Estimating dispersion for each CpG site, this will take a while ...

    ## Warning in mclapply(1:nrow(X2), foo, mc.cores = ncores): scheduled cores 127,
    ## 130, 137, 150, 154, 167, 171, 179, 182, 183, 184, 185, 186, 187, 188, 189 did
    ## not deliver results, all values of the jobs will be affected

    ## Warning in shrk.phi[ix] <- shrk.phi2: number of items to replace is not a
    ## multiple of replacement length

    ## Warning in mclapply(1:nrow(X2), foo, mc.cores = ncores): scheduled cores 124,
    ## 128, 138, 152, 153, 154, 156, 157, 158, 161, 162, 163, 165, 167, 169, 170, 171,
    ## 172, 173, 174, 175, 177, 178, 179, 180, 181, 183, 184, 185, 186 did not deliver
    ## results, all values of the jobs will be affected

    ## Warning in shrk.phi[ix] <- shrk.phi2: number of items to replace is not a
    ## multiple of replacement length

    ## Computing test statistics ...

    ##              used   (Mb) gc trigger    (Mb)   max used    (Mb)
    ## Ncells   10574794  564.8   27015671  1442.8   27015671  1442.8
    ## Vcells 1020658173 7787.1 2973613454 22686.9 2973610954 22686.9

``` r
res_lar <- run_stage("Larvae"); gc()
```

    ## Smoothing ...
    ## Estimating dispersion for each CpG site, this will take a while ...
    ## Computing test statistics ...

    ##              used   (Mb) gc trigger    (Mb)   max used    (Mb)
    ## Ncells   10574857  564.8   27015671  1442.8   27015671  1442.8
    ## Vcells 1020746151 7787.7 2973613454 22686.9 2973611614 22686.9

``` r
cat("Zygote: DMLs =", nrow(res_zyg$dml), " DMRs =", nrow(res_zyg$dmr), "\n")
```

    ## Zygote: DMLs = 5859  DMRs = 4393

``` r
cat("Larvae: DMLs =", nrow(res_lar$dml), " DMRs =", nrow(res_lar$dmr), "\n")
```

    ## Larvae: DMLs = 4183  DMRs = 3963

Use existing `.rds` saves if available

``` r
# summarise_stage <- function(stage) {
#   rds <- file.path("../output/10-diff-methyl-DSS", paste0(stage, "_DMLtest.rds"))
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

# 7 Save outputs

``` r
# Per-stage DML BEDs (diff = mu1 - mu2 = exposed - control)
save_bed(res_zyg$dml, "../output/10-diff-methyl-DSS/zygote_DML.bed", "diff")
save_bed(res_lar$dml, "../output/10-diff-methyl-DSS/larvae_DML.bed", "diff")

# Region-level BEDs (diff.Methy is exposed - control at the region level)
save_bed(res_zyg$dmr, "../output/10-diff-methyl-DSS/zygote_DMR.bed", "diff.Methy")
save_bed(res_lar$dmr, "../output/10-diff-methyl-DSS/larvae_DMR.bed", "diff.Methy")
```

## 7.1 Directional concordance with parental gamete DMLs

`07.2` only tested *positional* overlap, but we’re also interested in
directional agreement (hyper/hypo).

As a priliminary check for shared signal, read in the parental beds
(sperm/egg DMLs, column 5 = parental methylation diff), join to the
offspring DMLs on `chr` + `pos`, and compare the sign of the difference.

``` r
zyg_dmls <- read_tsv("../output/10-diff-methyl-DSS/zygote_DML.bed",
                  col_names = c("chr", "start", "end", "zyg_diff"))
```

    ## Rows: 5859 Columns: 4
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (1): chr
    ## dbl (3): start, end, zyg_diff
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
lar_dmls <- read_tsv("../output/10-diff-methyl-DSS/larvae_DML.bed",
                  col_names = c("chr", "start", "end", "lar_diff"))
```

    ## Rows: 4183 Columns: 4
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (1): chr
    ## dbl (3): start, end, lar_diff
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
sperm <- read_tsv("../data/adult_male_dml.bed",
                  col_names = c("chr", "start", "end", "strand", "par_diff")) %>%
  transmute(chr, start, end = end - 1, par_diff = par_diff/100)
```

    ## Rows: 4175 Columns: 5
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (2): chr, strand
    ## dbl (3): start, end, par_diff
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
egg <- read_tsv("../data/adult_female_dml.bed",
                  col_names = c("chr", "start", "end", "strand", "par_diff")) %>%
  transmute(chr, start, end = end - 1, par_diff = par_diff/100)
```

    ## Rows: 128 Columns: 5
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (2): chr, strand
    ## dbl (3): start, end, par_diff
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
shared_zyg_egg <- inner_join(zyg_dmls, egg, by = c("chr", "start", "end")) %>%
  mutate(concordant = sign(zyg_diff) == sign(par_diff))

shared_lar_egg<- inner_join(lar_dmls, egg, by = c("chr", "start", "end")) %>%
  mutate(concordant = sign(lar_diff) == sign(par_diff))

summarise(shared_zyg_egg, n_shared = n(), n_concordant = sum(concordant),
          pct_concordant = mean(concordant) * 100)
```

    ## # A tibble: 1 x 3
    ##   n_shared n_concordant pct_concordant
    ##      <int>        <int>          <dbl>
    ## 1        1            1            100

``` r
summarise(shared_lar_egg, n_shared = n(), n_concordant = sum(concordant),
          pct_concordant = mean(concordant) * 100)
```

    ## # A tibble: 1 x 3
    ##   n_shared n_concordant pct_concordant
    ##      <int>        <int>          <dbl>
    ## 1        2            2            100

``` r
shared_zyg_sperm <- inner_join(zyg_dmls, sperm, by = c("chr", "start", "end")) %>%
  mutate(concordant = sign(zyg_diff) == sign(par_diff))

shared_lar_sperm <- inner_join(lar_dmls, sperm, by = c("chr", "start", "end")) %>%
  mutate(concordant = sign(lar_diff) == sign(par_diff))

summarise(shared_zyg_sperm, n_shared = n(), n_concordant = sum(concordant),
          pct_concordant = mean(concordant) * 100)
```

    ## # A tibble: 1 x 3
    ##   n_shared n_concordant pct_concordant
    ##      <int>        <int>          <dbl>
    ## 1       96           96            100

``` r
summarise(shared_lar_sperm, n_shared = n(), n_concordant = sum(concordant),
          pct_concordant = mean(concordant) * 100)
```

    ## # A tibble: 1 x 3
    ##   n_shared n_concordant pct_concordant
    ##      <int>        <int>          <dbl>
    ## 1       70           70            100

``` r
write_tsv(dplyr::select(shared_zyg_egg, -concordant), "../output/10-diff-methyl-DSS/shared_zygote_egg_DML.bed", col_names=TRUE)
write_tsv(dplyr::select(shared_lar_egg, -concordant), "../output/10-diff-methyl-DSS/shared_larvae_egg_DML.bed", col_names=TRUE)

write_tsv(dplyr::select(shared_zyg_sperm, -concordant), "../output/10-diff-methyl-DSS/shared_zygote_sperm_DML.bed", col_names=TRUE)
write_tsv(dplyr::select(shared_lar_sperm, -concordant), "../output/10-diff-methyl-DSS/shared_larvae_sperm_DML.bed", col_names=TRUE)
```

For a real answer, however, we’d want to (a) re-evaluate parental
differential methylation signal using DSS as well (instead of relying on
the DMLs that were previously identified by Venkataraman et al using
`methylKit`), and (b) test whether parent-offspring DML overlap exceeds
chance (“enriched”)
