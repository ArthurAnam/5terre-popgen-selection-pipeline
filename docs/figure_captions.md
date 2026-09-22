# Figure conventions and captions

## Reproducible figure policy

For each scientific figure, the workflow should retain a high-resolution PNG, a vector PDF when supported, and a self-contained interactive HTML whenever the underlying data benefit from hover inspection.

Interactive individual-level figures must expose sample ID and plotted values on hover. Pairwise analyses must expose both sample IDs. Aggregate curves such as LD decay should expose bin coordinate, estimate, and observation count.

The HTML is an inspection aid. Manuscript figures remain the static publication-oriented PNG/PDF outputs, and every conclusion must remain traceable to the underlying versioned tables.

Captions are centralized in config/figure_registry.tsv.

## KING-Robust kinship versus IBS0

Pairwise KING-Robust kinship estimates versus IBS0 among the 46 post-QC Cinque Terre individuals. Each point represents one pair of individuals; hover text reports both sample IDs, IBS0, kinship, and relationship class. Horizontal reference lines mark the exploratory third-degree threshold (0.0442) and the second-degree lower bound (0.0884).

## Joint PCA

Joint smartpca analysis of the 46 Cinque Terre and 503 1000 Genomes EUR individuals using the primary MAF>=0.05, LD-pruned marker panel. Each point represents one individual; interactive hover reports sample ID, population, and PC coordinates. High-LD-mask panels are retained as sensitivity views.

## FROH from ROH >=1.5 Mb

Individual FROH calculated from all autosomal ROH >=1.5 Mb using a 2.77-Gb autosomal denominator. Hover text reports sample ID, population, FROH, total ROH length, and ROH count.

## FROH contributed by ROH >=5 Mb

Individual contribution to FROH from very long ROH >=5 Mb. Hover text identifies each individual and its long-ROH burden.

## ROH count versus total ROH length

Relationship between the number of ROH >=1.5 Mb and their summed physical length for each individual. Hover labels report sample, population, N_ROH, total ROH length, maximum ROH, and FROH.

## Cinque Terre LD decay, MAF >=0.05

Mean pairwise r2 as a function of physical distance in the final QCed Cinque Terre cohort using the MAF>=0.05 marker panel selected for LASSI calibration. LD is summarized in 1-kb bins from deterministic genome-wide anchor SNPs. The mean baseline r2 for 0.5-1.5 kb is 0.4145; one third is 0.1382. The first bin below that threshold is centered at 55.5 kb, with a five-consecutive-bin stability crossing at 57.5 kb.

## Cinque Terre LD decay sensitivity, MAF >=0.01

Prespecified lower-MAF sensitivity for CT LD-decay calibration. The first one-third-baseline crossing is 63.5 kb and the five-bin stable crossing is 67.5 kb. This result is retained as a robustness check only; the LASSI window is calibrated from MAF>=0.05.

## CEU LD decay, MAF >=0.05

Population-specific CEU LD decay used to calibrate the CEU LASSI window. The
mean 0.5-1.5-kb baseline r2 is 0.4126 and the direct/stable one-third-baseline
crossing is 53.5 kb. This physical scale converts to a median 116 SNPs and
winsize/winstep 116/12.

## TSI LD decay, MAF >=0.05

Population-specific TSI LD decay used to calibrate the TSI LASSI window. The
mean 0.5-1.5-kb baseline r2 is 0.4086 and the direct/stable crossing is 52.5 kb.
This converts to a median 112 SNPs and winsize/winstep 112/11.

## IBS LD decay, MAF >=0.05

Population-specific IBS LD decay used to calibrate the IBS LASSI window. The
mean 0.5-1.5-kb baseline r2 is 0.4113 and the direct/stable crossing is 53.5 kb.
This converts to a median 116 SNPs and winsize/winstep 116/12.

## Population-specific LD decay comparison

Comparison of CT, CEU, TSI and IBS MAF>=0.05 LD-decay curves. The selected
one-third-baseline crossings are 55.5, 53.5, 52.5 and 53.5 kb, respectively.
Although the physical decay scales are similar, their empirical SNP-delimited
LASSI windows differ (99, 116, 112 and 116 SNPs), reflecting population-specific
marker density after within-population MAF filtering.


## FIN LD decay, MAF >=0.05

Population-specific FIN LD decay is retained as European LD context. FIN is not
part of the primary LASSI scan, so its decay curve is descriptive and does not
define a production LASSI winsize.

## GBR LD decay, MAF >=0.05

Population-specific GBR LD decay is retained as European LD context. GBR is not
part of the primary LASSI scan, so its decay curve is descriptive and does not
define a production LASSI winsize.

## Full European LD-decay comparison

The combined interactive panel includes CT, CEU, FIN, GBR, IBS and TSI under
the same MAF>=0.05 LD-decay protocol. This gives a complete European reference
context for the CT decay curve while preserving the narrower CT/CEU/TSI/IBS
design for LASSI itself.
