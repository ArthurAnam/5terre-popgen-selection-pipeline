# saltiLASSI pre-production methodological freeze

The primary selection scan is saltiLASSI Lambda, with original LASSI T retained as a secondary methodological benchmark.

The LD-derived physical scale is defined population-wise as the minimum distance at which mean r2 falls below one third of the mean r2 for SNP pairs separated by approximately 1 kb:

d_LD,p = min{d: mean_r2_p(d) < (1/3) * mean_r2_p(1 kb)}.

In this workflow, the 1-kb reference is operationalized as the mean for pairs separated by 0.5-1.5 kb and LD is summarized in 1-kb bins. The first below-threshold bin is the primary estimate. A five-consecutive-bin crossing is retained only as a robustness diagnostic and is not used to choose the production window.

The project uses the already frozen within-population MAF>=0.05 panels. This differs from the 2022 empirical CEU/YRI saltiLASSI analysis, which retained all sample-polymorphic biallelic SNPs. Each population's observed windows and genome-wide null HFS use the same ascertainment, so the analysis remains internally coherent; absolute Lambda values are not treated as directly comparable to published CEU/YRI scores.

The local HFS windows are CT 99 SNP, CEU 116, TSI 112 and IBS 116. saltiLASSI step sizes are set to approximately 50% of each window (50, 58, 56 and 58 SNP), transferring the overlap proportion of the published human 201/100 design without treating it as a universal recommendation.

Physical distance is primary. max-extend-bp is fixed at 100000, the default in both lassip v1.1.1 and the pinned v1.2.1.

The null HFS is computed separately for each population as an equal-weight mean of the K-truncated HFS across all supplied autosomal windows. With K=10, lassip requires p_K >= 1/(100*K)=0.001.

Low-mappability regions are excluded using the GRCh37/hg19 CRG100 100-mer alignability track with mean threshold 0.9.

Because CT does not have a sufficiently calibrated recent demographic model for a defensible demography-matched whole-genome neutral threshold, CT uses a pre-specified empirical-outlier framework. The top 1% of genome-wide Lambda is the primary candidate threshold; top 0.1% and top 5% are descriptive. Consecutive above-threshold windows are merged. These regions are called candidate/outlier regions, not simulation-based significant regions.

Original LASSI T remains secondary. CEU is an empirical benchmark because Harris and DeGiorgio 2020 used 117/12 SNP windows, while the independent MAF>=0.05 calibration here yields 116/12. Comparison focuses on candidate loci and ranking/concordance rather than exact T-score replication.
