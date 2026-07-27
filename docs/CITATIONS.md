# Citations — RNA-seq Pipeline

Please cite the tools and methods you use. Core references for this pipeline:

## Quantification & alignment

- **Salmon**  
  Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C.  
  *Salmon provides fast and bias-aware quantification of transcript expression.*  
  Nature Methods. 2017;14:417–419.  
  https://doi.org/10.1038/nmeth.4197

- **tximport**  
  Soneson C, Love MI, Robinson MD.  
  *Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences.*  
  F1000Research. 2015;4:1521.  
  https://doi.org/10.12688/f1000research.7563.2

- **STAR**  
  Dobin A, Davis CA, Schlesinger F, et al.  
  *STAR: ultrafast universal RNA-seq aligner.*  
  Bioinformatics. 2013;29(1):15–21.  
  https://doi.org/10.1093/bioinformatics/bts635

- **featureCounts**  
  Liao Y, Smyth GK, Shi W.  
  *featureCounts: an efficient general purpose program for assigning sequence reads to genomic features.*  
  Bioinformatics. 2014;30(7):923–930.  
  https://doi.org/10.1093/bioinformatics/btt656

## Differential expression

- **DESeq2**  
  Love MI, Huber W, Anders S.  
  *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.*  
  Genome Biology. 2014;15:550.  
  https://doi.org/10.1186/s13059-014-0550-8

- **apeglm** (LFC shrinkage)  
  Zhu A, Ibrahim JG, Love MI.  
  *Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences.*  
  Bioinformatics. 2019;35(12):2084–2092.  
  https://doi.org/10.1093/bioinformatics/bty895

## Unwanted variation & batch correction

- **RUVSeq**  
  Risso D, Ngai J, Speed TP, Dudoit S.  
  *Normalization of RNA-seq data using factor analysis of control genes or samples.*  
  Nature Biotechnology. 2014;32:896–902.  
  https://doi.org/10.1038/nbt.2931

- **SVA / svaseq**  
  Leek JT, Storey JD. *Capturing heterogeneity in gene expression studies by surrogate variable analysis.* PLoS Genetics. 2007.  
  Leek JT. *svaseq: removing batch effects and other unwanted noise from sequencing data.* Nucleic Acids Research. 2014;42(21):e161.  
  https://doi.org/10.1093/nar/gku864

- **ComBat / ComBat-seq**  
  Johnson WE, Li C, Rabinovic A. *Adjusting batch effects in microarray expression data using empirical Bayes methods.* Biostatistics. 2007.  
  Zhang Y, Parmigiani G, Johnson WE. *ComBat-seq: batch effect adjustment for RNA-seq count data.* NAR Genomics and Bioinformatics. 2020.

## QC, trimming & utilities

- **fastp** — Chen S, et al. Bioinformatics. 2018.  
- **FastQC** — Andrews S. Babraham Bioinformatics.  
- **MultiQC** — Ewels P, et al. Bioinformatics. 2016.  
- **ERCC spike-ins** — External RNA Controls Consortium.

## Suggested acknowledgment

> RNA-seq data were processed with a pipeline based on Salmon quantification (Patro et al., 2017), tximport (Soneson et al., 2015), and DESeq2 (Love et al., 2014), with optional RUVSeq (Risso et al., 2014) or svaseq (Leek, 2014) correction for unwanted variation.
