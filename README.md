# Epigenomic profiling of mouse nucleus accumbens at single-cell resolution
This repository contains companion code for Bhatia &amp; Yang et al. (2023) and includes code for gene regulatory network (GRN) discovery, and analyses of single-nucleus chromatin accessibility and gene expression data for the mouse nucleus accumbens.

Raw FASTQ files, processed Seurat object and associated meta data with cell type annotations can be found at [GSE232774](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232774)

## Software requirements
- [Signac](https://stuartlab.org/signac/)
- [ChIPseeker](http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html)
- [SCENIC+](https://scenicplus.readthedocs.io/en/latest/install.html)
  - and associated dependencies
- [pandas](https://pandas.pydata.org/docs/getting_started/install.html)

## Required data
- snATAC-seq fragments file (fragments.tsv.gz)
  - and associated metadata (nucleus_accumbens_metadata.csv)
- scRNA-seq data (.h5ad file)


## References
González-Blas, C. B., Minnoye, L., Papasokrati, D., Aibar, S., Hulselmans, G., Christiaens, V., Davie, K., Wouters, J., & Aerts, S. (2019). cisTopic: Cis-regulatory topic modeling on single-cell ATAC-seq data. Nature Methods, 16(5), Article 5. https://doi.org/10.1038/s41592-019-0367-1

González-Blas, C. B., Winter, S. D., Hulselmans, G., Hecker, N., Matetovici, I., Christiaens, V., Poovathingal, S., Wouters, J., Aibar, S., & Aerts, S. (2022). SCENIC+: Single-cell multiomic inference of enhancers and gene regulatory networks (p. 2022.08.19.504505). bioRxiv. https://doi.org/10.1101/2022.08.19.504505

Stuart T, Srivastava A, Madad S, Lareau CA, Satija R. Single-cell chromatin state analysis with Signac. Nat Methods. 2021 Nov;18(11):1333-1341. doi: 10.1038/s41592-021-01282-5. PMID: 34725479; PMCID: PMC9255697.

Yu G, Wang LG, He QY. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. Bioinformatics. 2015 Jul 15;31(14):2382-3. doi: 10.1093/bioinformatics/btv145.



