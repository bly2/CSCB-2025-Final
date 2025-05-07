# CSCB-2025-Final

## Krishna Bhambhani, Brandon Ly, Varen Talwar

A method to infer CNAs from scRNA-seq data. Relies on the [Python port of InferCNV method](https://infercnvpy.readthedocs.io/en/latest/index.html). scRNA-seq should include expression matrix, cell type annotations, and genomic position information with standardized chromosome names before inferring CNAs.

Summary of methods employed:
- Gaussian Mixture Model adopted from CopyKAT to identify diploid and aneuploid cells before CNV detection
- Hidden Markov Model adopted from InferCNV to assign CNV states deletion, neutral, and amplification to aneuploid cells compared to reference diploid cells
- To be compared to standard workflows for InferCNVpy and maybe CopyKAT (which may would require heavy downsampling)

Initial benchmarking to be done on PBMC dataset with simulated CNVs:
- Copy number gain for CD4 T cells on ChrX
- Partial/heterozygous loss for CD14 monocytes on Chr22
- Total/homozygous loss for CD14 monocytes on Chr6

Further testing of methods to be done on 3 selected PSC datasets:

inferCNV of the Trinity CTAT Project.  https://github.com/broadinstitute/inferCNV
inferCNVpy. https://github.com/icbi-lab/infercnvpy
CopyKAT of Navin Lab at MD Anderson. https://github.com/navinlabcode/copykat
