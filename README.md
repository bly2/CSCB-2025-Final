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

Further testing of methods to be done on 3 selected PSC datasets.

### How to install and use

To install, please download the cscb_methods folder from our repository and store it in your working directory alongside your datasets. Import our methods:

``` py
from cscb_methods import *
```

Please kindly refer to ```benchmark_PBMC.ipynb``` for the workflow. In summary:

- Ensure chromosomes have standardized names and contain position data from Biomart.
- Perform QC steps, normalize, log transform, and/or downsample as necessary.
- 3 state HMM function will return a new AnnData with results in ```adata.obs['hmm_cnv]```. It needs to contain a column for diploid/aneuploid predictions in .obs and an input for cell type to focus on.

inferCNV of the Trinity CTAT Project.  https://github.com/broadinstitute/inferCNV

inferCNVpy. https://github.com/icbi-lab/infercnvpy

CopyKAT of Navin Lab at MD Anderson. https://github.com/navinlabcode/copykat
