# CSCB-2025-Final

## Krishna Bhambhani, Brandon Ly, Varen Talwar

A method to infer CNAs from scRNA-seq data. Relies on the [Python port of InferCNV method](https://infercnvpy.readthedocs.io/en/latest/index.html). scRNA-seq should include expression matrix, cell type annotations, and genomic position information with standardized chromosome names before inferring CNAs.

Notes on usage:
- Ensure that scRNA-seq datasets go through qc and preprocessing steps if needed.
- You may downsample your dataset if desired but should not be necessary.
- Cell type annotations need to be in .obs. If they are missing, the user should classify cells using their method of choice.
- Genomic position information includes chromosome, start, end, and strand (whether forward (1) or reverse (-1) strand). Only chromosome, start, and end are required to be in .obs
- Chromosome names should be standardized to chr1, chr2, ..., chr22, chrX, chrY, chrMT. Any genes located in chromosomes with non-standard names should be removed.