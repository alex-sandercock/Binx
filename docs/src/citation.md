# Citation

If Binx is useful in your research, please cite both the original methods and the Binx implementation.

## Citing Binx

```
Sandercock, A.M. (2025). Binx: A Rust-based CLI tool for polyploid and diploid 
genomic analysis. GitHub repository: https://github.com/alex-sandercock/Binx
```

BibTeX:
```bibtex
@software{sandercock2025binx,
  author = {Sandercock, Alexander M.},
  title = {Binx: A Rust-based CLI tool for polyploid and diploid genomic analysis},
  year = {2025},
  url = {https://github.com/alex-sandercock/Binx}
}
```

> **Note:** A formal publication with DOI is planned. This section will be updated when available.

## Citing the Methods

### rrBLUP (Mixed Model)

If using the GWAS or kinship functionality:

```
Endelman, J.B. (2011). Ridge regression and other kernels for genomic selection 
with R package rrBLUP. The Plant Genome 4:250-255.
```

```bibtex
@article{endelman2011rrblup,
  author = {Endelman, Jeffrey B.},
  title = {Ridge regression and other kernels for genomic selection with R package rrBLUP},
  journal = {The Plant Genome},
  volume = {4},
  pages = {250--255},
  year = {2011},
  doi = {10.3835/plantgenome2011.08.0024}
}
```

### GWASpoly (Polyploid GWAS)

If using polyploid genetic models:

```
Rosyara, U.R., De Jong, W.S., Douches, D.S., & Endelman, J.B. (2016). Software 
for genome-wide association studies in autopolyploids and its application to 
potato. The Plant Genome 9(2).
```

```bibtex
@article{rosyara2016gwaspoly,
  author = {Rosyara, Umesh R. and De Jong, Walter S. and Douches, David S. and Endelman, Jeffrey B.},
  title = {Software for genome-wide association studies in autopolyploids and its application to potato},
  journal = {The Plant Genome},
  volume = {9},
  number = {2},
  year = {2016},
  doi = {10.3835/plantgenome2015.08.0073}
}
```

### Updog (Dosage Estimation)

If using the dosage estimation functionality:

```
Gerard, D., Ferrão, L.F.V., Garcia, A.A.F., & Stephens, M. (2018). Genotyping 
polyploids from messy sequencing data. Genetics 210(3):789-807.
```

```bibtex
@article{gerard2018updog,
  author = {Gerard, David and Ferrão, Luis Felipe Ventorim and Garcia, Antonio Augusto Franco and Stephens, Matthew},
  title = {Genotyping polyploids from messy sequencing data},
  journal = {Genetics},
  volume = {210},
  number = {3},
  pages = {789--807},
  year = {2018},
  doi = {10.1534/genetics.118.301468}
}
```

## Example Acknowledgment

In your methods section:

> GWAS was performed using Binx v0.1.0 (Sandercock, 2025), which implements 
> the GWASpoly framework (Rosyara et al., 2016) and rrBLUP mixed model 
> (Endelman, 2011) in Rust. Multiple genetic models were tested including 
> additive and general models appropriate for autotetraploid inheritance.
