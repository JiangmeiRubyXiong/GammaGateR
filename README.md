# GammaGateR <img src="logoGammaGateR.png" align="right" height = "150" />

Install cfGMM prior to installing GammaGateR:
```
devtools::install_github("JiangmeiRubyXiong/cfGMM")
```

Install GammaGateR: (vignette build might take a while, so I recommend set it to be false for initial download, and read the html in the vignette fold directly):
```
devtools::install_github("JiangmeiRubyXiong/GammaGateR", build_vignettes = FALSE)
```

To run the example in the vignette, Bioconductor installation is needed to load the publicly available dataset. The `VectraPolarisData` package on Bioconductor contains multiplex immunohistochemistry (mIHC) in a sample of 153 patients with small cell lung cancer(Johnson et al. 2021; Wrobel and Ghosh 2022). To download the dataset, follow the instruction at https://bioconductor.org/packages/release/data/experiment/html/VectraPolarisData.html. After installation, run the following to store data as the cell object. The data is are cell-level, meaning that each row of the data represents a single cell.
