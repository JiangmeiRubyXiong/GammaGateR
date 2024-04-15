# GammaGateR <img src="logoGammaGateR.png" align="right" height = "150" />

Install cfGMM prior to installing GammaGateR:
```
devtools::install_github("JiangmeiRubyXiong/cfGMM")
```

Install GammaGateR: (The vignette build takes a while, so I recommend setting it to be false for the initial download. The compiled vignette can be viewed [here](https://statimagcoll.github.io/GammaGateR): 
```
devtools::install_github("JiangmeiRubyXiong/GammaGateR",build_vignettes = FALSE)
```

To run the example in the vignette, Bioconductor installation is needed to load the publicly available dataset. The VectraPolarisData package on Bioconductor contains multiplex immunohistochemistry (mIHC) in a sample of 153 patients with small cell lung cancer (Johnson et al. 2021; Wrobel and Ghosh 2022). To download the dataset, follow the instructions [here](https://bioconductor.org/packages/release/data/experiment/html/VectraPolarisData.html).

