## powerASEeQTL
The powerASEeQTL R package runs power calculations for allele specific expression (ASE) and eQTL analysis using RNA-seq data and creates plots of these results.

## Installation
```
install.packages("devtools") # devtools must be installed first
install.packages("MASS") #this package must be installed first
install.packages("VGAM") #this package must be installed first

devtools::install_github("SharonLutz/powerASEeQTL")
```

## Example
The code below creates a plot of the power calculations for the eQTL and ASE analyses for RNA-seq data with the following parameters:
```
library(powerASEeQTL)

#case 1: methods for eQTL and ASE (faster)
powerASEeQTL(n = 500, mu = 500, n.simu = 500,methods=c("eQTL.linear","eQTL.negBin","eQTL.quasipoisson","ASE"),
folds = seq(1.5, 2, by = 0.25),alpha = 0.001, phi = 1.1, theta = 0.1, maf = 0.2, sim = "sim1")

#case 2: all methods (slower)
powerASEeQTL(n = 500, mu = 500, n.simu = 500,methods=c("eQTL.linear","eQTL.negBin","eQTL.quasipoisson","ASE", "eQTL.ASE"),
folds = seq(1.5, 2, by = 0.25),alpha = 0.001, phi = 1.1, theta = 0.1, maf = 0.2, sim = "sim1")
```

## Output
For this analysis, we get the following plot where the y-axis is power and the x-axis is the fold change.

<img src="https://github.com/SharonLutz/powerASEeQTL/blob/master/powerASEeQTLn500Mu500sim1.png" width="500">

## Reference
Fingerlin TE, Thwing A, **Lutz SM**. (2019) Power Calculations for eQTL Mapping and Allele Specific Expression with RNA-seq Data. (Submitted)

