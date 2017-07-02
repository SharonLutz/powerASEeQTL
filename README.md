## powerASEeQTL
The powerASEeQTL R package runs power calculations for allele specific expression (ASE) and eQTL analysis using RNA-seq data and creates plots of these results.
#### Installation
```
install.packages("devtools") # devtools must be installed first
install.packages("MASS") #this package must be installed first
install.packages("VGAM") #this package must be installed first

devtools::install_github("SharonLutz/powerASEeQTL")
```
#### Example
The code below creates a plot of the power calculations for the eQTL and ASE anlayses for RNA-seq data with the following parameters:
```
library(powerASEeQTL)
?powerASEeQTL # For details on this function and how to customize its graphical output

powerASEeQTL(n = 100, mu = 500, n.simu = 200,folds = seq(1.5, 2.5, by = 0.5), 
alpha = 0.001, phi = 1, theta = 0.1, maf = 0.2)

```

#### Output
For this analysis we get the following matrix for the power calculations and corresponding plot:
```
  linear negativeB poisson  TReC ASEchi2 ASEbinom TReCASE
1   0.01      0.01   0.030 0.010   0.015    0.035   0.050
2   0.13      0.17   0.225 0.180   0.185    0.210   0.445
3   0.26      0.35   0.445 0.375   0.375    0.445   0.805
```
<img src="https://github.com/SharonLutz/powerASEeQTL/blob/master/powerRNAseqN100Mu500.pdf" width="200">


#### Reference
Lutz SM, Thwing A, Fingerlin TE. (2017) Power Calculations for eQTL Mapping and Allele Specific Expression with RNA-seq Data. (In Preparation)
