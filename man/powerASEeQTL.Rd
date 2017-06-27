\name{powerASEeQTL}
\alias{powerASEeQTL}
\title{powerASEeQTL}
\description{
	This function implements a power analysis for RNAseq data and then plots it.
}
\usage{ powerASEeQTL(n, mu = 500, n.simu = 200, methods=c("linear", "negativeB", "poisson","TReC","ASEchi2","ASEbinom", "TReCASE"), 
legend = TRUE, color = TRUE, folds = seq(1.5, 1.7, by = 0.2), alpha = 0.001, phi = 1, theta = 0.1, maf = 0.2, 
numcol = 1, title = "", subtitle = "", titlecolor = "black", subtitlecolor = "black", titlesize = 1, subtitlesize = 1, 
legendbox = "o", labelsize = 1, labelcolor = "black", linewidth = 2, tilt = 0, SEED = 1)
}
\arguments{
  \item{n}{
  denotes the sample size.
}
  \item{mu}{
  denotes the total number of reads.
}
  \item{n.simu}{
  denotes the number of simulations.
}
  \item{methods}{
  denotes a vector of the desired methods where "linear" is a
linear regression using an inverse normal transformation for the total
number of read counts. "negativeB" is a negative binomial regression
for the total number of read counts. "poisson" is a poisson regression
for the total number of read counts. "TReC" is an approach by Wei Sun
(2012) A Statistical Framework for eQTL Mapping Using RNA-seq Data,
Biometrics, 2012 Mar;68(1):1-11. "ASEchi2" is a chi-square test for
allele specific expression (ASE). "ASEbinom" is a binomial test for
ASE. "TreCASE" is an approach by Wei Sun that incorporates both eQTL
and ASE analysis.
}
  \item{legend}{
  is TRUE if you want the legend to be displayed for the graph or FALSE
if you don't want the legend displayed on the graph.
}
  \item{color}{
  is TRUE if you want the graph in color or FALSE if you want it in grayscale.
}
  \item{folds}{
  denotes the fold change.
}
  \item{alpha}{
  denotes the significance level.
}
  \item{phi}{
  denotes the dispersion parameter for the negative binomial
distribution for the total number of read counts.
}
  \item{theta}{
  specifies the extent of ASE.
}
  \item{maf}{
  denotes the minor allele frequency.
}
  \item{numcol}{
  specifies number of columns in the legend, default is 1, 5 will be a row.
}
  \item{title}{
  specifies the main title of the graph, default is none.
}
  \item{subtitle}{
  specifies the subtitle of the graph, default is none
}
  \item{titlecolor}{
  denotes what color the title will be
}
  \item{subtitlecolor}{
  denotes what color the subtitle will be
}
  \item{titlesize}{
  denotes what size the title will be, default is 1
}
  \item{subtitlesize}{
  denotes what size the subtitle will be, default is 1
}
  \item{legendbox}{
  legendbox = TRUE will put a box around the legend, default is FALSE
}
  \item{labelsize}{
  denotes what size the axis labels will be, default is 1
}
  \item{labelcolor}{
  denotes what color the axis labels will be
}
  \item{linewidth}{
  denotes the boldness of the lines on the graph
}
  \item{tilt}{
  denotes whether or not the y axis labels are tilted, default is 0, 1 will rotate labels 90 degrees
}
  \item{SEED}{
  sets the seed, default = 1.
}
}

\value{
 A pdf plot of the power analysis of an RNAseq.
}

\author{
Sharon Lutz, Annie Thwing
}

\examples{
}
\keyword{ package }
\keyword{ power }
