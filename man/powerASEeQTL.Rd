\name{powerASEeQTL}
\alias{powerASEeQTL}
\title{powerASEeQTL}
\description{
	This function implements a power analysis for RNAseq data and then plots it.
}
\usage{ powerASEeQTL(n, mu = 500, n.simu = 200, methods = c("Linear Model", "Negative Binomial", "TReC", "TReCASE", "ASEchi2","ASEbinom","Poisson"), 
legend = TRUE, color = TRUE, folds = seq(1.5, 1.7, by = 0.2), alpha = 0.001, phi = 1, theta = 0.1, maf = 0.2, 
numcol = 1, title = "", subtitle = "", titlecolor = "black", subtitlecolor = "black", titlesize = 1, subtitlesize = 1, 
legendbox = "o", labelsize = 1, labelcolor = "black", linewidth = 2, tilt = 0)
}
\arguments{
  \item{n}{
n is the sample size
}
  \item{mu}{
mu is the total number of reads
}
  \item{n.simu}{
n.simu is the number of simulations
}
  \item{methods}{
methods is a vector of the desired methods
}
  \item{legend}{
legend is TRUE if you want the legend to be displayed or FALSE if you want it hidden
}
  \item{color}{
color is TRUE if you want the graph in color or FALSE if you want it in grayscale
}
  \item{folds}{

}
  \item{alpha}{
alpha is the significance level
}
  \item{phi}{

}
  \item{theta}{

}
  \item{maf}{
maf is the minor allele frequency
}
  \item{numcol}{
numcol specifies number of columns in the legend, default is 1, 5 will be a row
}
  \item{title}{
title is a string which will set the title, default is none
}
  \item{subtitle}{
title is a string which will set the subtitle, default is none
}
  \item{titlecolor}{
titlecolor indicates what color the title will be
}
  \item{subtitlecolor}{
subtitlecolor indicates what color the subtitle will be
}
  \item{titlesize}{
titlecolor indicates what size the title will be, default is 1
}
  \item{subtitlesize}{
subtitlesize indicates what size the subtitle will be, default is 1
}
  \item{legendbox}{
legendbox = TRUE will put a box around the legend, default is FALSE
}
  \item{labelsize}{
labelsize indicates what size the axis labels will be, default is 1
}
  \item{labelcolor}{
labelcolor indicates what color the axis labels will be
}
  \item{linewidth}{
linewidth indicates the boldness of the lines on the graph
}
  \item{tilt}{
tilt indicates whether or not the y axis labels are tilted, default is 0, 1 will rotate labels 90 degrees
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
