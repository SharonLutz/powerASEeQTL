\name{logLTReC}
\alias{logLTReC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
logLTReC
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logLTReC(bxj, y, x, mu, b0, phi, fam)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bxj}{
%%     ~~Describe \code{bxj} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{b0}{
%%     ~~Describe \code{b0} here~~
}
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{fam}{
%%     ~~Describe \code{fam} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (bxj, y, x, mu, b0, phi, fam) 
{
    mu1 = mu
    ww2 = which(x == 2)
    ww1 = which(x == 1)
    mu1[ww2] = mu[ww2] * exp(bxj - b0)
    mu1[ww1] = mu[ww1] * (1 + exp(bxj))/(1 + exp(b0))
    if (fam == "negbin") {
        logL = loglikNB(phi, mu1, y)
    }
    else {
        logL = 0
        for (i in 1:length(y)) {
            logL = logL + dpois(y[i], mu1[i], log = TRUE)
        }
    }
    logL
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
