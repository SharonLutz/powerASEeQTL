\name{loglikNB}
\alias{loglikNB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
loglikNB
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loglikNB(phi, mu, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
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
function (phi, mu, y) 
{
    logL = 0
    vphi = 1/phi
    lik0 = vphi * log(vphi) - lgamma(vphi)
    for (i in 1:length(y)) {
        yi = y[i]
        mui = mu[i]
        if (yi == 0) {
            logL = logL + vphi * log(vphi) - vphi * log(vphi + 
                mui)
        }
        else {
            logL = logL + lgamma(yi + vphi) - lgamma(yi + 1) + 
                yi * log(mui)
            logL = logL - (vphi + yi) * log(vphi + mui) + lik0
        }
    }
    logL
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
