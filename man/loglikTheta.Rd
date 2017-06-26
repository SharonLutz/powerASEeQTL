\name{loglikTheta}
\alias{loglikTheta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
loglikTheta
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loglikTheta(theta, pi, nA, nTotal, zeta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{pi}{
%%     ~~Describe \code{pi} here~~
}
  \item{nA}{
%%     ~~Describe \code{nA} here~~
}
  \item{nTotal}{
%%     ~~Describe \code{nTotal} here~~
}
  \item{zeta}{
%%     ~~Describe \code{zeta} here~~
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
function (theta, pi, nA, nTotal, zeta) 
{
    sumL = 0
    for (i in 1:length(nA)) {
        ni = nTotal[i]
        ni0 = nA[i]
        if (zeta[i]) {
            sumL = sumL + logBB(ni, ni0, pi, theta)
        }
        else {
            sumL = sumL + logBB(ni, ni0, 0.5, theta)
        }
    }
    sumL
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
