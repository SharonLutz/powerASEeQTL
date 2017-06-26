\name{gradLogH0}
\alias{gradLogH0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gradLogH0
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gradLogH0(par, nA, nTotal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
%%     ~~Describe \code{par} here~~
}
  \item{nA}{
%%     ~~Describe \code{nA} here~~
}
  \item{nTotal}{
%%     ~~Describe \code{nTotal} here~~
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
function (par, nA, nTotal) 
{
    theta = par[1]
    pi0 = 0.5
    gradTh = 0
    for (i in 1:length(nA)) {
        ni = nTotal[i]
        ni0 = nA[i]
        gri = gradBB(ni, ni0, pi0, theta)
        gradTh = gradTh + gri[1]
    }
    gradTh
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
