\name{logBB}
\alias{logBB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
logBB
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logBB(ni, ni0, pi, theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ni}{
%%     ~~Describe \code{ni} here~~
}
  \item{ni0}{
%%     ~~Describe \code{ni0} here~~
}
  \item{pi}{
%%     ~~Describe \code{pi} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
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
function (ni, ni0, pi, theta) 
{
    tmp0 = lchoose(ni, ni0)
    if (ni0 > 0) {
        seq1 = seq(0, ni0 - 1, by = 1)
        for (k in seq1) {
            tmp0 = tmp0 + log(pi + k * theta)
        }
    }
    if (ni0 < ni) {
        seq2 = seq(0, ni - ni0 - 1, by = 1)
        for (k in seq2) {
            tmp0 = tmp0 + log(1 - pi + k * theta)
        }
    }
    seq3 = 1:(ni - 1)
    for (k in seq3) {
        tmp0 = tmp0 - log(1 + k * theta)
    }
    tmp0
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
