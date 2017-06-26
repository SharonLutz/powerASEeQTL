\name{grad.theta}
\alias{grad.theta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
grad.theta
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
grad.theta(theta, nA, nTotal, zeta, pi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
%%     ~~Describe \code{theta} here~~
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
  \item{pi}{
%%     ~~Describe \code{pi} here~~
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
function (theta, nA, nTotal, zeta, pi) 
{
    grad = 0
    for (i in 1:length(nA)) {
        if (nA[i] > 0) {
            ks = 0:(nA[i] - 1)
            if (zeta[i] > 0) {
                grad = grad + sum(ks/(pi + ks * theta))
            }
            else {
                grad = grad + sum(ks/(0.5 + ks * theta))
            }
        }
        if (nA[i] < nTotal[i]) {
            ks = 0:(nTotal[i] - nA[i] - 1)
            if (zeta[i] > 0) {
                grad = grad + sum(ks/(1 - pi + ks * theta))
            }
            else {
                grad = grad + sum(ks/(0.5 + ks * theta))
            }
        }
        ks = 1:(nTotal[i] - 1)
        grad = grad - sum(ks/(1 + ks * theta))
    }
    grad
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
