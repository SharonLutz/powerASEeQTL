\name{loglikJoin}
\alias{loglikJoin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{loglikJoin
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loglikJoin(bxj, y, x, nA, nTotal, zeta, mu, b0, phi, theta)
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
  \item{nA}{
%%     ~~Describe \code{nA} here~~
}
  \item{nTotal}{
%%     ~~Describe \code{nTotal} here~~
}
  \item{zeta}{
%%     ~~Describe \code{zeta} here~~
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
function (bxj, y, x, nA, nTotal, zeta, mu, b0, phi, theta) 
{
    mu1 = mu
    ww2 = which(x == 2)
    ww1 = which(x == 1)
    mu1[ww2] = exp(log(mu[ww2]) + bxj - b0)
    mu1[ww1] = exp(log(mu[ww1]) + log(1 + exp(bxj)) - log(1 + 
        exp(b0)))
    pi1 = exp(bxj)/(1 + exp(bxj))
    par = c(theta, pi1)
    logTReC = loglikNB(phi, mu1, y)
    logASE = logH1(par, nA, nTotal, zeta)
    logTReC + logASE
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
