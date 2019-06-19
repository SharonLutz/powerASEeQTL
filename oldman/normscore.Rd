\name{normscore}
\alias{normscore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{normscore
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normscore(vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
%%     ~~Describe \code{vec} here~~
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
function (vec) 
{
    len = length(na.omit(vec)) + 1
    rank = rank(na.omit(vec))
    ties = (rank - floor(rank)) > 0
    new.vec = vec[!is.na(vec)]
    new.vec[!ties] = qnorm(rank[!ties]/len)
    new.vec[ties] = 0.5 * (qnorm((rank[ties] + 0.5)/len) + qnorm((rank[ties] - 
        0.5)/len))
    vec[!is.na(vec)] = new.vec
    vec
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
