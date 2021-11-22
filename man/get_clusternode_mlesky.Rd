\name{get_clusternode_mlesky}
\alias{get_clusternode_mlesky}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_clusternode_mlesky(u = 406318, scanner_env = readRDS("scanner-env-2021-03-03.rds"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{
%%     ~~Describe \code{u} here~~
}
  \item{scanner_env}{
%%     ~~Describe \code{scanner_env} here~~
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
function (u = 406318, scanner_env = readRDS("scanner-env-2021-03-03.rds")) 
{
    library(mlesky)
    library(treedater)
    e1 = as.environment(scanner_env)
    attach(e1)
    mr = 0.00059158
    utre <- keep.tip(tre, descendantSids[[u]])
    sample_times <- sts[utre$tip.label]
    tr <- di2multi(utre, tol = 1e-05)
    tr = unroot(multi2di(tr))
    tr$edge.length <- pmax(1/29000/5, tr$edge.length)
    tr3 <- dater(unroot(tr), sts[tr$tip.label], s = 29000, omega0 = mr, 
        numStartConditions = 0, meanRateLimits = c(mr, mr + 1e-06), 
        ncpu = 6)
    msg <- mlskygrid(tr3, tau = NULL, tau_lower = 0.001, tau_upper = 10, 
        sampleTimes = sts[tr3$tip.label], res = 10, ncpu = 3)
    list(mlesky = msg, timetree = tr3, tree = tr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
