\name{ind.tau}
\alias{ind.tau}
\title{
Simulates an individual based case tree model of infectious disease transmission
}
\description{
The model produces times generations of infection for a given offspring distribution and serial interval distribution.
}
\usage{
ind.tau(Tg = 5, R0 = 2.6, Nmax = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tg}{
Mean of the serial interval distribution
}
  \item{R0}{
The basic reproductive number
}
  \item{Nmax}{
The maximium number of cases in the epidemic
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
function (Tg = 5, R0 = 2.6, Nmax = 30) 
{
    rtn_times <- numeric(mode = "numeric", length = Nmax)
    rtn_infectors <- numeric(mode = "numeric", length = Nmax)
    rtn_times[1] <- 0
    rtn_infectors[1] <- 0
    completed_infector <- 0
    max_infector <- 1
    while (completed_infector < max_infector) {
    }
    list(t = rtn_times, i = rtn_infectors)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
