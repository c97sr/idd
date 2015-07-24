\name{comp.seir}
\alias{comp.seir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solves the SEIR compartmental model
}
\description{
Gives both deterministic and stochastic solutions for the calssic SEIR
epidemic model.
}
\usage{
comp.seir(De = 1.48, Tg = 2.6, R0 = 1.8, N = 6800000, I0 = 10, dt = 1, R1 = 1, t1 = 999, R2 = 1, t2 = 9999, noTimeSteps = 10, deterministic = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{De}{
%%     ~~Describe \code{De} here~~
}
  \item{Tg}{
%%     ~~Describe \code{Tg} here~~
}
  \item{R0}{
%%     ~~Describe \code{R0} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{I0}{
%%     ~~Describe \code{I0} here~~
}
  \item{dt}{
%%     ~~Describe \code{dt} here~~
}
  \item{R1}{
%%     ~~Describe \code{R1} here~~
}
  \item{t1}{
%%     ~~Describe \code{t1} here~~
}
  \item{R2}{
%%     ~~Describe \code{R2} here~~
}
  \item{t2}{
%%     ~~Describe \code{t2} here~~
}
  \item{noTimeSteps}{
%%     ~~Describe \code{noTimeSteps} here~~
}
  \item{deterministic}{
%%     ~~Describe \code{deterministic} here~~
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
function (De = 1.48, Tg = 2.6, R0 = 1.8, N = 6800000, I0 = 10,
    dt = 1, R1 = 1, t1 = 999, R2 = 1, t2 = 9999, noTimeSteps = 10,
    deterministic = TRUE)
{
    rtn_inf_inc <- vector(mode = "numeric", length = noTimeSteps +
        1)
    Di <- Tg - De
    state_S <- N - I0
    state_E <- 0
    state_I <- I0
    state_R <- 0
    rtn_inf_inc[1] <- I0
    for (i in 1:noTimeSteps) {
        if (i < t1 * dt) {
            beta = R0/Di
        }
        else if (i < t2 * dt) {
            beta = R1/Di
        }
        else {
            beta = R2/Di
        }
        hazInf <- beta * state_I/N
        hazBecInf <- 1/De
        hazRec <- 1/Di
        pInf <- 1 - exp(-dt * hazInf)
        pBecInf <- 1 - exp(-dt * hazBecInf)
        pRec <- 1 - exp(-dt * hazRec)
        if (deterministic) {
            nInf <- state_S * pInf
            nBecInf <- state_E * pBecInf
            nRec <- state_I * pRec
        }
        else {
            nInf <- rbinom(1, state_S, pInf)
            nBecInf <- rbinom(1, state_E, pBecInf)
            nRec <- rbinom(1, state_I, pRec)
        }
        state_S <- state_S - nInf
        state_E <- state_E + nInf - nBecInf
        state_I <- state_I + nBecInf - nRec
        state_R <- state_R - nRec
        rtn_inf_inc[i + 1] <- nInf
    }
    list(inf_inc = rtn_inf_inc, time = (0:noTimeSteps) * dt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
