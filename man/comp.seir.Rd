\name{comp.seir}
%\alias{comp.seir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solves the SEIR compartmental model
}
\description{
Gives both deterministic and stochastic solutions for the classic SEIR
epidemic model.
}
\usage{
comp.seir(
  De = 1.48, 
  Tg = 2.6, 
  R0 = 1.8, 
  N = 6800000, 
  I0 = 10, 
  dt = 1, 
  R1 = 1, 
  t1 = 999, 
  R2 = 1, 
  t2 = 9999, 
  noTimeSteps = 10,
  noReals = 1, 
  deterministic = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{De}{
The average duration of the E stage, the period during which the host is infected but not yet infectious.
}
  \item{Tg}{
The generation time of the model. The duration of infectiousness Di is inferred from the generation time and the duration of the E stage De, according to Tg = De + Di.
}
  \item{R0}{
The basic reproductive number before time t1. The value of beta is inferred from the value of R0 according to R0 = beta Di.
}
  \item{N}{
The size of the population.
}
  \item{I0}{
The number initially infectious.
}
  \item{dt}{
The duration of time assumed to occur for each timestep.
}
  \item{R1}{
The basic reproductive number after time t1 but before time t2. The value of beta is inferred from the value of R1 according to R1 = beta Di. 
}
  \item{t1}{
The time step at which the reproduction number changes from R0 to R1.
}
  \item{R2}{
The basic reproductive number after time t2. The value of beta is inferred from the value of R2 according to R2 = beta Di. 
}
  \item{t2}{
The time step ay which the reproduction number changes from R1 to R2.
}
  \item{noTimeSteps}{
The total number of time steps that the solution is obtained for.
}
  \item{noReals}{
Number of realizations generated for stochastic solutions. Ignored for deterministic solutions.
}
  \item{deterministic}{
If FALSE, the model solution is deterministic, using a positive definite Euler-type integration. If FALSE, the model is solved using a fixed time step stochastic algorithm. See details.
}
}
\details{
The model in terms of state variables and the haard that events occur.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{inc}{Description here}
}
\references{
This is a very commonly used model, but is usually specified as a deterministic model. See the supporting info of the paper below for a nice description of the stochastic formulation.

Transmission dynamics of the etiological agent of SARS in Hong Kong: impact of public health interventions. Riley S, Fraser C, Donnelly CA, Ghani AC, Abu-Raddad LJ, Hedley AJ, Leung GM, Ho LM, Lam TH, Thach TQ, Chau P, Chan KP, Lo SV, Leung PY, Tsang T, Ho W, Lee KH, Lau EM, Ferguson NM, Anderson RM. Science. 2003 Jun 20;300(5627):1961-6. 

}
\author{
Steven Riley.
}
%%\note{%%~~further notes~~}

%%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}

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
