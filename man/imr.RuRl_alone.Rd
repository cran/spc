\name{imr.RuRl_alone}
\alias{imr.RuRl_alone}
\alias{imr.RuRl_alone_s3}
\alias{imr.RuRl_alone_tail}
\alias{imr.Ru_Rlgiven}
\alias{imr.Rl_Rugiven}
\title{Compute control limits of MR charts for normal data}
\description{Computation of control limits of standalone MR charts.}
\usage{imr.RuRl_alone(L0, N=30, qm=30, M0=12, eps=1e-3)

imr.RuRl_alone_s3(L0, N=30, qm=30, M0=12)

imr.RuRl_alone_tail(L0, N=30, qm=30, M0=12)

imr.Ru_Rlgiven(Rl, L0, N=30, qm=30, M0=12)

imr.Rl_Rugiven(Ru, L0, N=30, qm=30, M0=12)}
\arguments{
\item{L0}{pre-defined in-control ARL, that is, determine \code{Ru} and \code{Rl} so that the mean
number of observations until a false alarm is \code{L0}.}
\item{N}{controls the dimension of the linear equation system and consequently the accuracy of the result. See details.}
\item{qm}{number of quadrature nodes (and weights) to determine the definite collocation integrals.}
\item{M0}{mimics Inf --- by setting \code{M0} to some large value (having a standard normal distribution in mind),
the algorithm for IMR charts could be used as well for the standalone MR chart.}
\item{eps}{resolution parameter, which controls the approximation of the ARL slope at the in-control level of
the monitored standard deviation. It ensures the pattern that is called ARL unbiasedness. A small value is recommended.}
\item{Rl}{lower control limit multiple for moving range chart.}
\item{Ru}{upper control limit multiple for moving range chart.}
}
\details{
Crowder (1987a) provided some math to determine the ARL of the so-called individual moving range (IMR) chart,
which consists of the mean X chart and the standard deviation MR chart. 
Making the alarm threshold, \code{M0}, huge (default value here is 12) for the X chart allows us to utilize Crowder's
setup for standalone MR charts. For details about the IMR numerics see \code{imr.arl}.
The three different versions of \code{imr.RuRl_alone} determine limits that form an ARL unbiased design, follow the restriction
\code{Rl} = \code{1/Ru^3} and feature equal probability tails for the MR's half-normal distribution,
respectively in the order given above).
The other two functions are helper routines for \code{imr.RuRl_alone}.
Note that the elegant approach given in Acosta-Mejia/Pignatiello (2000) is only an approximation,
because the MR series is not Markovian.
}
\value{Returns control limit factors (alias multiples).}
\references{
S. V. Crowder (1987a)
Computation of ARL for Combined Individual Measurement and Moving Range Charts,
\emph{Journal of Quality Technology} 19(2), 98-102.
  
S. V. Crowder (1987b)
A Program for the Computation of ARL for Combined Individual Measurement and Moving Range Charts,
\emph{Journal of Quality Technology} 19(2), 103-106.

D. Radson, L. C. Alwan (1995)
Detecting Variance Reductions Using the Moving Range,
\emph{Quality Engineering} 8(1), 165-178.

C. A. Acosta-Mejia, J. J. Pignatiello (2000)
Monitoring process dispersion without subgrouping,
\emph{Journal of Quality Technology} 32(2), 89-102.
}
\author{Sven Knoth}
\seealso{later.}
\examples{
## Radson, Alwan (1995), Table 2 (Monte Carlo based), half-normal, known parameter case
## two-sided MR-alone chart, hence the ARL results has to be decreased by 1
## Here: a large M0=12 (default of the functions above) is deployed to mimic Inf
alpha <- 0.00915
Ru <- sqrt(2) * qnorm(1-alpha/4)
Rl <- sqrt(2) * qnorm(0.5+alpha/4)
M0 <- 12
\dontrun{
ARL0 <- imr.arl(M0, Ru, 0, 1, vsided="two", Rl=Rl)
RRR1995 <- imr.RuRl_alone_tail(ARL0)
RRRs <- imr.RuRl_alone_s3(ARL0)
RRR <- imr.RuRl_alone(ARL0)
results <- rbind(c(Rl, Ru), RRR1995, RRRs, RRR)
results}
}
\keyword{ts}
