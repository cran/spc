\name{imr.arl}
\alias{imr.arl}
\alias{imr.Ru_Mgiven}
\alias{imr.Rl_Mgiven}
\alias{imr.MandRu}
\alias{imr.MandRuRl}
\title{Compute ARLs and control limit factors for I-MR combos in case of normal data}
\description{Computation of the (zero-state) Average Run Length (ARL) at given mean \code{mu} and \code{sigma} etc.}
\usage{imr.arl(M, Ru, mu, sigma, vsided="upper", Rl=0, cmode="coll", N=30, qm=30)

imr.Ru_Mgiven(M, L0, N=30, qm=30)

imr.Rl_Mgiven(M, L0, N=30, qm=30)

imr.MandRu(L0, N=30, qm=30)

imr.MandRuRl(L0, N=30, qm=30)}
\arguments{
\item{M}{control limit multiple for mean chart.}
\item{Ru}{upper control limit multiple for moving range chart.}
\item{mu}{actual mean.}
\item{sigma}{actual standard deviation.}
\item{vsided}{switches between the more common "upper" and the less popular "two"(-sided) case of the MR chart.
Setting \code{vsided} to "two" and \code{Ru} sufficiently large (at least \code{2*M}), creates an I-MR chart
with a lower limit only for the MR part.}
\item{Rl}{lower control limit multiple for moving range chart (not needed in the upper case, i.e. if \code{vsided}="upper").}
\item{cmode}{selects the numerical algorithm. The default \code{"coll"} picks the piecewise collocation, which is
the most accurate method. Selecting \code{"Crowder"}, the algorithm from Crowder (1987b) is chosen (re-implemented in R).
Taking a label from \code{"gl"}, \code{"rectangular"}, \code{"trapezoid"}, \code{"simpson"} or \code{"simpson3_8"},
one decides for the quite common Nystroem procedure to numerically solve the considered integral equation.
It is astonishing that Crowder's modified Nystroem design with the trapezoidal quadrature works so well.
However, it is clearly dominated by the piecewise collocation algorithm.}
\item{N}{Controls the dimension of the linear equation system and consequently the accuracy of the result. See details.}
\item{qm}{Number of quadrature nodes (and weights) to determine the collocation definite integrals.}
\item{L0}{pre-defined in-control ARL, that is, determine \code{Ru}, \code{Rl}, or \code{M} and \code{Ru} or
all of them (essentially ending in a lower limit MR chart) so that the mean
number of observations until a false alarm is \code{L0}.}
}
\details{
Crowder (1987a) provided some math to determine the ARL of the so-called individual moving range (IMR) chart.
The given integral equation was approximated by a linear equation system applying trapezoidal quadratures.
Interestingly, Crowder did not recognize the specific behavior of the solution for \code{Ru} >= \code{M} (which is
the more common case), where the resulting function L() is constant in the central part of the
considered domain. In addition, by performing collocation on two (\code{Ru} > \code{M})
or three (\code{Ru} < \code{M}) subintervals piecewise, one obtains highly accurate
ARL numbers. Note that \code{imr.MandRu} yields \code{M} and \code{Ru} for the upper MR trace, whereas
\code{imr.MandRuRl} provides in addition the lower factor \code{Rl} for IMR consisting of two two-sided control charts.
Note that the underlying ARL unbiased condition suppresses the upper limit \code{Ru} in all considered cases so far.
This is not completely surprising, because the mean chart is already quite sensitive for increases in the variance.
The two functions \code{imr.Ru_Mgiven} and \code{imr.Rl_Mgiven} deliver the single upper and lower limit,
respectively, if a one-sided MR design is utilized and the control lmit factor \code{M} of
the mean chart is set already. Note that for \code{Ru} > \code{2*M}, the upper MR limit is
not operative anymore. If it was initially an upper MR chart, then it reduces to a single mean chart.
If it was originally a two-sided MR design, then it becomes a two-sided mean/lower variance chart combo.
Within the latter scheme, the mean chart signals variance increases (a well-known pattern), whereas
the MR subchart delivers only decreasing variance signals. However, these simple Shewhart charts
exhibit in all configurations week variance decreases detection behavior.
Eventually, we should note that the scientific control chart community mainly recommends to
ignore MR charts, see, for example, Vardeman and Jobe (2016), whereas standards (such as ISO), commercial
SPC software and many training manuals provide the IMR scheme with completely wrong upper limits for the MR chart.}
\value{Returns either the ARL or control limit factors (alias multiples).}
\references{
S. V. Crowder (1987a)
Computation of ARL for Combined Individual Measurement and Moving Range Charts,
\emph{Journal of Quality Technology} 19(2), 98-102.
  
S. V. Crowder (1987b)
A Program for the Computation of ARL for Combined Individual Measurement and Moving Range Charts,
\emph{Journal of Quality Technology} 19(2), 103-106.

K. C. B. Roes, R. J. M. M. Does, Y. Schurink,
Shewhart-Type Control Charts for Individual Observations,
\emph{Journal of Quality Technology} 25(3), 188-198.

S. E. Rigdon, E. N. Cruthis, C. W. Champ (1994)
Design Strategies for Individuals and Moving Range Control Charts,
\emph{Journal of Quality Technology} 26(4), 274-287.

D. Radson, L. C. Alwan (1995)
Detecting Variance Reductions Using the Moving Range,
\emph{Quality Engineering} 8(1), 165-178.

S. R. Adke, X. Hong (1997)
A Supplementary Test Based on the Control Chart for Individuals,
\emph{Journal of Quality Technology} 29(1), 16-20.

R. W. Amin, R. A. Ethridge (1998)
A Note on Individual and Moving Range Control Charts,
\emph{Journal of Quality Technology} 30(1), 70-74.

C. A. Acosta-Mejia, J. J. Pignatiello (2000)
Monitoring process dispersion without subgrouping,
\emph{Journal of Quality Technology} 32(2), 89-102.

N. B. Marks, T. C. Krehbiel (2011)
Design And Application Of Individuals And Moving Range Control Charts,
\emph{Journal of Applied Business Research ({JABR})} 25(5), 31-40.

D. Rahardja (2014)
Comparison of Individual and Moving Range Chart Combinations to Individual
Charts in Terms of ARL after Designing for a Common ``All OK'' ARL,
\emph{Journal of Modern Applied Statistical Methods} 13(2), 364-378.

S. B. Vardeman, J. M. Jobe (2016)
\emph{Statistical Methods for Quality Assurance},
Springer, 2nd edition.
}
\author{Sven Knoth}
\seealso{later.}
\examples{
## Crowder (1987b), Output Listing 1, trapezoidal quadrature (less accurate)

M <- 2
Ru <- 3
mu <- seq(0, 2, by=0.25)
LL <- LL2 <- rep(NA, length(mu))
for ( i in 1:length(mu) ) {
  LL[i] <- round( imr.arl(M, Ru, mu[i], 1), digits=4)
  LL2[i] <- round( imr.arl(M, Ru, mu[i], 1, cmode="Crowder", N=80), digits=4)
}
LL1987b <- c(18.2164, 16.3541, 12.4282, 8.7559, 6.1071, 4.3582, 3.2260, 2.4878, 1.9989)
print( data.frame(mu, LL2, LL1987b, LL) )

## Crowder (1987a), Table 1, trapezoidal quadrature (less accurate)

M <- 4
Ru <- 3
mu <- seq(0, 2, by=0.25)
LL <- rep(NA, length(mu))
for ( i in 1:length(mu) ) LL[i] <- round( imr.arl(M, Ru, mu[i], 1), digits=4)
LL1987a <- c(34.44, 34.28, 34.07, 33.81, 33.45, 32.82, 31.50, 28.85, 24.49)
print( data.frame(mu, LL1987a, LL) )

## Rigdon, Cruthis, Champ (1994), Table 1, Monte Carlo based

M <- 2.992
Ru <- 4.139
icARL <- imr.arl(M, Ru, 0, 1)
icARL1994 <- 200
print( data.frame(icARL1994, icARL) )

M <- 3.268
Ru <- 4.556
icARL <- imr.arl(M, Ru, 0, 1)
icARL1994 <- 500
print( data.frame(icARL1994, icARL) )

## ..., Table 2, Monte Carlo based

M <- 2.992
Ru <- 4.139
tau <- c(seq(1, 1.3, by=0.05), seq(1.4, 2, by=0.1))
LL <- rep(NA, length(tau))
for ( i in 1:length(tau) ) LL[i] <- round( imr.arl(M, Ru, 0, tau[i]), digits=2)
LL1994 <- c(200.54, 132.25, 90.84, 65.66, 49.35, 38.92, 31.11, 21.35, 15.47,
12.04, 9.81, 8.21, 7.03, 6.14)
print( data.frame(tau, LL1994, LL) )

## Radson, Alwan (1995), Table 2 (Monte Carlo based), half-normal, known parameter case
## two-sided (!) MR-alone (!) chart, hence the ARL results has to be decreased by 1
## Here: a large M (=12) is deployed to mimic Inf
alpha <- 0.00915
Ru <- sqrt(2) * qnorm(1-alpha/4)
Rl <- sqrt(2) * qnorm(0.5+alpha/4)
k <- 1.5 - (0:7)/10
LL <- rep(NA, length(k))
for ( i in 1:length(k) )
  LL[i] <- round( imr.arl(12, Ru, 0, k[i], vsided="two", Rl=Rl), digits=2) - 1
RA1995 <- c(18.61, 24.51, 34.21, 49.74, 75.08, 113.14, 150.15, 164.54)
print( data.frame(k, RA1995, LL) )

## Amin, Ethridge (1998), Table 2, column sigma/sigma_0 = 1.00

M <- 3.27
Ru <- 4.56
#M <- 3.268
#Ru <- 4.556
mu <- seq(0, 2, by=0.25)
LL <- rep(NA, length(mu))
for ( i in 1:length(mu) ) LL[i] <- round( imr.arl(M, Ru, mu[i], 1), digits=1)
LL1998 <- c(505.3, 427.6, 276.7, 156.2, 85.0, 46.9, 26.9, 16.1, 10.1)
print( data.frame(mu, LL1998, LL) )

## ..., column sigma/sigma_0 = 1.05

for ( i in 1:length(mu) ) LL[i] <- round( imr.arl(M, Ru, mu[i], 1.05), digits=1)
LL1998 <- c(296.8, 251.6, 169.6, 101.6, 58.9, 34.5, 20.9, 13.2, 8.7)
print( data.frame(mu, LL1998, LL) )

## Acosta-Mejia, Pignatiello (2000), Table 2
## AMP utilized Markov chain approximation
## However, the MR series is not Markovian!
## MR-alone (!) chart, hence the ARL results has to be decreased by 1
## Here: a large M (=8) is deployed to mimic Inf
Ru <- 3.93
sigma <- c(1, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.75)
LL <- rep(NA, length(sigma))
for ( i in 1:length(sigma) ) LL[i] <- round( imr.arl(8, Ru, 0, sigma[i], N=30), digits=1) - 1
AMP2000 <- c(201.0, 136.8, 97.9, 73.0, 56.3, 36.4, 25.6, 19.1, 11.0)
print( data.frame(sigma, AMP2000, LL) )

## Mark, Krehbiel (2011), Table 2, deployment of Crowder (1987b), nominal ic ARL 500

M <- c(3.09, 3.20, 3.30, 3.50, 4.00)
Ru <- c(6.00, 4.67, 4.53, 4.42, 4.36)
LL0 <- rep(NA, length(M))
for ( i in 1:length(M) ) LL0[i] <- round( imr.arl(M[i], Ru[i], 0, 1), digits=1)
print( data.frame(M, Ru, LL0) )
}
\keyword{ts}
