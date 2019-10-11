\name{pois.cusum.arl}
\alias{pois.cusum.arl}
\title{Compute ARLs of Poisson CUSUM control charts}
\description{Computation of the (zero-state) Average Run Length (ARL) at given mean \code{mu}.}
\usage{pois.cusum.arl(mu, km, hm, m, i0=0, sided="upper", rando=FALSE,
gamma=0, km2=0, hm2=0, m2=0, i02=0, gamma2=0)}
\arguments{
\item{mu}{actual mean.}
\item{km}{enumerator of rational approximation of reference value \code{k}.}
\item{hm}{enumerator of rational approximation of reference value \code{h}.}
\item{m}{denominator of rational approximation of reference value.}
\item{i0}{head start value as integer multiple of \code{1/m}; should be an element of \code{0:hm}.}
\item{sided}{distinguishes between different one- and two-sided CUSUM control chart by choosing
\code{"upper"}, \code{"lower"} and \code{"two"}, respectively.}
\item{rando}{Switch for activating randomization in order to allow continuous ARL control.}
\item{gamma}{Randomization probability. If the CUSUM statistic is equal to the threshold \code{h},
an control chart alarm is triggered with probability \code{gamma}.}
\item{km2,hm2,m2,i02,gamma2}{corresponding values of the second CUSUM chart (to building a two-sided CUSUM scheme).}
}
\details{
The monitored data follow a Poisson distribution with \code{mu}.
The ARL values of the resulting EWMA control chart are determined via Markov chain calculations.
We follow the algorithm given in Lucas (1985) expanded with some arithmetic 'tricks' (e.g., by deploying
Toeplitz matrix algebra).
}
\value{Returns a single value which resembles the ARL.}
\references{
J. M. Lucas (1985)
Counted data CUSUM's,
\emph{Technometrics} 27(2), 129-144.

C. H. White and J. B. Keats (1996)
ARLs and Higher-Order Run-Length Moments for the Poisson CUSUM,
\emph{Journal of Quality Technology} 28(3), 363-369.

C. H. White, J. B. Keats and J. Stanley (1997)
Poisson CUSUM versus c chart for defect data,
\emph{Quality Engineering} 9(4), 673-679.

G. Rossi and L. Lampugnani and M. Marchi (1999),
An approximate CUSUM procedure for surveillance of health events,
\emph{Statistics in Medicine} 18(16), 2111-2122.

S. W. Han, K.-L. Tsui, B. Ariyajunya, and S. B. Kim (2010),
A comparison of CUSUM, EWMA, and temporal scan statistics for detection of increases in poisson rates,
\emph{Quality and Reliability Engineering International} 26(3), 279-289.

M. B. Perry and J. J. Pignatiello Jr. (2011)
Estimating the time of step change with Poisson CUSUM and EWMA control charts,
\emph{International Journal of Production Research} 49(10), 2857-2871.

}
\author{Sven Knoth}
\seealso{later.}
\examples{
## Lucas 1985, upper chart (Tables 2 and 3)
k   <- .25
h   <- 10
m   <- 4
km  <- m * k
hm  <- m * h
mu0 <- 1 * k
ARL <- pois.cusum.arl(mu0, km, hm-1, m)
# Lucas reported 438 (in Table 2, first block, row 10.0 .25 .0 ..., column 1.0
# Recall that Lucas and other trigger an alarm, if the CUSUM statistic is greater than
# or equal to the alarm threshold h
print(ARL)

ARL <- pois.cusum.arl(mu0, km, hm-1, m, i0=round((hm-1)/2))
# Lucas reported 333 (in Table 3, first block, row 10.0 .25 .0 ..., column 1.0
print(ARL)

## Lucas 1985, lower chart (Tables 4 and 5)
ARL <- pois.cusum.arl(mu0, km, hm-1, m, sided="lower")
# Lucas reported 437 (in Table 4, first block, row 10.0 .25 .0 ..., column 1.0
print(ARL)

ARL <- pois.cusum.arl(mu0, km, hm-1, m, i0=round((hm-1)/2), sided="lower")
# Lucas reported 318 (in Table 5, first block, row 10.0 .25 .0 ..., column 1.0
print(ARL)
}
\keyword{ts}
