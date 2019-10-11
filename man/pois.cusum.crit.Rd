\name{pois.cusum.crit}
\alias{pois.cusum.crit}
\title{Compute alarm thresholds and randomization constants of Poisson CUSUM control charts}
\description{Computation of the CUSUM upper limit and, if needed, of the randomization probability, given mean \code{mu0}.}
\usage{pois.cusum.crit(mu0, km, A, m, i0=0, sided="upper", rando=FALSE)}
\arguments{
\item{mu0}{actual in-control mean.}
\item{km}{enumerator of rational approximation of reference value \code{k}.}
\item{A}{target in-control ARL (average run length).}
\item{m}{denominator of rational approximation of reference value.}
\item{i0}{head start value as integer multiple of \code{1/m}; should be an element of \code{0:100} (a more reasonable upper limit will be established soon).}
\item{sided}{distinguishes between different one- and two-sided CUSUM control chart by choosing
\code{"upper"}, \code{"lower"} and \code{"two"}, respectively.}
\item{rando}{Switch for activating randomization in order to allow continuous ARL control.}
}
\details{
The monitored data follow a Poisson distribution with \code{mu} (here the in-control level \code{mu0}).
The ARL values of the resulting EWMA control chart are determined via Markov chain calculations.
With some grid search, we obtain the smallest value for the integer threshold component \code{hm} so that
the resulting ARL is not smaller than \code{A}. If equality is needed, then activating \code{rando=TRUE}
yields the corresponding randomization probability \code{gamma}.
More details will follow in a paper that will be submitted in 2020.
}
\value{Returns two single values, integer threshold \code{hm} resulting in the final
alarm threshold \code{h=hm/m}, and the randomization probability.}
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
## Lucas 1985
mu0 <- 0.25
km <- 1
A <- 430
m  <- 4
#cv <- pois.cusum.crit(mu0, km, A, m)
cv <- c(40, 0)
# Lucas reported h = 10 alias hm = 40 (in Table 2, first block, row 10.0 .25 .0 ..., column 1.0
# Recall that Lucas and other trigger an alarm, if the CUSUM statistic is greater than
# or equal to the alarm threshold h
print(cv)
}
\keyword{ts}
