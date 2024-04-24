\name{dphat}
\encoding{utf-8}
\alias{dphat}
\alias{pphat}
\alias{qphat}
\title{Percent defective for normal samples}
\description{Density, distribution function and quantile function
for the sample percent defective calculated on normal samples
with mean equal to \code{mu} and standard deviation equal to \code{sigma}.}
\usage{dphat(x, n, mu=0, sigma=1, type="known", LSL=-3, USL=3, nodes=30)

pphat(q, n, mu=0, sigma=1, type="known", LSL=-3, USL=3, nodes=30)

qphat(p, n, mu=0, sigma=1, type="known", LSL=-3, USL=3, nodes=30)}
\arguments{
\item{x, q}{vector of quantiles.}
\item{p}{vector of probabilities.}
\item{n}{sample size.}
\item{mu, sigma}{parameters of the underlying normal distribution.}
\item{type}{choose whether the standard deviation is given and fixed (\code{"known"}) or estimated and potententially monitored (\code{"estimated"}).}
\item{LSL,USL}{lower and upper specification limit, respectively.}
\item{nodes}{number of quadrature nodes needed for \code{type="estimated"}.}
}
\details{Bruhn-Suhr/Krumbholz (1990) derived the cumulative distribution function
of the sample percent defective calculated on normal samples to applying them for a new variables sampling plan.
These results were heavily used in Krumbholz/Zöller (1995) for Shewhart and in Knoth/Steinmetz (2013) for EWMA control charts.
For algorithmic details see, essentially, Bruhn-Suhr/Krumbholz (1990).
Two design variants are treated: The simple case, \code{type="known"}, with known normal variance and the presumably much
more relevant and considerably intricate case, \code{type="estimated"}, where both parameters of
the normal distribution are unknown. Basically, given lower and upper specification limits and the normal distribution,
one estimates the expected yield based on a normal sample of size \code{n}.
}
\value{Returns vector of pdf, cdf or qf values for the statistic phat.}
\references{
M. Bruhn-Suhr and W. Krumbholz (1990),
A new variables sampling plan for normally distributed lots with unknown standard deviation and double specification limits,
\emph{Statistical Papers} 31(1), 195-207.

W. Krumbholz and A. Zöller (1995),
\code{p}-Karten vom Shewhartschen Typ für die messende Prüfung,
\emph{Allgemeines Statistisches Archiv} 79, 347-360.

S. Knoth and S. Steinmetz (2013),
EWMA \code{p} charts  under sampling by variables,
\emph{International Journal of Production Research} 51(13), 3795-3807.
}
\author{Sven Knoth}
\seealso{
\code{phat.ewma.arl} for routines using the herewith considered phat statistic.}
\examples{
# Figures 1 (c) and (d) from Knoth/Steinmetz (2013)
n      <-  5
LSL    <- -3
USL    <-  3

par(mar=c(5, 5, 1, 1) + 0.1)

p.star <- 2*pnorm( (LSL-USL)/2 ) # for p <= p.star pdf and cdf vanish

p_ <- seq(p.star+1e-10, 0.07, 0.0001) # define support of Figure 1

# Figure 1 (c)
pp_ <- pphat(p_, n)
plot(p_, pp_, type="l", xlab="p", ylab=expression(P( hat(p) <= p )),
     xlim=c(0, 0.06), ylim=c(0,1), lwd=2)
abline(h=0:1, v=p.star, col="grey")

# Figure 1 (d)
dp_ <- dphat(p_, n)
plot(p_, dp_, type="l", xlab="p", ylab="f(p)", xlim=c(0, 0.06),
     ylim=c(0,50), lwd=2)
abline(h=0, v=p.star, col="grey")
}
\keyword{ts}
