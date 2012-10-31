\name{p.ewma.arl}
\alias{p.ewma.arl}
\title{Compute ARLs of binomial EWMA p control charts}
\description{Computation of the (zero-state) Average Run Length (ARL) at given rate \code{p}.}
\usage{p.ewma.arl(lambda, ucl, n, p, z0, d.res=1, r.mode="ieee.round", i.mode="integer")}
\arguments{
\item{lambda}{smoothing parameter of the EWMA p control chart.}
\item{ucl}{upper control limit of the EWMA p control chart.}
\item{n}{batch size.}
\item{p}{(failure/success) rate.}
\item{z0}{so-called headstart (give fast initial response).}
\item{d.res}{resolution (see details).}
\item{r.mode}{round mode -- allowed modes are \code{"gan.floor"}, \code{"floor"}, \code{"ceil"}, \code{"ieee.round"}, \code{"round"}, \code{"mix"}.}
\item{i.mode}{type of interval center -- \code{"integer"} or \code{"half"} integer.}
}
\details{
The monitored data follow a binomial distribution with size \code{n} and failure/success probability \code{p}.
The ARL values of the resulting EWMA control chart are determined by a Markov chain approximation. Here, the original EWMA values are approximated by
multiples of one over \code{d.res}. Different ways of rounding (see \code{r.mode}) to the next multiple are implemented.
Besides Gan's paper nothing is published about the numerical subtleties.
}
\value{Return single value which resemble the ARL.}
\references{
F. F. Gan (1990),
Monitoring observations generated from a binomial distribution using modified exponentially weighted moving average control chart,
\emph{J. Stat. Comput. Simulation} 37, 45-60.

S. Knoth and S. Steinmetz (2012),
EWMA \code{p} charts  under sampling by variables,
accepted in \emph{International Journal of Production Research}.
}
\author{Sven Knoth}
\seealso{later.}
\examples{
## Gan (1990)

# Table 1

n <- 150
p0 <- .1
z0 <- n*p0

lambda <- c(1, .51, .165)
hu <- c(27, 22, 18)

p.value <- .1 + (0:20)/200

p.EWMA.arl <- Vectorize(p.ewma.arl, "p")

arl1.value <- round(p.EWMA.arl(lambda[1], hu[1], n, p.value, z0, r.mode="round"), digits=2)
arl2.value <- round(p.EWMA.arl(lambda[2], hu[2], n, p.value, z0, r.mode="round"), digits=2)
arl3.value <- round(p.EWMA.arl(lambda[3], hu[3], n, p.value, z0, r.mode="round"), digits=2)

arls <- matrix(c(arl1.value, arl2.value, arl3.value), ncol=length(lambda))
rownames(arls) <- p.value
colnames(arls) <- paste("lambda =", lambda)
arls

## Knoth/Steinmetz (2012)

n <- 5
p0 <- 0.02
z0 <- n*p0
lambda <- 0.3
ucl <- 0.649169922 ## should provide in-control ARL of about 370.4; determined with d.res = 2^14 = 16384

res.list <- 2^(1:12)
arl.list <- NULL
for ( res in res.list ) {
  arl <- p.ewma.arl(lambda, ucl, n, p0, z0, d.res=res)
  arl.list <- c(arl.list, arl)
}
cbind(res.list, arl.list)
}
\keyword{ts}
