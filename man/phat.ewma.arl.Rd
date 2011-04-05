\name{phat.ewma.arl}
\alias{phat.ewma.arl}
\alias{phat.ewma.crit}
\alias{phat.ewma.lambda}
\title{Compute ARLs of EWMA phat control charts}
\description{Computation of the (zero-state) Average Run Length (ARL), upper control limit (ucl)
for given in-control ARL, and lambda for minimal out-of control ARL at given shift.}
\usage{phat.ewma.arl(lambda, ucl, mu, n, z0, sigma=1, LSL=-3, USL=3, N=15, qm=15)

phat.ewma.crit(lambda, L0, mu, n, z0, sigma=1, LSL=-3, USL=3, N=15, qm=15)

phat.ewma.lambda(L0, mu, n, z0, sigma=1, max_l=1, min_l=.001, LSL=-3, USL=3, qm=15)
}
\arguments{
\item{lambda}{smoothing parameter of the EWMA control chart.}
\item{ucl}{upper control limit of the EWMA phat control chart.}
\item{L0}{pre-defined in-control ARL (Average Run Length).}
\item{mu}{true mean or mean where the ARL should be minimized (then the in-control mean is simply 0).}
\item{n}{Batch size.}
\item{z0}{so-called headstart (give fast initial response).}
\item{sigma}{actual standard deviation of the data -- the in-control value is 1.}
\item{max_l, min_l}{maximal and minimal value for optimal lambda search.}
\item{LSL,USL}{lower and upper specification limit, respectively.}
\item{N}{size of collocation base, dimension of the resulting linear equation system is equal to \code{N}.}
\item{qm}{number of nodes for collocation quadratures.}
}
\details{
More details will follow.
}
\value{Return single values which resemble the ARL, the critical value, and the optimal lambda, respectively.}
\references{
S. Knoth, S. Steinmetz (2011),
EWMA \emph{p} charts under Sampling by Variables -- Ideas, Numerics and Properties,
\emph{Submitted to Technometrics}.
}
\author{Sven Knoth}
\seealso{
\code{sewma.arl} for a further collocation based ARL calculation routine.}
\examples{
## S. Knoth, S. Steinmetz (2011),

# Table 1

lambdas <- c(.5, .25, .2, .1)
L0 <- 370.4
n <- 5
LSL <- -3
USL <- 3

phat.ewma.CRIT <- Vectorize("phat.ewma.crit", "lambda")
p.star <- pnorm( LSL ) + pnorm( -USL ) ## lower bound of the chart
ucls <- phat.ewma.CRIT(lambdas, L0, 0, n, p.star, LSL=LSL, USL=USL)
print(cbind(lambdas, ucls))

# Table 2

mus <- c((0:4)/4, 1.5, 2, 3)
phat.ewma.ARL <- Vectorize("phat.ewma.arl", "mu")
arls <- NULL
for ( i in 1:length(lambdas) ) {
  arls <- cbind(arls, round(phat.ewma.ARL(lambdas[i], ucls[i], mus, n, p.star, LSL=LSL, USL=USL), digits=2))
}
arls <- data.frame(arls, row.names=NULL)
names(arls) <- lambdas
print(arls)

# Table 3

mus <- c(.25, .5, 1, 2)
phat.ewma.LAMBDA <- Vectorize("phat.ewma.lambda", "mu")
lambdas <- phat.ewma.LAMBDA(L0, mus, n, p.star, LSL=LSL, USL=USL)
print(cbind(mus, lambdas))
}
\keyword{ts}
