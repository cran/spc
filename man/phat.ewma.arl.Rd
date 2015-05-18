\name{phat.ewma.arl}
\alias{phat.ewma.arl}
\alias{phat.ewma.crit}
\alias{phat.ewma.lambda}
\title{Compute ARLs of EWMA phat control charts}
\description{Computation of the (zero-state) Average Run Length (ARL), upper control limit (ucl)
for given in-control ARL, and lambda for minimal out-of control ARL at given shift.}
\usage{phat.ewma.arl(lambda, ucl, mu, n, z0, sigma=1, type="known", LSL=-3, USL=3, N=15,
qm=25, ntype="coll")

phat.ewma.crit(lambda, L0, mu, n, z0, sigma=1, type="known", LSL=-3, USL=3, N=15, qm=25)

phat.ewma.lambda(L0, mu, n, z0, sigma=1, type="known", max_l=1, min_l=.001, LSL=-3, USL=3,
qm=25)
}
\arguments{
\item{lambda}{smoothing parameter of the EWMA control chart.}
\item{ucl}{upper control limit of the EWMA phat control chart.}
\item{L0}{pre-defined in-control ARL (Average Run Length).}
\item{mu}{true mean or mean where the ARL should be minimized (then the in-control mean is simply 0).}
\item{n}{subgroup size.}
\item{z0}{so-called headstart (gives fast initial response).}
\item{type}{choose whether the standard deviation is given and fixed (\code{"known"}) or estimated and potentially monitored (\code{"estimated"}).}
\item{sigma}{actual standard deviation of the data -- the in-control value is 1.}
\item{max_l, min_l}{maximal and minimal value for optimal lambda search.}
\item{LSL,USL}{lower and upper specification limit, respectively.}
\item{N}{size of collocation base, dimension of the resulting linear equation system is equal to \code{N}.}
\item{qm}{number of nodes for collocation quadratures.}
\item{ntype}{switch between the default method \code{coll} (collocation) and the classic one \code{markov} (Markov chain approximation)
for calculating the ARL numerically.}
}
\details{
The three implemented functions allow to apply a new type control chart. Basically, lower and upper
specification limits are given. The monitoring
vehicle then is the empirical probability that an item will not follow these specification given the
sequence of sample means. If
the related EWMA sequence violates the control limits, then the alarm indicates a significant process
deterioration. For details see the
paper mentioned in the references. To be able to construct the control charts, see the first example.
}
\value{Return single values which resemble the ARL, the critical value, and the optimal lambda, respectively.}
\references{
S. Knoth and S. Steinmetz (2013),
EWMA \code{p} charts  under sampling by variables,
\emph{International Journal of Production Research} 51, 3795-3807.
}
\author{Sven Knoth}
\seealso{
\code{sewma.arl} for a further collocation based ARL calculation routine.}
\examples{
## Simple example to demonstrate the chart.

# some functions
h.mu <- function(mu) pnorm(LSL-mu) + pnorm(mu-USL)
ewma <- function(x, lambda=0.1, z0=0) filter(lambda*x, 1-lambda, m="r", init=z0)

# parameters
LSL <- -3       # lower specification limit
USL <-  3	# upper specification limit
n <- 5		# batch size
lambda <- 0.1	# EWMA smoothing parameter
L0 <- 1000	# in-control Average Run Length (ARL)
z0 <- h.mu(0)	# start at minimal defect level
ucl <- phat.ewma.crit(lambda, L0, 0, n, z0, LSL=LSL, USL=USL)

# data
x0 <- matrix(rnorm(50*n), ncol=5)	# in-control data
x1 <- matrix(rnorm(50*n, mean=0.5), ncol=5)# out-of-control data
x <- rbind(x0,x1)			# all data

# create chart
xbar <- apply(x, 1, mean)
phat <- h.mu(xbar)
z <- ewma(phat, lambda=lambda, z0=z0)
plot(1:length(z), z, type="l", xlab="batch", ylim=c(0,.02))
abline(h=z0, col="grey", lwd=.7)
abline(h=ucl, col="red")


## S. Knoth, S. Steinmetz (2013)

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
  arls <- cbind(arls, round(phat.ewma.ARL(lambdas[i], ucls[i], mus,
                n, p.star, LSL=LSL, USL=USL), digits=2))
}
arls <- data.frame(arls, row.names=NULL)
names(arls) <- lambdas
print(arls)

# Table 3

\dontrun{
mus <- c(.25, .5, 1, 2)
phat.ewma.LAMBDA <- Vectorize("phat.ewma.lambda", "mu")
lambdas <- phat.ewma.LAMBDA(L0, mus, n, p.star, LSL=LSL, USL=USL)
print(cbind(mus, lambdas))}
}
\keyword{ts}
