\name{xewma.q}
\alias{xewma.q}
\title{Compute RL quantiles of EWMA control charts}
\description{Computation of quantiles of the Run Length (RL)
for EWMA control charts monitoring normal mean.}
\usage{xewma.q(l, c, mu, p, zr=0, hs=0, sided="one", r=40)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{mu}{true mean.}
\item{p}{quantile level.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (give fast initial response).}
\item{sided}{distinguish between one- and two-sided EWMA control chart 
by choosing \code{"one"} and \code{"two"}, respectively.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-sided) or \code{r}
(two-sided).}
}
\details{
Instead of the popular ARL (Average Run Length) quantiles of the EWMA
stopping time (Run Length) are determined. The algorithm is based on
Waldmann's survival function iteration procedure.
}
\value{Returns a single value which resembles the ARL.}
\references{
K.-H. Waldmann (1986),
Bounds for the distribution of the run length of geometric moving 
average charts, \emph{Appl. Statist. 35}, 151-158.

F. F. Gan (1993),
An optimal design of EWMA control charts based on the median run length,
\emph{J. Stat. Comput. Simulation 45}, 169-184.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} for zero-state ARL computation of EWMA control charts.
}
\examples{
## Gan (1993), two-sided EWMA with fixed control limits
## some values of his Table 1 -- any median RL should be 500
XEWMA.Q <- Vectorize("xewma.q", c("l", "c"))
G.lambda <- c(.05, .1, .15, .2, .25)
G.h <- c(.441, .675, .863, 1.027, 1.177)
MEDIAN <- XEWMA.Q(G.lambda, G.h/sqrt(G.lambda/(2-G.lambda)), 0, .5, sided="two")
print(cbind(G.lambda, MEDIAN))

## increase accuracy of thresholds

# (i) calculate threshold for given in-control median value by
#     deplyoing secant rule
xewma.c.of.quantile <- function(l, L0, mu, p, zr=0, hs=0, sided="one", r=40) {
  c2 <- 0
  a2 <- 0
  while ( a2 < L0 ) {
    c2 <- c2 + .5
    a2 <- xewma.q(l, c2, mu, p, zr=zr, hs=hs, sided=sided, r=r)
  }
  c1 <- c2 - .5
  a1 <- xewma.q(l, c1, mu, p, zr=zr, hs=hs, sided=sided, r=r)
  a.error <- 1; c.error <- 1
  while ( a.error>1e-6 && c.error>1e-12 ) {
    c3 <- c1 + (L0 - a1)/(a2 - a1)*(c2 - c1)
    a3 <- xewma.q(l, c3, mu, p, zr=zr, hs=hs, sided=sided, r=r)
    c1 <- c2; c2 <- c3
    a1 <- a2; a2 <- a3
    a.error <- abs(a2 - L0); c.error <- abs(c2 - c1)
  }
  c3
}
XEWMA.c.of.quantile <- Vectorize("xewma.c.of.quantile", "l")

# (ii) re-calculate the thresholds and remove the standarization step
L0 <- 500
G.h.new <- XEWMA.c.of.quantile(G.lambda, L0, 0, .5, sided="two")
G.h.new <- round(G.h.new * sqrt(G.lambda/(2-G.lambda)), digits=5)

# (iii) compare Gan's original values and the new ones woth 5 digits
print(cbind(G.lambda, G.h.new, G.h))

# (iv) calculate the new medians
MEDIAN <- round(XEWMA.Q(G.lambda, G.h.new/sqrt(G.lambda/(2-G.lambda)), 0, .5, sided="two"), digits=1)
print(cbind(G.lambda, MEDIAN))
}
\keyword{ts}
