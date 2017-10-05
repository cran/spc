\name{mewma.psi}
\alias{mewma.psi}
\title{Compute steady-state density of the MEWMA statistic}
\description{Computation of the (zero-state) steady-state density function of the statistic deployed in
multivariate exponentially weighted moving average (MEWMA) charts monitoring multivariate normal mean.}
\usage{mewma.psi(l, cE, p, type="cond", hs=0, r=20)}
\arguments{
\item{l}{smoothing parameter lambda of the MEWMA control chart.}
\item{cE}{alarm threshold of the MEWMA control chart.}
\item{p}{dimension of multivariate normal distribution.}
\item{type}{switch between \code{"cond"} and \code{"cycl"} for differentiating between the conditional
(no false alarm) and the cyclical (after false alarm re-start in \code{hs}), respectively.}
\item{hs}{the re-starting point for the cyclical steady-state framework.}
\item{r}{number of quadrature nodes.}
}
\details{Basically, ideas from Knoth (2015, MEWMA numerics) and Knoth (2014, steady-state ARL concepts) are merged.
More details will follow.}
\value{Returns a function.}
\references{
Sven Knoth (2014),
The Case Against the Use of Synthetic Control Charts,
accepted in \emph{Journal of Quality Technology}.

Sven Knoth (2017),
ARL Numerics for MEWMA Charts,
\emph{Journal of Quality Technology 49}, 78-89.

Sven Knoth (2017),
Calculating the Steady-State ARL of MEWMA Control Charts,
unpublished manuscript.
}
\author{Sven Knoth}
\seealso{
\code{mewma.arl} for calculating the in-control ARL of MEWMA.
}
\examples{
lambda <- 0.1
L0 <- 1000
p <- 3
h4 <- mewma.crit(lambda, L0, p)
x_ <- seq(0, h4*lambda/(2-lambda), by=0.002)
psi <- mewma.psi(lambda, h4, p)
psi_ <- psi(x_)
#plot(x_, psi_, type="l", xlab="x", ylab=expression(psi(x)))
}
\keyword{ts}
