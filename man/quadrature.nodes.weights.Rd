\name{quadrature.nodes.weights}
\alias{quadrature.nodes.weights}
\title{Calculate quadrature nodes and weights}
\description{Computation of the nodes and weights to enable numerical quadrature.}
\usage{quadrature.nodes.weights(n, type="GL", x1=-1, x2=1)}
\arguments{
\item{n}{number of nodes (and weights).}
\item{type}{quadrature type -- currently Gauss-Legendre, \code{"GL"}, and Radau, \code{"Ra"}, are supported.}
\item{x1}{lower limit of the integration interval.}
\item{x2}{upper limit of the integration interval.}
}
\details{
A more detailed description will follow soon. The algorithm for the Gauss-Legendre quadrature was delivered by
Knut Petras to me, while the one for the Radau quadrature was taken from John Burkardt.
}
\value{Returns two vectors which hold the needed quadrature nodes and weights.}
\references{
H. Brass and K. Petras (2011),
\emph{Quadrature Theory. The Theory of Numerical Integration on a Compact Interval,}
Mathematical Surveys and Monographs, American Mathematical Society.
}
\author{Sven Knoth}
\seealso{
Many of the ARL routines use the Gauss-Legendre nodes.
}
\examples{
# GL
n <- 10
qnw <-quadrature.nodes.weights(n, type="GL")
qnw

# Radau
n <- 10
qnw <-quadrature.nodes.weights(n, type="Ra")
qnw
}
\keyword{ts}
