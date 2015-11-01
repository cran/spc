\name{xtcusum.arl}
\alias{xtcusum.arl}
\title{Compute ARLs of CUSUM control charts}
\description{Computation of the (zero-state) Average Run Length (ARL)
for different types of CUSUM control charts monitoring normal mean.}
\usage{xtcusum.arl(k, h, df, mu, hs = 0, sided="one", mode="tan", r=30)}
\arguments{
\item{k}{reference value of the CUSUM control chart.}
\item{h}{decision interval (alarm limit, threshold) of the CUSUM control chart.}
\item{df}{degrees of freedom -- parameter of the t distribution.}
\item{mu}{true mean.}
\item{hs}{so-called headstart (give fast initial response).}
\item{sided}{distinguish between one- and  two-sided CUSUM schemes by choosing \code{"one"} and \code{"two"}, respectively.}
\item{r}{number of quadrature nodes, dimension of the resulting linear equation system is equal to \code{r+1}.}
\item{mode}{Controls the type of variables substitution that might improve the numerical performance. Currently, \code{"identity"}, \code{"sin"}, \code{"sinh"}, and \code{"tan"} (default) are provided.}
}
\details{
\code{xtcusum.arl} determines the Average Run Length (ARL) by numerically
solving the related ARL integral equation by means of the Nystroem method
based on Gauss-Legendre quadrature.
}
\value{Returns a single value which resembles the ARL.}
\references{
A. L. Goel, S. M. Wu (1971),
Determination of A.R.L. and a contour nomogram for CUSUM charts to
control normal mean, \emph{Technometrics 13}, 221-230.

D. Brook, D. A. Evans (1972),
An approach to the probability distribution of cusum run length,
\emph{Biometrika 59}, 539-548.

J. M. Lucas, R. B. Crosier (1982),
Fast initial response for cusum quality-control schemes: 
Give your cusum a headstart, \emph{Technometrics 24}, 199-205.

L. C. Vance (1986), 
Average run lengths of cumulative sum control charts for controlling 
normal means, \emph{Journal of Quality Technology 18}, 189-193.

K.-H. Waldmann (1986),
Bounds for the distribution of the run length of one-sided and
two-sided CUSUM quality control schemes,
\emph{Technometrics 28}, 61-67.

R. B. Crosier (1986),
A new two-sided cumulative quality control scheme,
\emph{Technometrics 28}, 187-194.
}
\author{Sven Knoth}
\seealso{
\code{xtewma.arl} for zero-state ARL computation of EWMA control charts and \code{xtcusum.arl} for the zero-state ARL of CUSUM for normal data.
}
\examples{
## will follow
}
\keyword{ts}
