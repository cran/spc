\name{xtshewhart.ar1.arl}
\alias{xtshewhart.ar1.arl}
\title{Compute ARLs of modified Shewhart control charts for AR(1) data with Student t residuals}
\description{Computation of the (zero-state) Average Run Length (ARL)
for modified Shewhart charts deployed to the original AR(1) data where the residuals follow a Student t distribution.}
\usage{xtshewhart.ar1.arl(alpha, cS, df, delta=0, N1=50, N2=30, N3=2*N2, INFI=10, mode="tan")}
\arguments{
\item{alpha}{lag 1 correlation of the data.}
\item{cS}{critical value (alias to alarm limit) of the Shewhart control chart.}
\item{df}{degrees of freedom -- parameter of the t distribution.}
\item{delta}{potential shift in the data (in-control mean is zero.}
\item{N1}{number of quadrature nodes for solving the ARL integral equation, dimension of the resulting linear equation system is \code{N1}.}
\item{N2}{second number of quadrature nodes for combining the probability density function of the first observation following the margin distribution and the solution of the ARL integral equation.}
\item{N3}{third number of quadrature nodes for solving the left eigenfunction integral equation to determine the margin density (see Andel/Hrach, 2000),
dimension of the resulting linear equation system is \code{N3}.}
\item{INFI}{substitute of \code{Inf} -- the left eigenfunction integral equation is defined on the whole real axis; now it is reduced to \code{(-INFI,INFI)}.}
\item{mode}{Controls the type of variables substitution that might improve the numerical performance. Currently, \code{"identity"}, \code{"sin"}, \code{"sinh"}, and \code{"tan"} (default) are provided.}
}
\details{
Following the idea of Schmid (1995), \code{1-alpha} times the data turns out to be an
EWMA smoothing of the underlying AR(1) residuals. Hence, by combining the solution of
the EWMA ARL integral equation and the stationary distribution of the AR(1) data
(here Student t distribution is assumed) one gets easily the overall ARL.
}
\value{It returns a single value resembling the zero-state ARL of a modified Shewhart chart.}
\references{
J. Andel, K. Hrach (2000).
On calculation of stationary density of autoregressive processes.
\emph{Kybernetika, Institute of Information Theory and Automation AS CR 36}(3), 311-319.

H. Kramer, W. Schmid (2000).
The influence of parameter estimation on the ARL of Shewhart type charts for time series.
\emph{Statistical Papers 41}(2), 173-196.

W. Schmid (1995).
On the run length of a Shewhart chart for correlated data. 
\emph{Statistical Papers 36}(1), 111-130.
}
\author{Sven Knoth}
\seealso{
\code{xtewma.arl} for zero-state ARL computation of EWMA control charts in case of Student t distributed data.
}
\examples{
## will follow
}
\keyword{ts}
