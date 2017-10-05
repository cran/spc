\name{xshewhart.ar1.arl}
\alias{xshewhart.ar1.arl}
\title{Compute ARLs of modified Shewhart control charts for AR(1) data}
\description{Computation of the (zero-state) Average Run Length (ARL)
for modified Shewhart charts deployed to the original AR(1) data.}
\usage{xshewhart.ar1.arl(alpha, cS, delta=0, N1=50, N2=30)}
\arguments{
\item{alpha}{lag 1 correlation of the data.}
\item{cS}{critical value (alias to alarm limit) of the Shewhart control chart.}
\item{delta}{potential shift in the data (in-control mean is zero.}
\item{N1}{number of quadrature nodes for solving the ARL integral equation, dimension of the resulting linear equation system is \code{N1}.}
\item{N2}{second number of quadrature nodes for combining the probability density function of the first observation following the margin distribution and the solution of the ARL integral equation.}
}
\details{
Following the idea of Schmid (1995), \code{1- alpha} times the data turns out to be an
EWMA smoothing of the underlying AR(1) residuals. Hence, by combining the solution of
the EWMA ARL integral equation and the stationary distribution of the AR(1) data
(normal distribution is assumed) one gets easily the overall ARL.
}
\value{It returns a single value resembling the zero-state ARL of a modified Shewhart chart.}
\references{
S. Knoth, W. Schmid (2004).
Control charts for time series: A review.
In \emph{Frontiers in Statistical Quality Control 7},
edited by H.-J. Lenz, P.-T. Wilrich, 210-236, Physica-Verlag.

H. Kramer, W. Schmid (2000).
The influence of parameter estimation on the ARL of Shewhart type charts for time series.
\emph{Statistical Papers 41}(2), 173-196.

W. Schmid (1995).
On the run length of a Shewhart chart for correlated data. 
\emph{Statistical Papers 36}(1), 111-130.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} for zero-state ARL computation of EWMA control charts.
}
\examples{
## Table 1 in Kramer/Schmid (2000)

cS <- 3.09023
a  <- seq(0, 4, by=.5)
row1 <- row2 <- row3 <- NULL
for ( i in 1:length(a) ) {
  row1 <- c(row1, round(xshewhart.ar1.arl( 0.4, cS, delta=a[i]), digits=2))
  row2 <- c(row2, round(xshewhart.ar1.arl( 0.2, cS, delta=a[i]), digits=2))
  row3 <- c(row3, round(xshewhart.ar1.arl(-0.2, cS, delta=a[i]), digits=2))
}

results <- rbind(row1, row2, row3)
results

# original values are
# row1 515.44 215.48 61.85 21.63 9.19 4.58 2.61 1.71 1.29
# row2 502.56 204.97 56.72 19.13 7.95 3.97 2.33 1.59 1.25
# row3 502.56 201.41 54.05 17.42 6.89 3.37 2.03 1.46 1.20
}
\keyword{ts}
