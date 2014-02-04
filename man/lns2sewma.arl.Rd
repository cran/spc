\name{lns2ewma.arl}
\alias{lns2ewma.arl}
\title{Compute ARLs of EWMA ln \eqn{S^2}{S^2} control charts (variance charts)}
\description{Computation of the (zero-state) Average Run Length (ARL)
for different types of EWMA control charts 
(based on the log of the sample variance \eqn{S^2}) monitoring normal variance.}
\usage{lns2ewma.arl(l,cl,cu,sigma,df,hs=NULL,sided="upper",r=40)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{cl}{lower control limit of the EWMA control chart.}
\item{cu}{upper control limit of the EWMA control chart.}
\item{sigma}{true standard deviation.}
\item{df}{actual degrees of freedom, corresponds to subsample size (for known mean it is equal to the subsample size,
for unknown mean it is equal to subsample size minus one.}
\item{hs}{so-called headstart (enables fast initial response) -- the default value (hs=NULL) corresponds to the in-control
mean of ln \eqn{S^2}{S^2}.}
\item{sided}{distinguishes between one- and two-sided two-sided EWMA-\eqn{S^2}{S^2} control charts 
by choosing \code{"upper"} (upper chart with reflection at \code{cl}), \code{"lower"} (lower chart with reflection at \code{cu}),
and \code{"two"} (two-sided chart), respectively.} 
\item{r}{dimension of the resulting linear equation system: the larger the better.}
}
\details{
\code{lns2ewma.arl} determines the Average Run Length (ARL) by numerically
solving the related ARL integral equation by means of the Nystroem method
based on Gauss-Legendre quadrature.}
\value{Returns a single value which resembles the ARL.}
\references{
S. V. Crowder and M. D. Hamilton (1992),
An EWMA for monitoring a process standard deviation,
\emph{Journal of Quality Technology 24}, 12-21.

S. Knoth (2005),
Accurate ARL computation for EWMA-\eqn{S^2}{S^2} control charts,
\emph{Statistics and Computing 15}, 341-352.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} for zero-state ARL computation of EWMA control charts for monitoring normal mean.
}
\examples{
lns2ewma.ARL <- Vectorize("lns2ewma.arl", "sigma")

## Crowder/Hamilton (1992)
## moments of ln S^2
E_log_gamma <- function(df) log(2/df) + digamma(df/2)
V_log_gamma <- function(df) trigamma(df/2)
E_log_gamma_approx <- function(df) -1/df - 1/3/df^2 + 2/15/df^4
V_log_gamma_approx <- function(df) 2/df + 2/df^2 + 4/3/df^3 - 16/15/df^5

## results from Table 3 ( upper chart with reflection at 0 = log(sigma0=1) )
## original entries are (lambda = 0.05, K = 1.06, df=n-1=4)
# sigma   ARL
# 1       200
# 1.1      43
# 1.2      18
# 1.3      11
# 1.4       7.6
# 1.5       6.0
# 2         3.2

df <- 4
lambda <- .05
K <- 1.06
cu <- K * sqrt( lambda/(2-lambda) * V_log_gamma_approx(df) )

sigmas <- c(1 + (0:5)/10, 2)
arls <- round(lns2ewma.ARL(lambda, 0, cu, sigmas, df, hs=0, sided="upper"), digits=1)
data.frame(sigmas, arls)

## Knoth (2005)
## compare with Table 3 (p. 351)
lambda <- .05
df <- 4
K <- 1.05521
cu <- 1.05521 * sqrt( lambda/(2-lambda) * V_log_gamma_approx(df) )

## upper chart with reflection at sigma0=1 in Table 4
## original entries are
# sigma   ARL_0    ARL_-.267
# 1       200.0    200.0
# 1.1      43.04    41.55
# 1.2      18.10    19.92
# 1.3      10.75    13.11
# 1.4       7.63     9.93
# 1.5       5.97     8.11
# 2         3.17     4.67

M <- -0.267
cuM <- lns2ewma.crit(lambda, 200, df, cl=M, hs=M, r=60)[2]
arls1 <- round(lns2ewma.ARL(lambda, 0, cu, sigmas, df, hs=0, sided="upper"), digits=2)
arls2 <- round(lns2ewma.ARL(lambda, M, cuM, sigmas, df, hs=M, sided="upper", r=60), digits=2)
data.frame(sigmas, arls1, arls2)
}
\keyword{ts}
