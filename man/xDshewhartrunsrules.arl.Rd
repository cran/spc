\name{xDshewhartrunsrules.arl}
\alias{xDshewhartrunsrules.arl}
\alias{xDshewhartrunsrulesFixedm.arl}
\title{Compute ARLs of Shewhart control charts with and without runs rules
under drift}
\description{Computation of the zero-state Average Run Length (ARL)
under drift for Shewhart control charts with and without runs rules
monitoring normal mean.}
\usage{xDshewhartrunsrules.arl(delta, c = 1, m = NULL, type = "12")

xDshewhartrunsrulesFixedm.arl(delta, c = 1, m = 100, type = "12")
}
\arguments{
\item{delta}{true drift parameter.}
\item{c}{normalizing constant to ensure specific alarming behavior.}
\item{type}{controls the type of Shewhart chart used, seed details section.}
\item{m}{parameter of Gan's approach. If \code{m=NULL}, then \code{m} will increased until
the resulting ARL does not change anymore.}
}
\details{
Based on Gan (1991), the ARL is calculated for
Shewhart control charts with and without runs rules
under drift. The usual ARL function with mu=m*delta is determined and recursively via
m-1, m-2, ... 1 (or 0) the drift ARL determined.
\code{xDshewhartrunsrulesFixedm.arl} is the actual work horse, while
\code{xDshewhartrunsrules.arl} provides a convenience wrapper.
Note that Aerne et al. (1991) deployed a method that is
quite similar to Gan's algorithm. For \code{type} see
the help page of \code{xshewhartrunsrules.arl}.
}
\value{Returns a single value which resembles the ARL.}
\references{
F. F. Gan (1991),
EWMA control chart under linear drift,
\emph{J. Stat. Comput. Simulation 38}, 181-200.

L. A. Aerne, C. W. Champ and S. E. Rigdon (1991),
Evaluation of control charts under linear trend,
\emph{Commun. Stat., Theory Methods 20}, 3341-3349.
}
\author{Sven Knoth}
\seealso{
\code{xshewhartrunsrules.arl} for zero-state ARL computation of
Shewhart control charts with and without runs rules
for the classical step change model.
}
\examples{
## Aerne et al. (1991)
## Table I (continued)
## original numbers are
#     delta arl1of1 arl2of3 arl4of5  arl10
#  0.005623  136.67  120.90  105.34 107.08
#  0.007499  114.98  101.23   88.09  89.94
#  0.010000   96.03   84.22   73.31  75.23
#  0.013335   79.69   69.68   60.75  62.73
#  0.017783   65.75   57.38   50.18  52.18
#  0.023714   53.99   47.06   41.33  43.35
#  0.031623   44.15   38.47   33.99  36.00
#  0.042170   35.97   31.36   27.91  29.90
#  0.056234   29.21   25.51   22.91  24.86
#  0.074989   23.65   20.71   18.81  20.70
#  0.100000   19.11   16.79   15.45  17.29
#  0.133352   15.41   13.61   12.72  14.47
#  0.177828   12.41   11.03   10.50  12.14
#  0.237137    9.98    8.94    8.71  10.18
#  0.316228    8.02    7.25    7.26   8.45
#  0.421697    6.44    5.89    6.09   6.84
#  0.562341    5.17    4.80    5.15   5.48
#  0.749894    4.16    3.92    4.36   4.39
#  1.000000    3.35    3.22    3.63   3.52
c1of1 <- 3.069/3
c2of3 <- 2.1494/2
c4of5 <- 1.14
c10   <- 3.2425/3
DxDshewhartrunsrules.arl <- Vectorize(xDshewhartrunsrules.arl, "delta")
deltas <- 10^(-(18:0)/8)
arl1of1 <-
round(DxDshewhartrunsrules.arl(deltas, c=c1of1, type="1"), digits=2)
arl2of3 <-
round(DxDshewhartrunsrules.arl(deltas, c=c2of3, type="12"), digits=2)
arl4of5 <-
round(DxDshewhartrunsrules.arl(deltas, c=c4of5, type="13"), digits=2)
arl10 <- 
round(DxDshewhartrunsrules.arl(deltas, c=c10, type="SameSide10"), digits=2)
data.frame(delta=round(deltas, digits=6), arl1of1, arl2of3, arl4of5, arl10)
}
\keyword{ts}
