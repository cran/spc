\name{xcusum.crit.L0L1}
\alias{xcusum.crit.L0L1}
\title{Compute the CUSUM k and h for given in-control ARL L0 and out-of-control L1}
\description{Computation of the reference value k and the alarm threshold h
for one-sided CUSUM control charts monitoring normal mean, if the in-control ARL L0 and
the out-of-control L1 are given.}
\usage{xcusum.crit.L0L1(L0, L1, hs=0, sided="one", r=30, L1.eps=1e-6, k.eps=1e-8)}
\arguments{
\item{L0}{in-control ARL.}
\item{L1}{out-of-control ARL.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one-, two-sided and Crosier's modified
two-sided CUSUM schemoosing \code{"one"}, \code{"two"}, and \code{"Crosier"}, respectively.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-, two-sided) or \code{2r+1} (Crosier).}
\item{L1.eps}{error bound for the L1 error.}
\item{k.eps}{bound for the difference of two successive values of k.}
}
\details{
\code{xcusum.crit.L0L1} determines the reference value k and the alarm threshold h
for given in-control ARL \code{L0} and out-of-control ARL \code{L1}
by applying secant rule and using \code{xcusum.arl()} and \code{xcusum.crit()}.
These CUSUM design rules were firstly (and quite rarely afterwards) used by Ewan and Kemp.
}
\value{Returns two values which resemble the reference value \code{k} and the threshold \code{h}.}
\references{
W. D. Ewan and K. W. Kemp (1960),
Sampling inspection of continuous processes with no autocorrelation between successive results,
\emph{Biometrika 47}, 363-380.

K. W. Kemp (1962),
The Use of Cumulative Sums for Sampling Inspection Schemes,
\emph{Journal of the Royal Statistical Sociecty C, Applied Statistics, 10}, 16-31.
}
\author{Sven Knoth}
\seealso{\code{xcusum.arl} for zero-state ARL and \code{xcusum.crit} for threshold h computation.}
\examples{
## Table 2 from Ewan/Kemp (1960) -- one-sided CUSUM
#
# A.R.L. at A.Q.L.   A.R.L. at A.Q.L.     k      h
#       1000                3           1.12   2.40
#       1000                7           0.65   4.06
#        500                3           1.04   2.26
#        500                7           0.60   3.80
#        250                3           0.94   2.11
#        250                7           0.54   3.51
#
L0.set <- c(1000, 500, 250)
L1.set <- c(3, 7)
cat("\nL0\tL1\tk\th\n")
for ( L0 in L0.set ) {
  for ( L1 in L1.set ) {
    result <- round(xcusum.crit.L0L1(L0, L1), digits=2)
    cat(paste(L0, L1, result[1], result[2], sep="\t"), "\n")
  }
}
#
# two confirmation runs
xcusum.arl(0.54, 3.51, 0) # Ewan/Kemp
xcusum.arl(result[1], result[2], 0) # here
xcusum.arl(0.54, 3.51, 2*0.54) # Ewan/Kemp
xcusum.arl(result[1], result[2], 2*result[1]) # here
#
## Table II from Kemp (1962) -- two-sided CUSUM
#
#    Lr                  k
#             La=250   La=500   La=1000
#    2.5       1.05     1.17     1.27
#    3.0       0.94     1.035    1.13
#    4.0       0.78     0.85     0.92
#    5.0       0.68     0.74     0.80
#    6.0       0.60     0.655    0.71
#    7.5       0.52     0.57     0.62
#   10.0       0.43     0.48     0.52
#
L0.set <- c(250, 500, 1000)
L1.set <- c(2.5, 3:6, 7.5, 10)
cat("\nL1\tL0=250\tL0=500\tL0=1000\n")
for ( L1 in L1.set ) {
  cat(L1)
  for ( L0 in L0.set ) {
    result <- round(xcusum.crit.L0L1(L0, L1, sided="two"), digits=2)
    cat("\t", result[1])
  }
  cat("\n")
}
#
}
\keyword{ts}
