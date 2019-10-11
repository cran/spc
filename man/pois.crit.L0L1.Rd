\name{pois.cusum.crit.L0L1}
\alias{pois.cusum.crit.L0L1}
\title{Compute the CUSUM k and h for given in-control ARL L0 and out-of-control ARL L1, Poisson case}
\description{Computation of the reference value k and the alarm threshold h
for one-sided CUSUM control charts monitoring Poisson data, if the in-control ARL L0 and the out-of-control ARL L1 are given.}
\usage{pois.cusum.crit.L0L1(mu0, L0, L1, sided="upper", OUTPUT=FALSE)}
\arguments{
\item{mu0}{in-control Poisson mean.}
\item{L0}{in-control ARL.}
\item{L1}{out-of-control ARL.}
\item{sided}{distinguishes between \code{"upper"} and \code{"lower"} CUSUM designs.}
\item{OUTPUT}{controls whether iteration details are printed.}
}
\details{
\code{pois.cusum.crit.L0L1} determines the reference value k and the alarm threshold h
for given in-control ARL \code{L0} and out-of-control ARL \code{L1}
by applying grid search and using \code{pois.cusum.arl()} and \code{pois.cusum.crit()}.
These CUSUM design rules were firstly (and quite rarely afterwards) used by Ewan and Kemp.
In the Poisson case, Rossi et al. applied them while analyzing three different normal
approximations of the Poisson distribution. See the example which illustrates
the validity of all these approaches.
}
\value{Returns a data frame with results for the denominator \code{m} of the rational approximation,
\code{km} as (integer) enumerator of the reference value approximation, the resulting
out-of-control mean \code{mu1}, the final approximation \code{k} of the reference value,
the threshold values \code{hm} (integer) and \code{h} (\code{=hm/m}), and the randomization constant
\code{gamma} (the target in-control ARL is exactly matched).}
\references{
W. D. Ewan and K. W. Kemp (1960),
Sampling inspection of continuous processes with no autocorrelation between successive results,
\emph{Biometrika} 47 (3/4), 363-380.

K. W. Kemp (1962),
The Use of Cumulative Sums for Sampling Inspection Schemes,
\emph{Journal of the Royal Statistical Sociecty C, Applied Statistics} 11(1), 16-31.

G. Rossi, L. Lampugnani and M. Marchi (1999),
An approximate CUSUM procedure for surveillance of health events,
\emph{Statistics in Medicine} 18(16), 2111-2122.
}
\author{Sven Knoth}
\seealso{\code{pois.cusum.arl} for zero-state ARL and \code{pois.cusum.crit} for threshold h computation.}
\examples{
## Table 1 from Rossi et al. (1999) -- one-sided CUSUM
La <- 500 # in-control ARL
Lr <- 7 # out-of-control ARL
m_a <- 0.52 # in-control mean of the Poisson variate
\dontrun{kh <- xcusum.crit.L0L1(La, Lr, sided="one")} # instead of deploying EK1960, one could use more accurate numbers
EK_k <- 0.60 # EK1960 results in
EK_h <- 3.80 # Table 2 on p. 372
eZR <- 2*EK_h # reproduce normal ooc mean from reference value k
m_r <- 1.58 # EK1960 Table 3 on p. 377 for m_a = 0.52
R1 <- round( eZR/sqrt(m_a) + 1, digits=2)
R2 <- round( ( eZR/2/sqrt(m_a) + 1 )^2, digits=2)
R3 <- round(( sqrt(4 + 2*eZR/sqrt(m_a)) - 1 )^2, digits=2)
RS <- round( m_r / m_a, digits=2 )
\dontrun{K_hk <- pois.cusum.crit.L0L1(m_a, La, Lr)} # 'our' 'exact' approach
K_hk <- data.frame(m=1000, km=948, mu1=1.563777, k=0.948, hm=3832, h=3.832, gamma=0.1201901)
# get k for competing means mu0 (m_a) and mu1 (m_r)
k_m01 <- function(mu0, mu1) (mu1 - mu0) / (log(mu1) - log(mu0))
# get ooc mean mu1 (m_r) for given mu0 (m_a) and reference value k
m1_km0 <- function(mu0, k) {
  zero <- function(x) k - k_m01(mu0,x)
  upper <- mu0 + .5
  while ( zero(upper) > 0 ) upper <- upper + 0.5
  mu1 <- uniroot(zero, c(mu0*1.00000001, upper), tol=1e-9)$root
  mu1
}
K_m_r <- m1_km0(m_a, K_hk$k)
RK <- round( K_m_r / m_a, digits=2 )
cat(paste(m_a, R1, R2, R3, RS, RK, "\n", sep="\t"))
}
\keyword{ts}
