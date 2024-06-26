% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cppRSHalton_br}
\alias{cppRSHalton_br}
\title{Generate numbers from a Halton Sequence with a random start}
\usage{
cppRSHalton_br(
  n = 10L,
  bases = as.numeric(c()),
  seeds = as.numeric(c()),
  verbose = FALSE
)
}
\arguments{
\item{n}{Number of points required}

\item{bases}{Co-prime base for the Halton Sequence}

\item{seeds}{Random starting point in each dimension}

\item{verbose}{A boolean indicating whether informational messages are to be issued.}
}
\value{
Matrix with the columns, order of point, x in [0,1) and y in [0,1).
}
\description{
For efficiency, this function can generate points along a random start
Halton Sequence for a predefined Halton.
}
\details{
This function was first written in R by Paul van Dam-Bates for the
package BASMasterSample. Subsequently it was written in C/C++ by Phil Davies.
}
\examples{
# First 10 points in the Halton Sequence for base 2,3
 spbal::cppRSHalton_br(n = 10)
# First 10 points in the Halton Sequence for base 2,3 with
# starting point at the 15th and 22nd index.
 spbal::cppRSHalton_br(n = 10, seeds = c(14, 21))

}
