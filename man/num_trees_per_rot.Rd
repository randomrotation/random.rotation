% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_rotation.R
\name{num_trees_per_rot}
\alias{num_trees_per_rot}
\title{Compute number of trees per rotation}
\usage{
num_trees_per_rot(N, R, h, fn = weights_exp)
}
\arguments{
\item{N}{number of trees in total}

\item{R}{number of rotations in total}

\item{h}{value of tuning parameter}

\item{fn}{distribution function, currently in weights_{cut, exp, rre}}
}
\description{
Compute the correct number of trees to grow per rotation, given h and the distribution function
}
\examples{
num_trees_per_rot(5000, 100, 10, fn=weights_cut)
}
\keyword{random}
\keyword{rotation}
