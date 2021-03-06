% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_rotation.R
\name{uni_rank}
\alias{uni_rank}
\title{Apply map}
\usage{
uni_rank(v, v_map)
}
\arguments{
\item{v}{numeric input vector}

\item{v_map}{list of inputs and corresponding mapped values}
}
\description{
Use an in-sample rank map (v_map) to obtain the ranks of an out-of-sample input vector v
}
\examples{
uni_rank(c(1.7, 2.1, 3.5), zero_one_rank_map(0, 1, 4, 9))
}
\keyword{ranking}
