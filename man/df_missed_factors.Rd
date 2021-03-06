% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_rotation.R
\name{df_missed_factors}
\alias{df_missed_factors}
\title{Auto-detect non-factor columns that should be factors}
\usage{
df_missed_factors(df, tf = 0.075)
}
\arguments{
\item{df}{data frame}

\item{tf}{threshold for percentage of unique values under which it should be a factor}
}
\description{
Returns columns for which there are very few unique values (should be factor)
}
\examples{
df_missed_factors(df)
}
\keyword{pre-processing}
