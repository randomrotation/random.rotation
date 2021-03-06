% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_rotation.R
\name{pre_sample_preprocessing}
\alias{pre_sample_preprocessing}
\title{Pre-processing of raw input file (BEFORE in-/out-of-sample split)}
\usage{
pre_sample_preprocessing(
  df,
  target,
  tn_in = 0.6,
  th_in = 0.9,
  tl_in = 0.85,
  tf_in = 0.075
)
}
\arguments{
\item{df}{data frame}

\item{target}{column name or index of the response}

\item{tn_in}{threshold to determine if a column should be numeric}

\item{th_in}{threshold to determine if a string or factor column contains too many unique values and should be dropped}

\item{tl_in}{threshold to determine if a specific string or factor is repeated too many times and the column should be dropped}

\item{tf_in}{threshold to determine if a column should be converted to a factor}
}
\description{
This function takes a data frame and response column (target) and performs the following tasks:
(1) For string (character) columns it checks if many of them are numbers and converts the column to numeric if necessary
(2) For columns with very few unique values, it converts these columns to factors
(3) For string (character) or factor columns it checks if most rows are unique and drops columns for which this is the case
(4) For string (character) or factor columns it checks if one specific string or factor is very frequently present and drops those columns
}
\details{
Additional pre-processing functions could be added here.

It should be noted that these steps occur BEFORE the file is split into training and testing data and as such care must be taken to
operations that lead to a bias when used to train a classifier (e.g. ranking of all values, imputation based on all values, etc).

Also see: post_sample_preprocessing() for a function that only gets called AFTER the file is split into training and testing data.

The trade-off is speed vs bias. Putting everything here leads to a significant speed-up but with a potential bias and vice versa.
}
\examples{
pre_sample_preprocessing <- function(df, 5)
}
\keyword{pre-processing}
