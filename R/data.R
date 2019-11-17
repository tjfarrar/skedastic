#' Probability distribution for number of peaks in a continuous,
#' uncorrelated stochastic series
#'
#' A dataset containing the probability mass function for the distribution of
#' the number of peaks in a continuous, uncorrelated stochastic series. These
#' probabilities were generated from the \code{\link{dpeak}} function.
#' This function is computationally very slow for \eqn{n > 170}; thus the
#' functions of \code{skedastic} package that require peak probabilities
#' (\code{\link{ppeak}} and \code{\link{goldfeld_quandt}}) by default
#' obtain the probabilities from this data set rather than from
#' \code{\link{dpeak}}, provided that \eqn{n \leq 500}.
#'
#' @format A list of 500 objects. The \eqn{n}th object is a double vector
#' of length \eqn{n}, with elements representing the probability of \eqn{k}
#' peaks for \eqn{k=0,1,\ldots,n-1}.
"dpeakdat"
