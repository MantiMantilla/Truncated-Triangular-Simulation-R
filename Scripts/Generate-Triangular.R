## returns a function that returns a list with each of its elements defined as functions or vectors of one
## element that describe aspects of a triangular distribution (pdf, cdf, mean, median, mode, upper bound, lower bound, variance).

generate.triangular <- function(L, U, M) {
  ## create a function that returns the pdf of the triangular, given some x.
  pdf <- function(x) {
    if (x < L) {
      result <- 0
      return(result)
    } else if (x < M) {
      result <- 2 * (x - L) / ((U - L) * (M - L))
      return(result)
    } else if (x == M) {
      result <- 2 / (U - L)
      return(result)
    } else if (x <= U) {
      result <- 2 * (U - x) / ((U - L) * (U - M))
      return(result)
    } else if (U < x) {
      result <- 0
      return(result)
    }
  }
  
  ## create a function that returns the cdf of the triangular, given some x.
  cdf <- function(x) {
    if (x <= L) {
      result <- 0
      return(result)
    } else if (x <= M) {
      result <- ((x - L) ^ 2) / ((U - L) * (M - L))
      return(result)
    } else if (x < U) {
      result <- 1 - ((U - x) ^ 2) / ((U - L) * (U - M))
      return(result)
    } else if (U <= x) {
      result <- 1
      return(result)
    }
  }
  
  ## create a vector that describes the mean of the distribution
  tri.mean <- (L + U + M) / 3
  
  ## create a vector that describes the median of the distribution
  tri.median <- if (M >= (L + U) / 2) {
    L + sqrt(((U - L) * (M - L)) / 2)
  } else if (M < (L + U) / 2) {
    U - sqrt(((U - L) * (U - M)) / 2)
  }
  
  ## create a vector that describes the mode of the distribution
  tri.mode <- M
  
  ## create a vector that describes the upper bound of the distribution
  tri.upper <- U
  
  ## create a vector that describes the lower bound of the distribution
  tri.lower <- L
  
  ## create a vector that describes the variance of the distribution
  tri.var <-
    ((L ^ 2) + (U ^ 2) + (M ^ 2) - L * U - L * M - U * M) / 18
  
  return(
    list(
      pdf = pdf,
      cdf = cdf,
      tri.mean = tri.mean,
      tri.median = tri.median,
      tri.mode = tri.mode,
      tri.upper = tri.upper,
      tri.lower = tri.lower,
      tri.var = tri.var
    )
  )
}
