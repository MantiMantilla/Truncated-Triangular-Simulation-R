## returns a function that returns a list with each of its elements defined as functions or vectors of one
## element that describe aspects of a truncated triangular distribution (pdf, cdf, mean, median, mode, upper bound, lower bound, variance).

generate.truncated.triangular <- function(a, b, orig.tri.dist) {
  M <- orig.tri.dist$tri.mode
  L <- orig.tri.dist$tri.lower
  U <- orig.tri.dist$tri.upper
  ## Throw an error if the new bound do not fall witin the old bounds of the triangular pdf.
  if(a < L | b > U) {
    stop("The new bounds of the pdf do not fall within the original range")
  }
  
  ## create a function that returns the pdf of the truncated triangular, given some x.
  pdf <- function(x) {
    if(a < x & x <= b){
      pdf.g <- orig.tri.dist$pdf
      cumul.F <- orig.tri.dist$cdf
      result <- pdf.g(x) / (cumul.F(b) - cumul.F(a))
      return(result)
    } else {
      result <- 0
      return(result)
    }
  }
  
  ## create a function that returns the cdf of the truncated triangular, given some x.
  cdf <- function(x) {
    cumul.F <- orig.tri.dist$cdf
    result <- (cumul.F(x) - cumul.F(a)) / (cumul.F(b) - cumul.F(a))
  }
  
  ## create a vector that describes the mean of the distribution
  trun.tri.mean <- if (a <= b & b < M) {
    result.numerator <- (2 * (- 3 * L * ((b ^ 2) / 2 - (a ^ 2) / 2) + b ^ 3 - a ^ 3)) / (3 * (U - L) * (M - L))
    result.denominator <- orig.tri.dist$cdf(b) - orig.tri.dist$cdf(a)
    result.numerator / result.denominator
  } else if (a < M & b >= M) {
    result.numerator <- (- M ^ 3 * U - 2 * a ^ 3 * U + 3 * L * a ^ 2 * U + 3 * M * b ^ 2 * U - 3 * L * b ^ 2 * U + L * M ^ 3 + 2 * M * a ^ 3 - 3 * L * M * a ^ 2 - 2 * M * b ^ 3 + 2 * L * b ^ 3) / (3 * (U - L) * (M - L) * (U - M))
    result.denominator <- orig.tri.dist$cdf(b) - orig.tri.dist$cdf(a)
    result.numerator / result.denominator
  } else if (M <= a & a <= b) {
    result.numerator <- (3 * b ^ 2 * U - 3 * a ^ 2 * U - 2 * b ^ 3 + 2 * a ^ 3) / (3 * (U - L) * (U - M))
    result.denominator <- orig.tri.dist$cdf(b) - orig.tri.dist$cdf(a)
    result.numerator / result.denominator
  } 
  
  ## create a vector that describes the median of the distribution
  # trun.tri.median <- if (c >= (a + b) / 2) {
  #   a + sqrt(((b - a) * (c - a)) / 2)
  # } else if (c < (a + b) / 2) {
  #   b - sqrt(((b - a) * (b - c)) / 2)
  # }
  
  ## create a vector that describes the mode of the distribution
  trun.tri.mode <- if (a <= b & b <= M) {
    b
  } else if(a <= b & b >= M) {
    M
  } else if(a >= M & b >= a) {
    a
  }
  
  ## create a vector that describes the upper bound of the distribution's range
  trun.tri.upper <- b
  
  ## create a vector that describes the lower bound of the distribution's range
  trun.tri.lower <- a
  
  ## create a vector that describes the variance of the distribution
  # trun.tri.var <-
  
  return(
    list(
      pdf = pdf,
      cdf = cdf,
      trun.tri.mean = trun.tri.mean,
      # trun.tri.median = trun.tri.median,
      trun.tri.mode = trun.tri.mode,
      trun.tri.upper = trun.tri.upper,
      trun.tri.lower = trun.tri.lower#,
      # trun.tri.var = trun.tri.var
    )
  )
}

## TODO: implement median and variance vector values.