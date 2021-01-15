# Import relevant functions to generate triangular and truncated triangular distributions with cached properties (created as a list).
source(file.path("Scripts","Generate-Triangular.R"))
source(file.path("Scripts","Generate-Truncated-Triangular.R"))

# Create a list that contains the major properties of a triangular distribution with L = 2, U = 9, and M = 7
my.tri.dist <- generate.triangular(L = 2, U = 9, M = 7)

# Create a list that contains the major properties of a truncated triangular distribution between a = 3 and b = 8, based on the previous triangular.
my.trun.tri.dist <- generate.truncated.triangular(a = 3, b = 8, orig.tri.dist = my.tri.dist)

# Calculate the mean of the distribution
mu <- my.trun.tri.dist$trun.tri.mean

# simulate 1000 instances of a random variable with f(x) = truncated triangular.
vec <- rep(0, 1000)
success.count <- 0
while (success.count < 1000) {
  u.1 <- runif(n = 1, min = my.trun.tri.dist$trun.tri.lower, max = my.trun.tri.dist$trun.tri.upper)
  u.2 <- runif(n = 1, min = 0, max = my.trun.tri.dist$trun.tri.mode)
  if (u.2 <= my.trun.tri.dist$pdf(u.1)) {
    success.count <- success.count + 1
    vec[success.count] <- u.1
  }
}

hist(vec)
