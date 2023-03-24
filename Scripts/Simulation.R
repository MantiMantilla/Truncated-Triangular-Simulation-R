## Import relevant functions to generate triangular and truncated triangular distributions with cached properties (created as a list)
source(file.path("Scripts", "Generate-Triangular.R"))
source(file.path("Scripts", "Generate-Truncated-Triangular.R"))
test <- function(testcaseparameters, plotPDFs=FALSE, plotCDFs=FALSE, plotSIM=FALSE) {
  my.L <- testcaseparameters[1]
  my.a <- testcaseparameters[2]
  my.M <- testcaseparameters[3]
  my.b <- testcaseparameters[4]
  my.U <- testcaseparameters[5]
  print(paste("Test Case: ", "L=",my.L,", a=", my.a, ", M=",my.M,", b=",my.b,", U=", my.U, sep="", collapse=""))
  
  my.plotPDFs <- plotPDFs
  my.plotCDFs <- plotCDFs
  
  my.sim.n <- 20000
  my.sim.triFrom <- my.L
  my.sim.triTo <- my.U
  my.sim.trunFrom <- my.a - 0.0001
  my.sim.trunTo <- my.b + 0.0001
  my.sim.triFromCdf <- my.L
  my.sim.triToCdf <- my.U
  my.sim.trunFromCdf <- my.L
  my.sim.trunToCdf <- my.U
  my.sim.plot <- plotSIM
  
  ## Set a seed so results are replicable
  set.seed(10)
  
  ## Create a list that contains the major properties of a triangular distribution with L, U, and M
  my.tri.dist <- generate.triangular(L = my.L, U = my.U, M = my.M)
  
  ## Create a list that contains the major properties of a truncated triangular distribution between a and b, based on the previous triangular
  my.trun.tri.dist <-
    generate.truncated.triangular(a = my.a,
                                  b = my.b,
                                  orig.tri.dist = my.tri.dist)
  
  print("Triangular Distribution:")
  print(paste("  Analytical Mean:     ", my.tri.dist$tri.mean))
  print(paste("  Analytical Median:   ", my.tri.dist$tri.median))
  print(paste("  Analytical Variance: ", my.tri.dist$tri.var))
  print("Truncated Triangular Distribution:")
  print(paste("  Analytical Mean:     ", my.trun.tri.dist$trun.tri.mean))
  print(paste("  Analytical Median:   ", my.trun.tri.dist$trun.tri.median))
  print(paste("  Analytical Variance: ", my.trun.tri.dist$trun.tri.var))

  ## Plot and describe the PDFs
  if (my.plotPDFs) {
    plot.function(
        function(x) {
          sapply(x, FUN = my.trun.tri.dist$pdf)
        },
        from = my.sim.triFrom,
        to = my.sim.triTo,
        n = my.sim.n,
        ylab = "f( x | a < x < b )",
        xlab = "x",
        main = "Truncated Triangular Distribution - PDF"
    )
    
    plot.function(
      function(x) {
        sapply(x, FUN = my.tri.dist$pdf)
      },
      from = my.sim.triFrom,
      to = my.sim.triTo,
      n = my.sim.n,
      ylab = "f(x)",
      xlab = "x",
      main = "Triangular Distribution - PDF",
      add = TRUE,
      col = "red"
    )
  }
  ## Plot and describe the CDFs
  if (my.plotCDFs) {
    plot.function(
      function(x) {
        sapply(x, FUN = my.trun.tri.dist$cdf)
      },
      from = my.sim.trunFromCdf,
      to = my.sim.trunToCdf,
      n = my.sim.n,
      ylab = "F( x | a < x < b )",
      xlab = "x",
      main = "Truncated Triangular Distribution - CDF"
    )
    plot.function(
      function(x) {
        sapply(x, FUN = my.tri.dist$cdf)
      },
      from = my.sim.triFromCdf,
      to = my.sim.triToCdf,
      n = my.sim.n,
      ylab = "F(x)",
      xlab = "x",
      main = "Triangular Distribution - CDF",
      add = TRUE,
      col = "red"
    )
    plot.function(
      function(x) { sapply(x, FUN = function(y){0.5})},
      from = my.sim.triFromCdf,
      to = my.sim.triToCdf,
      n = my.sim.n,
      add = TRUE,
      col = "blue"
    )
  }
  
  ## Simulate n instances of a random variable with f(x) = truncated triangular by trial and error.
  vec <- rep(0, my.sim.n)
  success.count <- 0
  success.truntrimax <- if (my.trun.tri.dist$trun.tri.mode <= my.trun.tri.dist$trun.tri.lower) {
    # at lower it will be 0, go slightly above that
    (my.trun.tri.dist$pdf(my.trun.tri.dist$trun.tri.lower+0.00000001))
  } else if (my.trun.tri.dist$trun.tri.mode > my.trun.tri.dist$trun.tri.upper) {
    (my.trun.tri.dist$pdf(my.trun.tri.dist$trun.tri.upper))
  } else {
    (my.trun.tri.dist$pdf(my.trun.tri.dist$trun.tri.mode))
  }
  while (success.count < my.sim.n) {
    u.1 <-
      runif(
        n = 1,
        min = my.trun.tri.dist$trun.tri.lower,
        max = my.trun.tri.dist$trun.tri.upper
      )
    u.2 <-
      runif(
        n = 1,
        min = 0,
        max = success.truntrimax
      )
    if (u.2 <= my.trun.tri.dist$pdf(u.1)) {
      success.count <- success.count + 1
      vec[success.count] <- u.1
    }
  }
  
  ## Simulate n instances of a random variable with f(x) = truncated triangular by inverse cdf method (faster than by trial and error)
  rand.unif <- runif(my.sim.n, min = 0, max = 1)
  vec2 <- sapply(rand.unif, my.trun.tri.dist$inverse.cdf)
  
  print(paste("Simulation 1 Mean:     ", mean(vec)))
  print(paste("Simulation 1 Variance: ", var(vec)))
  print(paste("Simulation 2 Mean:     ", mean(vec2)))
  print(paste("Simulation 2 Variance: ", var(vec2)))

  if (my.trun.tri.dist$trun.tri.mean < my.a
      || my.trun.tri.dist$trun.tri.mean > my.b) {
    print(paste("FAILED test 'lower <= mean <= upper': ",
                my.a,
                "<=", my.trun.tri.dist$trun.tri.mean,
                "<=", my.b))
  }
  if (my.trun.tri.dist$trun.tri.var < 0) {
    print(paste("FAILED test 'variance is positive': ",
                my.trun.tri.dist$trun.tri.var))
  }
  
  ## following tests are not statistically significant, but are a quick check to verify that
  ## it is at least in the right ballpark for our test cases
  if (abs((mean(vec2) - my.trun.tri.dist$trun.tri.mean)/my.trun.tri.dist$trun.tri.mean) > 0.01) {
    print(paste("FAILED test 'simulated mean should be close to analytical mean': ",
                abs((mean(vec2) - my.trun.tri.dist$trun.tri.mean)/my.trun.tri.dist$trun.tri.mean)))
  }
  if (abs((var(vec2) - my.trun.tri.dist$trun.tri.var) / my.trun.tri.dist$trun.tri.mode) > 0.01) {
    print(paste("FAILED test 'simulated variance should be close to analytical variance': ",
                abs((var(vec2) - my.trun.tri.dist$trun.tri.var) / my.trun.tri.dist$trun.tri.mode)))
  }

  if (my.sim.plot) {
    ## Make the truncated PDF compatible with vectors for plotting
    trun.vector.pdf <- function(x) {
      sapply(x, FUN = my.trun.tri.dist$pdf)
    }
    
    ## Plot and describe the first simulation
    hist(vec,
         xlab = "x",
         breaks = 20,
         freq = FALSE)
    
    curve(
      trun.vector.pdf,
      from = my.sim.trunFrom,
      to = my.sim.trunTo,
      n = my.sim.n,
      col = "red",
      add = TRUE
    )
  
    ## Plot and describe the second simulation
    hist(vec2,
         xlab = "x",
         breaks = 20,
         freq = FALSE)
    curve(
      trun.vector.pdf,
      from = my.sim.trunFrom,
      to = my.sim.trunTo,
      n = my.sim.n,
      col = "red",
      add = TRUE
    )
  }
}

testcases <- list(
  #c(L,   a,    M,         b,    U)
  #c(0.5, 0.75, 0.50,      1.25, 2),
  #c(0.5, 0.75, 0.55,      1.25, 2),
  #c(0.5, 0.75, 0.749,     1.25, 2),
  #c(0.5, 0.75, 0.75,      1.25, 2),
  #c(0.5, 0.75, 0.7500001, 1.25, 2),
  #c(0.5, 0.75, 1.00,      1.25, 2),
  #c(0.5, 0.75, 1.249999,  1.25, 2),
  #c(0.5, 0.75, 1.25,      1.25, 2),
  #c(0.5, 0.75, 1.251,     1.25, 2),
  #c(0.5, 0.75, 1.50,      1.25, 2),
  #c(0.5, 0.75, 2.00,      1.25, 2),
  c(2,   3,    2.5,       8,    9),
  c(2,   3,    8.5,       8,    9),
  c(2,   3,    7,         8,    9)
)

for (testcase in testcases) {
  test(unlist(testcase))
}
