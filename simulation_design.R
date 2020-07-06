expit <- function(x) {
  exp(x) / (1 + exp(x))
}

data_gen_censor = function(n=5000, p=10,  PH=TRUE, censor = "20%", setting = 1, overlap = "strong") {
  
  # n <- 50000
  # p <- 10
  X <- matrix(rnorm(p * n, sd = 0.35), nrow = n, ncol = p / 2)
  x1 <- X[, 1]
  x2 <- X[, 2]
  x3 <- X[, 3]
  x4 <- X[, 4]
  x5 <- X[, 5]
 
  x6 <- sample(0:1, n, replace = TRUE, prob = c(0.5, 0.5))
  x7 <- sample(0:1, n, replace = TRUE, prob = c(0.5, 0.5))
  x8 <- sample(0:1, n, replace = TRUE, prob = c(0.5, 0.5))
  x9 <- sample(0:1, n, replace = TRUE, prob = c(0.5, 0.5))
  x10 <- sample(0:1, n, replace = TRUE, prob = c(0.5, 0.5))
  
  if (overlap == "strong") {
    p1 <- expit(
      -0.1 * x1 - .9 * x2 - .3 * x3  - .1 * x5 - 0.1 * x6 -0.2*x7 - 0.4 * x9 + .5 * x10
    )
  }
 
 
  if (overlap == "medium") {
    p1 <- expit(
      0.3-0.1 * 2.5 * x1 - .9 * 2.5 * x3 - .3 * 2.5 * x3  - .1 * 2.5 * x5 - 0.1 * 2.5 * x6 -0.3 * 2.5* x7 - 0.4 * 2.5 * x9 + .5 * 2.5 * x10
    )
  }
  
  if (overlap == "weak") {
    p1 <- expit(
      0.7 - 0.1 * 5 * x1 - .9 * 5 * x3 - .3 * 5 * x3  - .1 * 5 * x5 - 0.1 * 5 * x6 -0.3 * 5* x7 - 0.4 * 5 * x9 + .5 * 5 * x10
    )
  }
  
  
  p2 <- 1 - p1
  p1_category <- cut(p1, breaks=c(quantile(p1, probs = seq(0, 1, by = 0.02))), include.lowest = T, labels = 1:50) 
  W = NULL
  for (i in 1:n) {
    W[i] <- sample(c(0, 1),
                   size = 1,
                   replace = TRUE,
                   prob = c(p1[i], p2[i]))
  }
  table(W)

  if (PH == FALSE) {
    eta_1 <- exp(0.7 + 1.8 * x3 + 0.8 * x7)
    eta_2 <- exp(0.9 - 0.5 * x1 + 0.5 * x2)

    summary(eta_1)
  } 
  if (PH == TRUE) {
    eta_1 <- 2
    eta_2 <- 2
  }
  
  if (setting == 1) {
    LP1 <-
      .2 - 0.5 * x1 - 0.8 * x3 - 1.8 * x5 - 0.9 * x6 - 0.1 * x7
    LP2 <-
      -.2 + 0.1 * expit(x1) + 0.1 * sin(x3) - 0.1 * x5 ^ 2   - 0.3 * x6 - 0.2 * x7
  }
  
  if (setting == 2) {
    LP1 <-
      -.1+ 0.1 * x1 ^ 2 - 0.2 * sin(x3) + 0.2 * expit(x5) + 0.2 * x6 -0.3 * x7
    LP2 <-
      -.2 + 0.1 * expit(x1) + 0.1 * sin(x3) - 0.1 * x5 ^ 2   - 0.3 * x6 - 0.2 * x7
  }
  
  if (setting == 3) {
    LP1 <-
      -.1+ 0.1 * x1 ^ 2 - 0.2 * sin(x3) + 0.2 * expit(x5) + 0.2 * x6 -0.3 * x7
    LP2 <- .5 - 0.8* expit(x2) + 0.1 * sin(x3) - 0.1 * x4 ^ 2 + 0.2 * x4 - 0.1 * x5 ^ 2 - 0.3 * x6 
  }
  # summary(LP1)
  # summary(LP2)
  
  #independent censoring
  if (censor == "20%"){
    C <- rexp(n, rate = 0.007)
  }
  
  if (censor == "60%"){
    C <- rexp(n, rate = 0.03)
  }

  #generate U ~ unif(0,1)
  U = runif(n, 0, 1)
  #scale parameter lambda>0, exp(Linear predictor)
  lambda1 <- 1200
  lambda2 <- 2000
  
  #potential survival times
  T1 <- (lambda1 * (-log(U)) / exp(LP1)) ^ (1 / eta_1)
  summary(T1)
  # mean(T1>60)
  T2 <- (lambda2 * (-log(U)) / exp(LP2)) ^ (1 / eta_2)
  summary(T2)
  # mean(T2>60)
  mean(T1>72)
  mean(T2>72)
  # observed outcomes
  T <- cbind(T1, T2)
  TW <- cbind(T, W)
  Tobs <- apply(TW, 1, function(x)
    x[1:2][x[3] + 1]) #observed when trt is received
  # Tobs[1:10]
  # TW[1:10,]
  summary(Tobs)

  
  Tobs_C <- pmin(Tobs, C)
  summary(Tobs_C)
  
  #censoring rate
  censor_rate <- sum(Tobs > C) / n
  censor_rate

  #censoring indicator
  delta <- ifelse(Tobs > C, 0, 1)
  summary(Tobs_C)
  
  return(
    list(
      LP1 = LP1,
      LP2 = LP2,
      lambda1 = lambda1,
      lambda2 = lambda2,
      eta_1 = eta_1,
      eta_2 = eta_2,
      p1_category =  p1_category,
      X = cbind(X, x6, x7, x8, x9, x10, W, Time = Tobs_C, Event =delta),
      T1 = T1,
      T2 = T2
    )
  )
  
}

