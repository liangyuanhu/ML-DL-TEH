#============================================================================#
# The codes are mostly borrowed from AFTrees:: SurvivalProb                  #
#============================================================================#
AFTrees_SurvivalProb <- function (object, time.points = NULL, xind.train = NULL, xind.test = NULL, 
          credible = FALSE, test.only = FALSE, train.only = FALSE) 
{ 
  # ntrain <- dim(object$m.train)[2]
  if (sum(time.points <= 0)) {
    stop("All time points must be positive")
  }
  if (test.only & is.null(object$m.test)) {
    stop("Must have fitted for test set when using test.only=TRUE")
  }
  if (is.null(time.points)) {
    log.times <- log(object$times)
    qq <- c(min(log.times), quantile(log.times, probs = seq(0.025, 
                                                            0.975, length.out = 30)), max(log.times))
    time.points <- exp(qq)
  }
  if (is.null(object$m.test)) {
    train.only = TRUE
  }
  nsubjects <- ncol(object$m.train)
  nsamples <- nrow(object$locations)
  ngrid <- length(time.points)
  log.time.points <- log(time.points)
  if (test.only) {
    ntest <- ncol(object$m.test)
    if (is.null(xind.test)) {
      xind.test <- 1:ntest
    }
    if (length(time.points) > 1) {
      SS.test <- matrix(0, nrow = ntest, ncol = ngrid)
      for (i in xind.test) {
        for (k in 1:ngrid) {
          Amat <- (log.time.points[k] - object$locations - 
                     object$m.test[, i])/object$sigma
          SS.test[i, k] <- sum(pnorm(Amat, lower.tail = FALSE) * 
                                 object$mix.prop)/nsamples
        }
      }
      SS.test.mean <- colMeans(SS.test)
    }
    else {
      SS.test <- rep(0, ntest)
      for (i in xind.test) {
        Amat <- (log.time.points - object$locations - 
                   object$m.test[, i])/object$sigma
        SS.test[i] <- sum(pnorm(Amat, lower.tail = FALSE) * 
                            object$mix.prop)/nsamples
      }
    }
    SS.train <- SS.train.mean <- NULL
  }
  else if (train.only) {
    if (is.null(xind.train)) {
      ntrain <- dim(object$m.train)[2]
      xind.train <- 1:ntrain
    }
    if (length(time.points) > 1) {
      SS.train <- matrix(0, nrow = nsubjects, ncol = ngrid)
      for (i in xind.train) {
        for (k in 1:ngrid) {
          Amat <- (log.time.points[k] - object$locations - 
                     object$m.train[, i])/object$sigma
          SS.train[i, k] <- sum(pnorm(Amat, lower.tail = FALSE) * 
                                  object$mix.prop)/nsamples
        }
      }
      SS.train.mean <- colMeans(SS.train)
    }
    else {
      SS.train <- rep(0, nsubjects)
      for (i in xind.train) {
        Amat <- (log.time.points - object$locations - 
                   object$m.train[, i])/object$sigma
        SS.train[i] <- sum(pnorm(Amat, lower.tail = FALSE) * 
                             object$mix.prop)/nsamples
      }
      SS.train.mean <- mean(SS.train)
    }
    SS.test <- SS.test.mean <- NULL
  }
  else if (!train.only & !test.only) {
    ntest <- ncol(object$m.test)
    if (is.null(xind.train)) {
      xind.train <- 1:ntrain
    }
    if (is.null(xind.test)) {
      xind.test <- 1:ntest
    }
    if (length(time.points) > 1) {
      SS.test <- matrix(0, nrow = ntest, ncol = ngrid)
      SS.train <- matrix(0, nrow = ntrain, ncol = ngrid)
      for (i in xind.test) {
        for (k in 1:ngrid) {
          Amat <- (log.time.points[k] - object$locations - 
                     object$m.test[, i])/object$sigma
          SS.test[i, k] <- sum(pnorm(Amat, lower.tail = FALSE) * 
                                 object$mix.prop)/nsamples
        }
      }
      for (i in xind.train) {
        for (k in 1:ngrid) {
          Amat <- (log.time.points[k] - object$locations - 
                     object$m.train[, i])/object$sigma
          SS.train[i, k] <- sum(pnorm(Amat, lower.tail = FALSE) * 
                                  object$mix.prop)/nsamples
        }
      }
      SS.test.mean <- colMeans(SS.test)
      SS.train.mean <- colMeans(SS.train)
    }
    else {
      SS.test <- rep(0, ntest)
      for (i in xind.test) {
        Amat <- (log.time.points - object$locations - 
                   object$m.test[, i])/object$sigma
        SS.test[i] <- sum(pnorm(Amat, lower.tail = FALSE) * 
                            object$mix.prop)/nsamples
      }
      SS.train <- rep(0, nsubjects)
      for (i in xind.train) {
        Amat <- (log.time.points - object$locations - 
                   object$m.train[, i])/object$sigma
        SS.train[i] <- sum(pnorm(Amat, lower.tail = FALSE) * 
                             object$mix.prop)/nsamples
      }
      SS.test.mean <- mean(SS.test)
      SS.train.mean <- mean(SS.train)
    }
  }
  ans <- list()
  class(ans) <- "aftsurvcurves"
  ans$Surv.train <- SS.train
  ans$Surv.test <- SS.test
  ans$Surv.train.mean <- SS.train.mean
  ans$Surv.test.mean <- SS.test.mean
  ans$time.points <- time.points
  return(ans)
}
