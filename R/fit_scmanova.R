###############################################################################
#### functions needed to run a regularized MANOVA test on a simulated scenario
#### with semicontinuous high-dimensional data
###############################################################################


###### ESTIMATION ######
## function that returns the estimated parameters of the model 
## and run the regularized MANOVA test
#' @param x , n x p matrix with the data
#' @param n , vector with the number of observations for each group
#' @param lambda ,  minimum and maximum values between we search for the optimal lambda
#' @param lambda0 , minimum and maximum values between we search for the optimal lambda0
#' @param B , number of permutations used in the permutation test
#' @param penalty , penalty used in the Information criteria
#' @param ident , use or not the identity function in the regularization
# It outputs the estimated parameters of the test and the p-value of the permutation test

estimation <- function(
    x,
    n,
    lambda = c(0, 100),
    lambda0 = c(0, 100),
    B = 1000,
    penalty = function(n, p) log(n) + 0.5 * log(p),
    ident = TRUE,
    missing.model = NULL
) {
  # All thesis replication scripts are intentionally serial.
  # Build scMANOVA call.
  # If missing.model is NULL, we do NOT pass it.
  # Then scMANOVA uses the package default, i.e. the original paper model.
  sc_args <- list(
    x = x,
    n = n,
    lambda = lambda,
    lambda0 = lambda0,
    tol = 1e-8,
    p.value.perm = TRUE,
    fixed.lambda = FALSE,
    parallel = "no",
    ncpus = 1L,
    B = B,
    penalty = penalty,
    ident = ident
  )
  
  if (!is.null(missing.model)) {
    sc_args$missing.model <- missing.model
  }
  
  # Estimate the parameters and run the test.
  res <- try(do.call(semicontMANOVA::scMANOVA, sc_args))
  
  if (class(res)[1] == "try-error") {
    print("Try-error")
    rep(NA, 12)
  } else {
    sigmaRidge <- res$sigmaRidge
    sigma0Ridge <- res$sigma0Ridge
    logLik <- res$logLik
    logLik0 <- res$logLik0
    lambda <- res$lambda
    lambda0 <- res$lambda0
    df <- res$df
    df0 <- res$df0
    aic <- res$aic
    aic0 <- res$aic0
    statistic <- res$statistic
    p.valuePERM <- res$p.value
    
    c(
      logLik,
      logLik0,
      lambda,
      lambda0,
      df,
      df0,
      aic,
      aic0,
      statistic,
      p.valuePERM,
      dim(sigmaRidge)[1],
      dim(sigma0Ridge)[1]
    )
  }
}
########################


