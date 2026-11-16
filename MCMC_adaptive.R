MCMC_adaptive_1 <- function(obj, 
                               allocate_parameters,
                               goveqs,
                               goveqs_scaleup,
                               theta0, 
                               niter, 
                               mul_sigma, 
                               bounds, 
                               prm, 
                               ref, 
                               thetai, 
                               cov0,
                               lhd_list,
                               lhd_switch,
                               lhd_weight,
                               data_notif_bounds,
                               old_inc_year,
                               output_names,
                               hiv_model = FALSE) {
  
  # result$samples - parameter samples
  # result$values - log likelihood at each iteration
  # result$sim_outputs - incidence, mortality, notifications, prevalence, accept, etc.
  # result$accept_rate - overall acceptance rate
  # currently not saving lhd and sfin in sim_outputs but can add to sim_sto
  
  failed_result <- function() {c(list(lhd = -1e10),
                                 setNames(as.list(rep(NA_real_, length(output_names))),
                                          output_names))}
  
  d <- length(theta0)  # number of parameters
  b <- 0.05            # mixing weight for initial vs adaptive proposal
  # sd <- mul_sigma * (2.4^2) / d  # adaptive scaling factor
  sd <- 1 # use this generally
  
  # Initial covariance (if none provided): diagonal matrix
  if (is.null(cov0)) {
    cov0 <- 1e-2 * diag(d)
  }
  
  # Initialize storage
  xsto <- matrix(0, nrow = d, ncol = niter)  # sampled parameters
  outsto <- numeric(niter)                  # Vector to store objective function values (log-likelihood)
  
  # Set initial values
  xsto[, 1] <- theta0
  xsel <- theta0
  
  # Store initial values for objective function
  resX <- obj(theta = as.vector(theta0), 
              allocate_parameters,
              goveqs,
              goveqs_scaleup,
              bounds,
              prm, 
              ref, 
              thetai,
              lhd_list,
              lhd_switch,
              lhd_weight,
              data_notif_bounds,
              old_inc_year,
              output_names)
  FX <- resX[["lhd"]]
  outsto[1] <- FX
  acc <- 0 # accept counter
  
  first_row <- c(unlist(resX[output_names]), accept = 1)
  sim_sto <- matrix(NA, nrow = niter, ncol = length(first_row),
                    dimnames = list(NULL, names(first_row)))
  sim_sto[1, ] <- first_row
  
  
  # Start the MCMC loop
  for (t in 2:niter) {
    # Print to see progress of sampler
    if(t %% 1000 == 0) 
      print(paste("Iteration: ", t))  
    # Current parameter value
    X <- xsto[ , t-1]
    
    # STEP 1: Generate initial proposal Y0 (always needed)
    Y0 <- MASS::mvrnorm(1, mu = X, Sigma = 0.1^2 * cov0 * mul_sigma / d)
    
    if (t < 101) {
      # STEP 2a: Before iteration 101, use Y0 directly, and ensure non-negative
      Y <- pmax(Y0, 0)
    } else {
      # STEP 2b: From iteration 101 onward, use adaptive proposal
      # Estimate covariance matrix from last 100 samples
      ind0 <- max(1, t - 100)
      ind1 <- t - 1
      covmat <- cov(t(xsto[, ind0:ind1]))
      covmat <- (covmat + t(covmat)) / 2  # enforce symmetry
      
      # Adaptive proposal from N(X, sd * covmat)
      Y_adaptive <- tryCatch({
        MASS::mvrnorm(1, mu = X, Sigma = sd * covmat)
      }, error = function(e) {
        # fallback if proposal fails
        MASS::mvrnorm(1, mu = X, Sigma = sd * cov0)
      })
      
      # STEP 3: Mixture proposal, ensure non-negativity
      Y <- pmax((1 - b) * Y_adaptive + b * Y0, 0)
    }
    
    # STEP 4: Evaluate log-likelihood for proposal
    # Select parameter depending on model type
    mort_value <- Y[thetai$rf_mort_TB]
    
    if (hiv_model) {
      mort_value <- mort_value[[1]]
    }
    
    # Check constraint
    if (any(mort_value < prm$r$muTx[2] / prm$r$mort_TB0)) {
      print(paste0(t, " rejecting this parameter: ", paste(mort_value, collapse = ", "), " < ", prm$r$muTx[2] / prm$r$mort_TB0))
      resY <- failed_result()
    } else {
      resY <- tryCatch({
        `R.utils`::withTimeout({
          obj(theta = Y, 
              allocate_parameters,
              goveqs,
              goveqs_scaleup,
              bounds, prm, ref, 
              thetai, 
              lhd_list,
              lhd_switch,
              lhd_weight,
              data_notif_bounds,
              old_inc_year)
        }, timeout = 3, onTimeout = "error")  # 3 sec timeout
      }, error = function(e){
        failed_result()
      })
    }
    FY <- resY[["lhd"]]
    
    # STEP 5: Accept/reject based on Metropolis-Hastings ratio
    delta <- FY - FX
    if (!is.finite(delta)) delta <- -1e10
    accept_prob <- exp(delta)
    if (!is.finite(accept_prob)) accept_prob <- 0
    rnd <- runif(1)
    
    if (rnd < accept_prob) {
      # accept
      xsel <- Y
      FX <- FY
      acc <- acc + 1
      sim_sto[t, ] <- c(unlist(resY[output_names]), accept = 1)
    } else {
      # reject
      xsel <- xsto[, t - 1]
      sim_sto[t, ] <- sim_sto[t - 1, ]
      sim_sto[t, "accept"] <- 0  # Mark as rejected
    }
    
    # STEP 6: Store results
    xsto[, t] <- xsel
    outsto[t] <- FX
  }
  
  accept_rate <- acc / niter
  return(list(samples = xsto, 
              values = outsto, 
              accept_rate = accept_rate,
              sim_outputs = sim_sto))
}
