get_distribution_fns <- function(prctiles, distribution, visualizing = TRUE) {
  
  # Ensure percentiles are sorted
  prctiles <- sort(prctiles)
  
  # Initial estimates of mean and variance
  mn <- prctiles[2]
  var <- (prctiles[3] - prctiles[1])^2 / 4
  
  if (distribution == "lognormal") {
    # Define function for optimization (matching percentiles)
    # eg. plnorm(2, meanlog=0, sdlog=1) is the prob that a lognormally distributed 
    # random variable is <= 2
    obj_fun <- function(params) {
      mu <- params[1]
      sigma <- params[2]
      cdf_vals <- plnorm(prctiles, meanlog = mu, sdlog = sigma)
      sum((cdf_vals / c(0.025, 0.50, 0.975) - 1)^2)
    }
    
    # Initial guess based on method of moments (this is on the log scale, so log(mn) = mu_init)
    sigma_init <- sqrt(log(var / mn^2 + 1))
    mu_init <- log(mn) - sigma_init^2 / 2
    
    # Optimization to find best parameters
    opt_result <- optim(c(mu_init, sigma_init), obj_fun, method = "Nelder-Mead")
    
    if (opt_result$convergence != 0) {
      stop("Optimization did not converge")
    }
    
    mu <- opt_result$par[1]
    sigma <- opt_result$par[2]
    
    # Log-likelihood function
    log_likelihood <- function(x) {
      -(log(x) - mu)^2 / (2 * sigma^2) - log(x * sigma * sqrt(2 * pi))
    }
    
    # Distribution parameters
    distrib_params <- c(mu, sigma)
    
  } else if (distribution == "beta") {
    # Objective function for beta optimization
    obj_fun <- function(params) {
      a <- params[1]
      b <- params[2]
      cdf_vals <- pbeta(prctiles, shape1 = a, shape2 = b)
      sum((cdf_vals / c(0.025, 0.50, 0.975) - 1)^2)
    }
    
    # Initial guess using method of moments
    # You can prove this by replacing the value of mn and var from beta dist.
    tmp <- mn * (1 - mn) / var - 1
    a_init <- tmp * mn
    b_init <- tmp * (1 - mn)
    
    # Optimization with non-negative constraints
    opt_result <- optim(c(a_init, b_init), obj_fun, method = "Nelder-Mead")
    
    if (opt_result$convergence != 0) {
      stop("Optimization did not converge")
    }
    
    a <- opt_result$par[1]
    b <- opt_result$par[2]
    
    # Log-likelihood function
    log_likelihood <- function(x) {
      (a - 1) * log(x) + (b - 1) * log(1 - x) - lbeta(a, b)
    }
    
    # Distribution parameters
    distrib_params <- c(a, b)
  } else {
    stop("Invalid distribution. Choose 'lognormal' or 'beta'.")
  }
  
  # Visualization
  if (visualizing) {
    x_vals <- seq(prctiles[1] * 0.8, prctiles[3] * 1.2, length.out = 100)
    if (distribution == "lognormal") {
      y_vals <- dlnorm(x_vals, meanlog = mu, sdlog = sigma)
    } else {
      y_vals <- dbeta(x_vals, shape1 = a, shape2 = b)
    }
    
    df <- data.frame(x = x_vals, y = y_vals)
    plot_obj <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
      ggplot2::geom_line(color = "blue") +
      ggplot2::geom_vline(xintercept = prctiles, linetype = "dashed", color = "red") +
      ggplot2::labs(title = paste("Estimated", distribution, "PDF"),
           x = "x", y = "Density") +
      ggplot2::theme_minimal()
    
    print(plot_obj)
  }
  
  return(list(log_likelihood = log_likelihood, parameters = distrib_params))
}

# # Example Usage: (matches Nim's code)
# likelihood_f <- get_distribution_fns(c(164, 195, 228), "lognormal", visualizing = TRUE)[[1]]
# likelihood <- likelihood_f(179)
