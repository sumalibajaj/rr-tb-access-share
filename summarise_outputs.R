make_agg_global <- function(draw_list, countries){
  
  # Year- and country-specific population from WPP 2024
  # Note: WPP population is in thousands, so multiply by 1000
  pop_tbl <- bind_rows(read_excel("data/WPP2024_GEN_F01_DEMOGRAPHIC_INDICATORS_COMPACT.xlsx",
                                  sheet = "Estimates", skip = 16),
                       read_excel("data/WPP2024_GEN_F01_DEMOGRAPHIC_INDICATORS_COMPACT.xlsx",
                                  sheet = "Medium variant", skip = 16)) %>%
    dplyr::filter(`Type` == "Country/Area") %>%
    dplyr::rename(country = `Region, subregion, country or area *`,
                  time = Year,
                  population_k = `Total Population, as of 1 July (thousands)`) %>%
    dplyr::distinct(country, time, .keep_all = TRUE) %>%
    dplyr::mutate(population = as.numeric(population_k)*1000,
                  time = as.numeric(time)) %>%
    dplyr::select(country, time, population)
  
  # Add country column and combine
  data <- purrr::map2_dfr(draw_list, countries, ~ .x %>% dplyr::mutate(country = .y)) %>%
    dplyr::mutate(time = as.numeric(time)) %>%
    dplyr::left_join(pop_tbl, by = c("country", "time")) %>%
    dplyr::select(country, time, population, inc_rr, tb_mort, tb_mort_ds, tb_mort_rr, scenario, iter)
  
  
  long_data <- data %>%
    tidyr::pivot_longer(
      cols = c(inc_rr, tb_mort, tb_mort_ds, tb_mort_rr),
      names_to = "indicator",
      values_to = "rate"
    ) %>%
    dplyr::mutate(
      count = rate * population / 1e5
    )
  
  # Note - Past scenarios: keep all years, Future scenarios: start in 2026 (2025 of future scenarios is garbage)
  global_draws <- long_data %>%
    dplyr::mutate(scenario = if_else(scenario == "Past no p_MDR_rec", "Counterfactual: no RR services", scenario),
                  scenario = if_else(scenario == "Past no MDR recognition", "Counterfactual: no RR services", scenario),
                  scenario = ifelse(scenario == "Improved public DST", "Improved DST", scenario),
                  scenario = ifelse(scenario == "Improved public DST + SL outcomes", "Improved DST + SL outcomes", scenario)) %>%
    dplyr::filter(scenario %in% c("Past actual", "Counterfactual: no RR services") | time > 2025) %>%
    dplyr::group_by(iter, time, scenario, indicator) %>%
    dplyr::summarise(
      total_count = sum(count, na.rm = TRUE),
      total_pop   = sum(population, na.rm = TRUE),
      rate        = total_count / total_pop * 1e5,
      .groups = "drop"
    )
  
  # To get summary data by countries 
  country_draws <- long_data %>%
    dplyr::mutate(scenario = if_else(scenario == "Past no p_MDR_rec", "Counterfactual: no RR services", scenario),
                  scenario = if_else(scenario == "Past no MDR recognition", "Counterfactual: no RR services", scenario),
                  scenario = ifelse(scenario == "Improved public DST", "Improved DST", scenario),
                  scenario = ifelse(scenario == "Improved public DST + SL outcomes", "Improved DST + SL outcomes", scenario)) %>%
    dplyr::filter(scenario %in% c("Past actual", "Counterfactual: no RR services") | time > 2025) %>%
    dplyr::group_by(iter, time, scenario, indicator, country) %>%
    dplyr::summarise(
      total_count = sum(count, na.rm = TRUE),
      total_pop   = sum(population, na.rm = TRUE),
      rate        = total_count / total_pop * 1e5,
      .groups = "drop"
    )
  
  
  # PAST COUNTERFACTUALS
  # Aggregate globally
  global_past_data <- global_draws %>%
    dplyr::filter(scenario %in% c("Past actual", "Counterfactual: no RR services"),
                  indicator %in% c("inc_rr", "tb_mort_rr"),
                  time > 1979) %>%
    dplyr::mutate(time = time - 1) %>%
    dplyr::group_by(time, scenario, indicator) %>%
    dplyr::summarise(
      median = median(rate, na.rm = TRUE),
      lower  = quantile(rate, 0.025, na.rm = TRUE),
      upper  = quantile(rate, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  global_past_plot <- ggplot(global_past_data, aes(time, median, colour = scenario, fill = scenario)) +
    geom_line() +
    expand_limits(y = 0) +
    labs(y = "Rate per 100k population", x = "Year",
         fill = "Scenario", color = "Scenario") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
    facet_wrap(~indicator, scales = "free_y",
               labeller = labeller(indicator = c(inc_rr = "RR-TB incidence",
                                                 tb_mort_rr = "RR-TB mortality"))) +
    scale_color_manual(
      values = c("Past actual"   = "#DDA448",
                 "Counterfactual: no RR services" = "#ef476f")) +
    scale_fill_manual(values = c("Past actual"   = "#DDA448",
                                 "Counterfactual: no RR services" = "#ef476f")) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text.x = element_text(angle = 25, hjust = 0.7))
  
  # FUTURE INTERVENTIONS
  global_future_data <- global_draws %>%
    dplyr::filter(scenario %in% c("Baseline", "Improved DST", "Improved DST + SL outcomes"),
                  indicator %in% c("inc_rr", "tb_mort_rr"),
                  time > 2025) %>%
    dplyr::mutate(time = time - 1) %>%
    dplyr::group_by(time, scenario, indicator) %>%
    dplyr::summarise(
      median = median(rate, na.rm = TRUE),
      lower  = quantile(rate, 0.025, na.rm = TRUE),
      upper  = quantile(rate, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  global_future_plot <- ggplot(global_future_data, aes(time, median, colour = scenario, fill = scenario)) +
    geom_line() +
    expand_limits(y = 0) +
    labs(y = "Rate per 100k population", x = "Year",
         fill = "Scenario", color = "Scenario") +
    scale_x_continuous(breaks = seq(min(global_future_data$time), max(global_future_data$time), by = 2)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
    facet_wrap(~indicator, scales = "free_y",
               labeller = labeller(indicator = c(inc_rr = "RR-TB incidence",
                                                 tb_mort_rr = "RR-TB mortality"))) +
    scale_color_manual(values = c("Baseline" = "#DDA448",
                                  "Improved DST" = "#785589",
                                  "Improved DST + SL outcomes" = "#14591D")) +
    scale_fill_manual(values = c("Baseline" = "#DDA448",
                                 "Improved DST" = "#785589",
                                 "Improved DST + SL outcomes" = "#14591D")) +
    
    theme_minimal() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text.x = element_text(angle = 25, hjust = 0.7))
  
  
  make_summary <- function(global_draws, ref_scenario, comp_scenario, period_label) {
    
    tmp <- global_draws %>%
      dplyr::filter(
        scenario %in% c(ref_scenario, comp_scenario),
        indicator %in% c("inc_rr", "tb_mort_rr")
      ) %>%
      dplyr::group_by(iter, indicator, scenario) %>%
      dplyr::summarise(
        cumulative = sum(total_count, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(names_from = scenario, values_from = cumulative) %>%
      dplyr::mutate(
        cumulative_reference = .data[[ref_scenario]],
        cumulative_comparison = .data[[comp_scenario]],
        cumulative_averted = if (period_label == "Past") {
          cumulative_comparison - cumulative_reference
        } else {
          cumulative_reference - cumulative_comparison
        },
        percent_averted = if (period_label == "Past") {
          100 * cumulative_averted / pmax(cumulative_comparison, 1e-12)
        } else {
          100 * cumulative_averted / pmax(cumulative_reference, 1e-12)
        },
        times_higher = if (period_label == "Past") {
          cumulative_comparison / pmax(cumulative_reference, 1e-12)
        } else {
          cumulative_reference / pmax(cumulative_comparison, 1e-12)
        }
      )
    
    tmp %>%
      dplyr::group_by(indicator) %>%
      dplyr::summarise(
        ref_med = median(cumulative_reference, na.rm = TRUE),
        ref_low = quantile(cumulative_reference, 0.025, na.rm = TRUE),
        ref_up  = quantile(cumulative_reference, 0.975, na.rm = TRUE),
        
        comp_med = median(cumulative_comparison, na.rm = TRUE),
        comp_low = quantile(cumulative_comparison, 0.025, na.rm = TRUE),
        comp_up  = quantile(cumulative_comparison, 0.975, na.rm = TRUE),
        
        av_med = median(cumulative_averted, na.rm = TRUE),
        av_low = quantile(cumulative_averted, 0.025, na.rm = TRUE),
        av_up  = quantile(cumulative_averted, 0.975, na.rm = TRUE),
        
        pct_med = median(percent_averted, na.rm = TRUE),
        pct_low = quantile(percent_averted, 0.025, na.rm = TRUE),
        pct_up  = quantile(percent_averted, 0.975, na.rm = TRUE),
        
        th_med = median(times_higher, na.rm = TRUE),
        th_low = quantile(times_higher, 0.025, na.rm = TRUE),
        th_up  = quantile(times_higher, 0.975, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::transmute(
        period = period_label,
        indicator,
        reference_scenario = ref_scenario,
        comparison_scenario = comp_scenario,
        
        cumulative_reference = ref_med,
        cumulative_reference_lower = ref_low,
        cumulative_reference_upper = ref_up,
        
        cumulative_comparison = comp_med,
        cumulative_comparison_lower = comp_low,
        cumulative_comparison_upper = comp_up,
        
        cumulative_averted = av_med,
        cumulative_averted_lower = av_low,
        cumulative_averted_upper = av_up,
        
        percent_averted = pct_med,
        percent_averted_lower = pct_low,
        percent_averted_upper = pct_up,
        
        times_higher = th_med,
        times_higher_lower = th_low,
        times_higher_upper = th_up
      )
  }
  
  global_past_summary_table <- make_summary(
    global_draws,
    ref_scenario = "Past actual",
    comp_scenario = "Counterfactual: no RR services",
    period_label = "Past"
  )
  
  global_future_summary_table <- dplyr::bind_rows(
    make_summary(global_draws, "Baseline", "Improved DST", "Future"),
    make_summary(global_draws, "Baseline", "Improved DST + SL outcomes", "Future")
  )
  
  global_summary_table <- rbind(global_past_summary_table, global_future_summary_table)
  global_summary_table$country = "global"
  
  # country specific summary table
  country_past_summary_table <- dplyr::bind_rows(
    lapply(unique(country_draws$country), function(ctry) {
      make_summary(
        global_draws = country_draws %>% dplyr::filter(country == ctry),
        ref_scenario = "Past actual",
        comp_scenario = "Counterfactual: no RR services",
        period_label = "Past"
      ) %>%
        dplyr::mutate(country = ctry)
    })
  )
  
  country_future_summary_table <- dplyr::bind_rows(
    lapply(unique(country_draws$country), function(ctry) {
      dplyr::bind_rows(
        make_summary(
          global_draws = country_draws %>% dplyr::filter(country == ctry),
          ref_scenario = "Baseline",
          comp_scenario = "Improved DST",
          period_label = "Future"
        ),
        make_summary(
          global_draws = country_draws %>% dplyr::filter(country == ctry),
          ref_scenario = "Baseline",
          comp_scenario = "Improved DST + SL outcomes",
          period_label = "Future"
        )
      ) %>%
        dplyr::mutate(country = ctry)
    })
  )
  
  country_summary_table <- dplyr::bind_rows(country_past_summary_table, country_future_summary_table)
  
  summary_table <- dplyr::bind_rows(global_summary_table, country_summary_table)
  
  return(list(global_past_data = global_past_data, global_past_plot = global_past_plot,
              global_future_data = global_future_data, global_future_plot = global_future_plot,
              summary_table = summary_table))
}


plot_mcmc_traces_1 <- function(mcmc_outputs, thetai, thin = 1, burnin = 0.5) {
  # Filter out NULL chains
  mcmc_outputs <- mcmc_outputs[!sapply(mcmc_outputs, is.null)]
  
  # If all chains are NULL, stop and return a message
  if (length(mcmc_outputs) == 0) {
    stop("All MCMC chains are NULL - nothing to plot.")
  }
  
  # Combine chains with chain ID
  all_chains <- lapply(seq_along(mcmc_outputs), function(i) {
    samples <- mcmc_outputs[[i]]$samples
    # Keep only every 'thin' iteration
    keep_iter <- seq(1, ncol(samples), by = thin)
    samples_thin <- samples[, keep_iter, drop = FALSE]
    samples_thin %>% t() %>% as.data.frame() %>%
      mutate(iter = keep_iter, chain = paste0("Chain_", i))
  }) %>% bind_rows()
  
  # Assign proper column names
  # Generate names
  col_names <- unlist(mapply(function(name, values) {
    if (length(values) > 1) {
      paste0(name, seq_along(values))
    } else {
      name
    }
  }, names(thetai), thetai))
  colnames(all_chains)[1:(ncol(all_chains)-2)] <- col_names
  
  # Convert to long format for faceting
  long_chains <- all_chains %>%
    pivot_longer(cols = -c(iter, chain), names_to = "parameter", values_to = "value") %>%
    group_by(chain, parameter) %>%
    filter(iter > burnin * max(iter))  # discard burn-in
  
  # Generate trace plot
  p <- ggplot(long_chains, aes(x = iter, y = value, color = chain)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ parameter, scales = "free_y") +
    theme_minimal() +
    labs(x = "Iteration", y = "Parameter Value",
         title = "MCMC parameter trace across chains")
  p
}

plot_mcmc_ll_1 <- function(mcmc_outputs, thin = 1, burnin = 0.5) {
  # Filter out NULL chains
  mcmc_outputs <- mcmc_outputs[!sapply(mcmc_outputs, is.null)]
  
  # If all chains are NULL, stop and return a message
  if (length(mcmc_outputs) == 0) {
    stop("All MCMC chains are NULL - nothing to plot.")
  }
  
  # Combine chains with chain ID
  all_chains <- lapply(seq_along(mcmc_outputs), function(i) {
    samples <- mcmc_outputs[[i]]$values %>% as.data.frame() %>%t()
    # Keep only every 'thin' iteration
    keep_iter <- seq(1, ncol(samples), by = thin)
    samples_thin <- samples[, keep_iter, drop = FALSE]
    samples_thin %>% t() %>% as.data.frame() %>%
      mutate(iter = keep_iter, chain = paste0("Chain_", i))
  }) %>% bind_rows() %>%
    filter(iter > burnin * max(iter))  # discard burn-in
  
  # Assign proper column names
  # Generate names
  colnames(all_chains)[1] <- "log_likelihood"
  
  # Generate trace plot
  p <- ggplot(all_chains, aes(x = iter, y = log_likelihood, color = chain)) +
    geom_line(alpha = 0.7) +
    theme_minimal() +
    labs(x = "Iteration", y = "Log likelihood Value",
         title = "MCMC log likelihood trace across chains")
  p
}

summarise_posterior_with_rhat_1 <- function(mcmc_outputs, thetai, burnin_frac = 0.5, thin = 1){
  
  # Step 1: Process each chain individually (burn-in and thinning)
  chains_post <- lapply(mcmc_outputs, function(x) {
    chain <- t(x$samples)  # Transpose to iterations × parameters
    burnin <- floor(nrow(chain) * burnin_frac)
    chain_post <- chain[-seq_len(burnin), , drop = FALSE]
    chain_post <- chain_post[seq(1, nrow(chain_post), by = thin), , drop = FALSE]
    return(chain_post)  # rows = iterations, cols = parameters
  })
  
  # Step 2: Get dimensions
  m <- length(chains_post)       # number of chains
  n <- nrow(chains_post[[1]])    # iterations per chain
  d <- ncol(chains_post[[1]])    # number of parameters
  
  # Step 3: Name the parameters using thetai
  col_names <- unlist(mapply(function(name, values) {
    if (length(values) > 1) paste0(name, seq_along(values)) else name
  }, names(thetai), thetai))
  
  for (i in seq_along(chains_post)) {
    colnames(chains_post[[i]]) <- col_names
  }
  
  # Step 4: Compute Rhat for each parameter
  rhat_vals <- numeric(d)
  for (j in 1:d) {
    chains_j <- sapply(chains_post, function(mat) mat[, j])  # matrix: iterations x chains
    chain_means <- colMeans(chains_j)
    overall_mean <- mean(chain_means)
    
    B <- n * var(chain_means)  # Between-chain variance
    W <- mean(apply(chains_j, 2, var))  # Within-chain variance
    var_hat <- ((n - 1) / n) * W + (1 / n) * B
    rhat_vals[j] <- sqrt(var_hat / W)
  }
  
  # Step 5: Combine all chains into one dataframe for summary stats
  all_samples <- do.call(rbind, lapply(seq_along(chains_post), function(i) {
    df <- as.data.frame(chains_post[[i]])
    df$chain <- i   # Add chain index (1, 2, 3, …)
    df
  }))
  colnames(all_samples)[1:d] <- col_names
  
  # Step 6: Summary statistics
  summary_df <- all_samples %>%
    summarise(across(-chain, list(
      mean = ~mean(.x, na.rm = TRUE),
      lower = ~quantile(.x, 0.025, na.rm = TRUE),
      upper = ~quantile(.x, 0.975, na.rm = TRUE)
    )))
  
  # Step 7: Restructure summary and attach Rhat
  summary_long <- tibble(parameter = col_names) %>%
    mutate(
      mean = unlist(summary_df[1, seq(1, ncol(summary_df), by = 3)]),
      lower = unlist(summary_df[1, seq(2, ncol(summary_df), by = 3)]),
      upper = unlist(summary_df[1, seq(3, ncol(summary_df), by = 3)]),
      Rhat = rhat_vals
    )
  
  return(list(
    posterior_samples = all_samples,
    summary = summary_long
  ))
}

plot_posterior_samples_1 <- function(posterior_samples, bounds){
  
  # Convert to long format
  posterior_long <- posterior_samples %>%
    tidyr::pivot_longer(cols = -chain, names_to = "parameter", values_to = "value")
  
  # Add parameter names to bounds
  bounds <- bounds %>%
    as.data.frame() %>%
    rename(lower = V1, upper = V2) %>%
    mutate(parameter = unique(posterior_long$parameter))
  
  # Merge with posterior samples
  posterior_long <- posterior_long %>%
    left_join(bounds, by = "parameter")
  
  # Plot histograms, one per parameter
  p <- ggplot(posterior_long, aes(x = value, fill = factor(chain))) +
    geom_histogram(bins = 50, color = "white") +
    facet_wrap(~ parameter, scales = "free", ncol = 3) +
    theme_minimal() +
    labs(title = "Posterior Distributions", x = "Value", y = "Count") +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    geom_blank(aes(x = lower)) +
    geom_blank(aes(x = upper)) +
    coord_cartesian()
  p
  
}


plot_ppc_1 <- function(mcmc_outputs, burnin_frac = 0.5, thin= 1, data) {
  # Remove burn-in from each chain, then combine samples from all chains
  all_quantities <- do.call(rbind, lapply(mcmc_outputs, function(x) {
    chain_quantities <- (x$sim_outputs)  # transpose if needed
    burnin <- floor(nrow(chain_quantities) * burnin_frac)
    chain_quantities_post_burnin <- chain_quantities[-seq_len(burnin), , drop = FALSE]
    return(chain_quantities_post_burnin)
  })) %>% as.data.frame() %>% 
    as.data.frame() %>%
    .[seq(1, nrow(.), by = thin), ]  # Keep every 200th row
  
  # Keep only accepted values
  all_quantities <- all_quantities %>%
    filter(accept == 1)
  
  # Quantities present in BOTH model outputs and observed data
  plot_quantities <- intersect(names(data), names(all_quantities))
  
  # Convert to long format with one row per parameter
  summary_long <- all_quantities %>%
    dplyr::select(dplyr::all_of(plot_quantities)) %>%
    tidyr::pivot_longer(cols = everything(),
                        names_to = "quantity",
                        values_to = "value") %>%
    dplyr::group_by(quantity) %>%
    dplyr::summarise(mean  = mean(value, na.rm = TRUE),
                     lower = quantile(value, 0.025, na.rm = TRUE),
                     upper = quantile(value, 0.975, na.rm = TRUE),
                     .groups = "drop")
  
  # Combine with observed data
  summary_data <- data.frame(quantity = names(data),
                             t(sapply(data, function(x) c(lower_obs = min(x), 
                                                          mean_obs = mean(x), 
                                                          upper_obs = max(x))))) 
  
  # Combine observed and modelled quantities
  merged_df <- merge(summary_data, summary_long, 
                     by = "quantity", all = TRUE)
  
  # Normalise by observed mean
  merged_df_norm <- merged_df |>
    mutate(mean_for_norm = mean_obs,
           lower_obs = lower_obs/mean_for_norm, 
           upper_obs = upper_obs/mean_for_norm,
           mean_obs = mean_obs/mean_for_norm,
           lower = lower/mean_for_norm,
           upper = upper/mean_for_norm,
           mean = mean/mean_for_norm)
  
  # Labels for plotting only
  plot_labels <- c(incd_old  = "Overall incidence old",
                   incd_new  = "Overall incidence new",
                   incd_rr   = "RR-TB incidence",
                   tb_mort   = "TB mortality",
                   notifs    = "Notifications",
                   notifs_rr = "RR-TB notifications")
  
  manual_colors <- c("Observed" = "#023e8a", 
                     "Posterior predicted" = "#ef233c")
  
  # Plot observed and predicted values
  p_ppc <- ggplot(merged_df_norm, aes(x = quantity)) +
    # Predicted mean with error bars
    geom_point(aes(y = mean, color = "Posterior predicted"),
               position = position_nudge(x = 0.15)) +
    geom_errorbar(aes(ymin = lower, ymax = upper, color = "Posterior predicted"),
                  width = 0.1, position = position_nudge(x = 0.15)) +
    # Observed values
    geom_point(aes(y = mean_obs, color = "Observed"),
               shape = 17, position = position_nudge(x = -0.15)) +
    geom_errorbar(aes(ymin = lower_obs, ymax = upper_obs, color = "Observed"),
                  width = 0.1, position = position_nudge(x = -0.15)) +
    scale_color_manual(values = manual_colors) +
    scale_x_discrete(labels = plot_labels) +
    guides(color = guide_legend(title = "")) +
    labs(x = "Quantity", y = "Scaled by observed mean") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text.x = element_text(angle = 25, hjust = 0.7))
  
  p_ppc_unnorm <- ggplot(merged_df,
                         aes(x = quantity)) +
    # Predicted mean with error bars
    geom_point(aes(y = mean, color = "Posterior predicted"),
               position = position_nudge(x = 0.15)) +
    geom_errorbar(aes(ymin = lower, ymax = upper, color = "Posterior predicted"),
                  width = 0.1, position = position_nudge(x = 0.15)) +
    # Observed values
    geom_point(aes(y = mean_obs, color = "Observed"),
               shape = 17, position = position_nudge(x = -0.15)) +
    geom_errorbar(aes(ymin = lower_obs, ymax = upper_obs, color = "Observed"),
                  width = 0.1, position = position_nudge(x = -0.15)) +
    scale_color_manual(values = manual_colors) +
    scale_x_discrete(labels = plot_labels) +
    guides(color = guide_legend(title = "")) +
    labs(x = "Quantity", y = "Scaled by observed mean") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text.x = element_text(angle = 25, hjust = 0.7))
  
  return(list(p_ppc = p_ppc, p_ppc_unnorm = p_ppc_unnorm))
}

