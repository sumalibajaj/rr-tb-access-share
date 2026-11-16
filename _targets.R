# Created by use_targets().
library(targets)
library(future)
plan(multisession, workers = 4)

# Set target options:
tar_option_set(
  packages = c("tibble", "dplyr", "ggplot2", "cowplot", "purrr",
               "Matrix", "deSolve", "MASS", "tidyr", "R.utils", "magrittr", 
               "scales", "numDeriv", "stringr"),
  format = "rds")

source("data_exploration.R")
source("get_hiv_data.R")
source("misc_plotting.R")
source("get_addresses.R")
source("get_distribution_fns.R")
source("allocate_parameters.R")
source("goveqs.R")
source("get_objective.R")
source("MCMC_adaptive.R")
source("summarise_outputs.R")
source("get_intervention.R")


list(
  tar_target(country_names, c(  "India",
                                 "China",
                                 "Philippines",
                                 "Indonesia",
                                 "Russian Federation",
                                 "Pakistan",
                                 "South Africa",
                                 "Myanmar",
                                 "Ukraine",
                                 "Viet Nam")),

  # Data exploration
  tar_target(prct_pulm_bac_conf_df, data.frame(country = country_names,
                                               prct_pulm_bac_conf = c(0.63, 0.73,
                                                                      0.46, 0.55,
                                                                      0.57, 0.50,
                                                                      0.62, 0.51,
                                                                      0.75, 0.85))),  
  
  tar_target(country_data_relations,
             extract_country_data_relations(country_names = country_names,
                                            prct_pulm_bac_conf_df = prct_pulm_bac_conf_df)),
  
  tar_target(country_dst_limits, extract_dst_limit(country_names = country_names)),
  
  tar_target(tx_outcomes_all, extract_tx_outcomes(country_names = country_names)),
  
  # Figures - data
  tar_target(combined_burden_int_figure, plot_combined_improvements(tx_outcomes_all = tx_outcomes_all)),
  
  tar_target(save_sl_tx_figure,
             ggsave("outputs/sl_tx_figure.pdf", plot = combined_burden_int_figure$sl_plot,
                    width = 6, height = 8), format = "file"),
  tar_target(save_sl_tx_figure_png,
             ggsave("outputs/sl_txint_figure.png", plot = combined_burden_int_figure$sl_plot,
                    width = 6, height = 8), format = "file"),
  
  tar_target(save_rr_figure,
             ggsave("outputs/rr_figure.pdf", plot = combined_burden_int_figure$burden_plot,
                    width = 8, height = 8), format = "file"),
  tar_target(save_rr_figure_png,
             ggsave("outputs/rr_figure.png", plot = combined_burden_int_figure$burden_plot,
                    width = 8, height = 8), format = "file"),
  
  tar_target(plot_incs_over_time, plot_inc(countries = country_names)$p),
  tar_target(save_plot_incs_over_time,
             ggsave("outputs/inc_over_time.pdf", plot = plot_incs_over_time,
                    width = 15, height = 10), format = "file"),
  
  # tar_target(hiv_sa_1, get_hiv_data("data/Estimates_2025_en.csv", country_name = "South Africa")), # being created later
  tar_target(hiv_ua_1, get_hiv_data_thembisa(country_name = "Ukraine")),
  tar_target(plot_hiv_inc, plot_hiv_incidence(hiv_sa = hiv_sa_1, hiv_ua = hiv_ua_1)),
  tar_target(save_plot_hiv_inc, ggsave("outputs/plot_hiv_inc.pdf", 
                                       plot = plot_hiv_inc,
                                       width = 10, height = 4), format = "file"),
  tar_target(save_plot_hiv_inc_png, ggsave("outputs/plot_hiv_inc.png", 
                                           plot = plot_hiv_inc,
                                           width = 10, height = 4), format = "file"),
  
  # Figure- model (can uncomment and run this once all models are run)
  tar_target(countries_eg, {c("India",
                                         "Philippines",
                                         "South Africa" ,
                                         "China")}),
  tar_target(combined_trace_ll_all, {list(trace_ll_plot_in_1,
                                          trace_ll_plot_ph_1,
                                          trace_ll_plot_sa_1,
                                          trace_ll_plot_ch_1)}),

  tar_target(plot_combined_trace_ll, {
    plots <- Map(function(p, country) {
      p + ggplot2::labs(title = country) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
    },
    combined_trace_ll_all, countries_eg)
    patchwork::wrap_plots(plots, ncol = 2)
  }
  ),
  tar_target(save_plot_combined_trace_ll, ggsave("outputs/plot_combined_trace_ll.pdf",
                                                 plot = plot_combined_trace_ll,
                                                 width = 8, height = 8), format = "file"),

  tar_target(combined_ppc_all, {list(ppc_plot_in_1,
                                     ppc_plot_ph_1,
                                     ppc_plot_sa_1,
                                     ppc_plot_ch_1)}),

  tar_target(plot_combined_ppc, {
    plots <- Map(function(p, country) {
      p + ggplot2::labs(title = country) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
    },
    combined_ppc_all, countries_eg)
    patchwork::wrap_plots(plots, ncol = 1)
  }
  ),
  tar_target(save_plot_combined_ppc, ggsave("outputs/plot_combined_ppc.pdf",
                                            plot = plot_combined_ppc,
                                            width = 10, height = 14), format = "file"),

  # Global aggregation plots and summary table
  tar_target(agg_global, {make_agg_global(draw_list = list(counterfactuals_draws_in_1,
                                                           counterfactuals_draws_ch_1,
                                                           counterfactuals_draws_ph_1,
                                                           counterfactuals_draws_sa_1),
                                          countries = countries_eg)
    }),
  tar_target(save_agg_global_past, ggsave("outputs/agg_global_past.pdf",
                                              plot = agg_global$global_past_plot,
                                              width = 8, height = 4), format = "file"),
  tar_target(save_agg_global_future, ggsave("outputs/agg_global_future.pdf",
                                                plot = agg_global$global_future_plot,
                                                width = 8, height = 4), format = "file"),

  tar_target(summary_agg_global, agg_global$summary_table),

  tar_target(save_summary_agg_global,{
    path <- "outputs/processed data/summary_table_summary_agg_global_unc.csv"
    write.csv(summary_agg_global, path)
    path },
    format = "file"),

  # For making summary table of posterior estimates
  tar_target(summary_list, {list(posterior_summary_in_1,
                                 posterior_summary_ch_1,
                                 posterior_summary_ph_1,
                                 posterior_summary_sa_1)
    }),
  
  # Set iterations, etc
  tar_target(iters_init_screen, 100),
  tar_target(iters_mcmc_stage1, 100),
  tar_target(iters_mcmc_stage2, 100),
  tar_target(intx_sams, 5),
  
  # ----------------------------------------------------------------------------
  # Philippines
  # M2
  # ----------------------------------------------------------------------------
  
  tar_target(states0_m2_1, c("U")),
  tar_target(states1_m2_1, c("Lf", "Ls", "I", "Rlo", "Rhi", "R")),
  tar_target(states2_m2_1, c("Tx", "Tx2")),
  tar_target(gps_m2_1, list(strains = c("ds", "rr"), sectors = c("pu", "pr_nen"))),
  
  # Step-by-step get_addresses() calls
  tar_target(res1_m2_1, get_addresses(list(states0_m2_1))),
  tar_target(res2_m2_1, get_addresses(list(states1_m2_1, gps_m2_1[["strains"]]),
                                       res1_m2_1[["i"]], res1_m2_1[["s"]],
                                       res1_m2_1[["d"]], res1_m2_1[["lim"]])),
  tar_target(res3_m2_1, get_addresses(list(states2_m2_1, gps_m2_1[["strains"]], gps_m2_1[["sectors"]]),
                                       res2_m2_1[["i"]], res2_m2_1[["s"]],
                                       res2_m2_1[["d"]], res2_m2_1[["lim"]])),
  tar_target(ref_m2_1, {
    # i
    i <- res3_m2_1[["i"]]
    lim <- res3_m2_1[["lim"]]
    
    # Add auxiliaries
    names <- c("inc", "noti", "mort_tb", "mort_tb_ds", "mort_tb_rr", "noti_rr")
    lgths <- c(2, 1, 1, 1, 1, 1)  # Number of indices for each auxiliary
    
    aux_list <- list()
    for (j in seq_along(names)) {
      inds <- (lim + 1):(lim + lgths[j])
      aux_list[[names[j]]] <- inds
      lim <- max(inds)
    }
    
    i[["aux"]] <- aux_list
    i[["nx"]] <- lim
    
    # s
    s <- res3_m2_1[["s"]]
    # If you want to create or update infectious/prevalent here:
    s[["infectious"]] <- unique(c(s[["I"]],
                                  intersect(s[["Tx"]], s[["rr"]])))
    s[["prevalent"]] <- unique(c(s[["infectious"]], s[["Tx"]], s[["Tx2"]]))
    
    # d
    d <- unlist(res3_m2_1[["d"]])
    
    list(i = i, s = s, d = d)
    
  }),
  
  tar_target(thetai_m2_1, {
    # Note that elements labelled 'rf' are multipliers on rates (useful when we
    # have a central estimate from literature, but want to incorporate uncertainty)
    list(
      r_beta = c(1, 2),
      r_Tx_init_new = 3,
      rf_Tx_init_old = 4,
      p_pu = 5,
      rf_mort_TB = 6,
      rf_self_cure = 7,
      rf_progression = 8,
      rf_stabilisation = 9,
      rf_reactivation = 10,
      rf_relapse = c(11, 12, 13),
      p_imm = 14,
      p_MDR_rec = 15,
      p_MDR_rec_pr_nen = 16,
      p_SL_trans_pr_nen = 17)
  }),
  
  
  tar_target(prm_ph_1, {
    r <- list(
      # Natural history
      progression0 = 0.0826, # Progression from Lf to I
      stabilisation0 = 0.872, # Stabilisation from Lf to Ls
      reactivation0 = 0.0006, # Reactivation from Ls to I
      self_cure0 = 1/6,
      mort_TB0 = 1/6,
      relapse0 = c(0.032, 0.14, 0.0015),
      mort = 1/66,
      MDR_acqu = 0.01,
      Tx = 2,
      Tx2 = 0.5
    )
    
    p <- list(
      SL_trans = c(0.8) # reporting sectors
    )
    
    # 1 row FL, 2 row SL:  1.Success 2.Failed 3.Died 4.LTFU
    Tx_outc <- tx_outcomes_all[["Philippines"]][c("FL", "SL"), ]
    tmp1 <- Tx_outc / rowSums(Tx_outc) # Row-wise proportions of outcomes
    tmp2 <- cbind(rowSums(tmp1[, 1:2]), tmp1[,3], tmp1[,4]) # Collapse success + failure into "completed"
    p$cure <- tmp1[,1] / tmp2[,1] # Find proportion succeeding, out of those completing treatment
    p$cure[3] <- p$cure[2] # for pr_nen
    
    rTx <- c(r[["Tx"]], r[["Tx2"]]) # Treatment rates
    r$muTx <- (tmp2[, 2] / tmp2[, 1]) * rTx # Mortality rates
    r$ltfu <- (tmp2[, 3] / tmp2[, 1]) * rTx # LTFU rates
    
    r$muTx[3] <- r$muTx[2] # for muTx2 of non-reporting sector (pr_nen)
    
    r$ltfu[3] <- r$ltfu[2] # for ltfu2 of non-reporting sector (pr_nen)
    
    list(p = p, r = r)
  }),
  
  tar_target(prct_bac_conf_ph_1, prct_pulm_bac_conf_df %>%
               dplyr::filter(country == "Philippines") %>%
               dplyr::pull(prct_pulm_bac_conf)),
  
  tar_target(bounds_ph_1, {
    bounds <- matrix(NA, ncol = 2, nrow = length(unlist(thetai_m2_1)))
    
    bounds[thetai_m2_1$r_beta, ] <- matrix(c(1e-6, 1e-6, 50, 50), nrow = 2)
    bounds[thetai_m2_1$r_Tx_init_new, ] <- c(0, 12)
    bounds[thetai_m2_1$rf_Tx_init_old, ] <- c(0, 1)
    bounds[thetai_m2_1$p_pu, ] <- c(0, 1)
    bounds[thetai_m2_1$rf_mort_TB, ] <- c(0, 1.5)
    bounds[thetai_m2_1$rf_self_cure, ] <- c(0.5, 1.5)
    bounds[thetai_m2_1$rf_progression, ] <- c(0.5, 1.5)
    bounds[thetai_m2_1$rf_stabilisation, ] <- c(0.5, 1.5)
    bounds[thetai_m2_1$rf_reactivation, ] <- c(0.5, 1.5)
    bounds[thetai_m2_1$rf_relapse, ] <- matrix(c(0.5, 0.5, 0.5, 1.5, 1.5, 1.5), nrow = 3)
    bounds[thetai_m2_1$p_imm, ] <- c(0.5, 0.9)
    bounds[thetai_m2_1$p_MDR_rec, ] <- c(0, prct_bac_conf_ph_1)
    bounds[thetai_m2_1$p_MDR_rec_pr_nen, ] <- c(0, 0.25) # let it vary between 0 and 0.25
    bounds[thetai_m2_1$p_SL_trans_pr_nen, ] <- c(0, 0.8) # let it vary between 0 and pu level
    
    return(bounds)
  }),
  
  tar_target(extract_data_ph_1, extract_calib_data(country_name = "Philippines", 
                                                   old_inc_year = 2023,
                                                   prct_pulm_bac_conf = prct_bac_conf_ph_1)),
  tar_target(data_ph_1, {
    list(
      incd_new = extract_data_ph_1$incd_new,
      incd_old = extract_data_ph_1$incd_old, 
      incd_rr = extract_data_ph_1$incd_rr,
      tb_mort = extract_data_ph_1$tb_mort,
      notifs = extract_data_ph_1$notifs_pu, # NOT calibrated; range between all notifs and bac confirmed notifs
      notifs_rr = extract_data_ph_1$notifs_rr_pu
    )
  }),
  
  tar_target(data_ph_notif_bounds_1, 
             list(notif_l = extract_data_ph_1$notifs_pu[1],
                  notif_u = extract_data_ph_1$notifs_pu[3])),
  
  tar_target(lhd_incd_old_ph_1, get_distribution_fns(data_ph_1$incd_old, 'lognormal')[[1]]),
  tar_target(lhd_incd_new_ph_1, get_distribution_fns(data_ph_1$incd_new, 'lognormal')[[1]]),
  tar_target(lhd_incd_rr_ph_1, get_distribution_fns(data_ph_1$incd_rr, 'lognormal')[[1]]),
  tar_target(lhd_tb_mort_ph_1, get_distribution_fns(data_ph_1$tb_mort, 'lognormal')[[1]]),
  tar_target(lhd_notifs_ph_1, get_distribution_fns(data_ph_1$notifs, 'lognormal')[[1]]),
  tar_target(lhd_notifs_rr_ph_1, get_distribution_fns(data_ph_1$notifs_rr, 'lognormal')[[1]]),
  
  tar_target(lhd_list_ph_1, list(lhd_incd_old_ph_1,
                                 lhd_incd_new_ph_1 ,
                                 lhd_incd_rr_ph_1, 
                                 lhd_tb_mort_ph_1,
                                 lhd_notifs_ph_1, 
                                 lhd_notifs_rr_ph_1)),
  
  # Not calibrating to all notifications (hard bounded)
  tar_target(lhd_switch_m2m3_1, c(1, 1, 1, 1, 0, 1)),
  
  tar_target(lhd_weight_m2m3_1, c(1, 1, 1, 1, 1, 1)),
  
  tar_target(output_names_m2m3_1, c("incd_old",
                                  "incd_new",
                                  "incd_rr",
                                  "tb_mort",
                                  "tb_mort_ds",
                                  "tb_mort_rr",
                                  "notifs",
                                  "prev",
                                  "notifs_rr")),
  
  tar_target(theta_initial_screening_results_ph_1, {

    print("Finding better space coverage with LHS")
    # 1. Better space coverage with LHS
    library(lhs)
    n_params <- nrow(bounds_ph_1)
    lhs_design <- randomLHS(iters_init_screen, n_params)
    
    theta_matrix <- t(apply(lhs_design, 1, function(x) {
      bounds_ph_1[,1] + x * (bounds_ph_1[,2] - bounds_ph_1[,1])
    }))
    
    # 2. Evaluate all candidates
    print("Evaluating all candidates")
    log_likelihoods <- numeric(nrow(theta_matrix))
    for (i in seq_len(nrow(theta_matrix))) {
      if (i %% 100 == 0 || i == 1) {
        cat("Evaluating theta", i, "of", nrow(theta_matrix), "\n")
      }
      theta <- theta_matrix[i, ]
      log_likelihoods[i] <- get_objective_m2_1(
        theta = theta, 
        allocate_parameters = allocate_parameters_m2_1,
        goveqs = goveqs_m2_1,
        goveqs_scaleup = goveqs_scaleup,
        bounds = bounds_ph_1, 
        prm = prm_ph_1, 
        ref = ref_m2_1, 
        thetai = thetai_m2_1,
        lhd_list = lhd_list_ph_1,
        lhd_switch = lhd_switch_m2m3_1,
        lhd_weight = lhd_weight_m2m3_1,
        data_notif_bounds = data_ph_notif_bounds_1,
        old_inc_year = 2023,
        output_names = output_names_m2m3_1)[["lhd"]]
    }
    
    # 3. Take top 10 (can optimise these but takes time)
    print("Taking top 10")
    results <- data.frame(theta_matrix)
    results$loglikelihood <- log_likelihoods
    sorted_results <- results[order(-results$loglikelihood), ]
    
    # Make sure to keep only mortTB > mort on SL Tx
    sorted_results <- sorted_results %>% filter(.[[thetai_m2_1$rf_mort_TB]] > prm_ph_1$r$muTx[2]/prm_ph_1$r$mort_TB0)
    
    best_theta <- as.numeric(sorted_results[1, 1:n_params])
    
    list(
      sorted_results = sorted_results,
      best_theta = best_theta,
      cov0 = cov(sorted_results[1:10, 1:n_params])
    )
  }),
  
  tar_target(
    top_theta_list_ph_1,
    {
      # Take only the first theta and repeat it 4 times
      first_theta <- theta_initial_screening_results_ph_1$best_theta
      out <- rep(list(first_theta), 4)
      names(out) <- paste0("theta_", seq_along(out))
      out
    }
  ),
  
  tar_target(
    output_mcmc_ph_1,
    {
      # --- Stage 1 ---
      stage1 <- MCMC_adaptive_1(
        obj = get_objective_m2_1,
        allocate_parameters = allocate_parameters_m2_1,
        goveqs = goveqs_m2_1,
        goveqs_scaleup = goveqs_scaleup,
        theta0 = unlist(top_theta_list_ph_1),
        niter = iters_mcmc_stage1,
        mul_sigma = 1.5,
        bounds = bounds_ph_1,
        prm = prm_ph_1,
        ref = ref_m2_1,
        thetai = thetai_m2_1,
        cov0 = theta_initial_screening_results_ph_1$cov0,
        lhd_list = lhd_list_ph_1,
        lhd_switch = lhd_switch_m2m3_1,
        lhd_weight = lhd_weight_m2m3_1,
        data_notif_bounds = data_ph_notif_bounds_1,
        old_inc_year = 2023,
        output_names = output_names_m2m3_1,
        hiv_model = FALSE)
      
      # --- Process Stage 1 Output ---
      best_idx <- which.max(stage1$values)
      theta_best <- stage1$samples[, best_idx]
      pars_idx_for_cov_est <- (round(ncol(stage1$samples)/2):ncol(stage1$samples))
      samples_tail <- stage1$samples[, pars_idx_for_cov_est]
      cov_est <- if (ncol(samples_tail) > 1) cov(t(samples_tail)) else 1e-2 * diag(length(theta_best))
      
      # --- Stage 2 ---
      stage2 <- MCMC_adaptive_1(
        obj = get_objective_m2_1,
        allocate_parameters = allocate_parameters_m2_1,
        goveqs = goveqs_m2_1,
        goveqs_scaleup = goveqs_scaleup,
        theta0 = theta_best,
        niter = iters_mcmc_stage2,
        mul_sigma = 1.5,
        bounds = bounds_ph_1,
        prm = prm_ph_1,
        ref = ref_m2_1,
        thetai = thetai_m2_1,
        cov0 = cov_est,
        lhd_list = lhd_list_ph_1,
        lhd_switch = lhd_switch_m2m3_1,
        lhd_weight = lhd_weight_m2m3_1,
        data_notif_bounds = data_ph_notif_bounds_1,
        old_inc_year = 2023,
        output_names = output_names_m2m3_1,
        hiv_model = FALSE)
      
      return(stage2)
    },
    pattern = map(top_theta_list_ph_1),
    iteration = "list"
  ),
  # Use this to run chains in parallel
  # tar_make_future(output_mcmc_ph_1, workers = 4)
  
  tar_target(trace_plot_ph_1, plot_mcmc_traces_1(output_mcmc_ph_1, thetai_m2_1)),
  tar_target(save_trace_ph_1,
             ggsave("outputs/ph_1_trace.pdf", plot = trace_plot_ph_1,
                    width = 12, height = 8), format = "file"),
  # Posterior densities
  tar_target(trace_ll_plot_ph_1, plot_mcmc_ll_1(output_mcmc_ph_1)),
  tar_target(save_trace_ll_ph_1,
             ggsave("outputs/ph_1_trace_ll.pdf", plot = trace_ll_plot_ph_1,
                    width = 12, height = 8), format = "file"),
  
  # Summary statistics of parameters
  tar_target(posterior_summary_ph_1,
             summarise_posterior_with_rhat_1(output_mcmc_ph_1, thetai_m2_1, burnin_frac = 0.5)$summary),
  # Posterior samples across all chains
  tar_target(posterior_samples_ph_1,
             summarise_posterior_with_rhat_1(output_mcmc_ph_1, thetai_m2_1, burnin_frac = 0.5)$posterior_samples),
  
  # Distribution of posterior samples
  tar_target(posterior_dist_plot_ph_1, plot_posterior_samples_1(posterior_samples_ph_1, bounds = bounds_ph_1)),
  tar_target(save_posterior_dist_plot_ph_1,
             ggsave("outputs/ph_1_posterior_dist.pdf", plot = posterior_dist_plot_ph_1,
                    width = 12, height = 8), format = "file"),
  
  # Posterior predictive checks
  tar_target(ppc_unnorm_plot_ph_1, plot_ppc_1(mcmc_outputs = output_mcmc_ph_1, data = data_ph_1)$p_ppc_unnorm),
  tar_target(save_ppc_unnorm_plot_ph_1,
             ggsave("outputs/ph_1_ppc_unnorm.pdf", plot = ppc_unnorm_plot_ph_1,
                    width = 12, height = 8), format = "file"),
  tar_target(ppc_plot_ph_1, plot_ppc_1(mcmc_outputs = output_mcmc_ph_1, data = data_ph_1)$p_ppc),
  tar_target(save_ppc_plot_ph_1,
             ggsave("outputs/ph_1_ppc.pdf", plot = ppc_plot_ph_1,
                    width = 12, height = 8), format = "file"),
  
  tar_target(counterfactuals_ph_1,
             apply_counterfactuals_m2_1(posterior_samples = posterior_samples_ph_1, 
                                               n_sams = intx_sams,
                                               allocate_parameters = allocate_parameters_m2_1, 
                                               goveqs = goveqs_m2_1, 
                                               goveqs_scaleup = goveqs_scaleup,
                                               bounds = bounds_ph_1,
                                               prm = prm_ph_1, 
                                               ref = ref_m2_1,
                                               thetai = thetai_m2_1,
                                               lhd_list = lhd_list_ph_1,
                                               lhd_switch = lhd_switch_m2_1,
                                               pop_country = pop_2024_ph),
  ),
  
  tar_target(counterfactuals_summary_ph_1, counterfactuals_ph_1$summary_dsol),
  tar_target(counterfactuals_draws_ph_1, counterfactuals_ph_1$draws),
  
  tar_target(past_cf_plot_ph_1, plot_counterfactuals_1(counterfactuals_summary_ph_1)$p_past),
  tar_target(save_past_cf_plot_ph_1,
             ggsave("outputs/ph_1_past_cf.pdf", plot = past_cf_plot_ph_1,
                    width = 8, height = 4), format = "file"),
  tar_target(future_cf_plot_ph_1, plot_counterfactuals_1(counterfactuals_summary_ph_1)$p_future),
  tar_target(save_future_cf_plot_ph_1,
             ggsave("outputs/ph_1_future_cf.pdf", plot = future_cf_plot_ph_1,
                    width = 12, height = 6), format = "file"),
  
  # ----------------------------------------------------------------------------
  # China
  # M3
  # ----------------------------------------------------------------------------
  tar_target(states0_m3_1, c("U")),
  tar_target(states1_m3_1, c("Lf", "Ls", "I", "Rlo", "Rhi", "R")),
  tar_target(states2_m3_1, c("Tx", "Tx2")),
  tar_target(gps_m3_1, list(strains = c("ds", "rr"))),
  
  # Step-by-step get_addresses() calls
  tar_target(res1_m3_1, get_addresses(list(states0_m3_1))),
  tar_target(res2_m3_1, get_addresses(list(states1_m3_1, gps_m3_1[["strains"]]),
                                      res1_m3_1[["i"]], res1_m3_1[["s"]],
                                      res1_m3_1[["d"]], res1_m3_1[["lim"]])),
  tar_target(res3_m3_1, get_addresses(list(states2_m3_1, gps_m3_1[["strains"]]),
                                      res2_m3_1[["i"]], res2_m3_1[["s"]],
                                      res2_m3_1[["d"]], res2_m3_1[["lim"]])),
  tar_target(ref_m3_1, {
    # i
    i <- res3_m3_1[["i"]]
    lim <- res3_m3_1[["lim"]]
    
    # Add auxiliaries
    names <- c("inc", "noti", "mort_tb", "mort_tb_ds", "mort_tb_rr", "noti_rr" )
    lgths <- c(2, 1, 1, 1, 1, 1)  # Number of indices for each auxiliary
    
    aux_list <- list()
    for (j in seq_along(names)) {
      inds <- (lim + 1):(lim + lgths[j])
      aux_list[[names[j]]] <- inds
      lim <- max(inds)
    }
    
    i[["aux"]] <- aux_list
    i[["nx"]] <- lim
    
    # s
    s <- res3_m3_1[["s"]]
    # If you want to create or update infectious/prevalent here:
    s[["infectious"]] <- unique(c(s[["I"]],
                                  intersect(s[["Tx"]], s[["rr"]])))
    s[["prevalent"]] <- unique(c(s[["infectious"]], s[["Tx"]], s[["Tx2"]]))
    
    # d
    d <- unlist(res3_m3_1[["d"]])
    
    list(i = i, s = s, d = d)
    
  }),
  
  tar_target(thetai_m3_1, {
    # Note that elements labelled 'rf' are multipliers on rates (useful when we
    # have a central estimate from literature, but want to incorporate uncertainty)
    list(
      r_beta = c(1, 2),
      r_Tx_init_new = 3,
      rf_Tx_init_old = 4,
      rf_mort_TB = 5,
      rf_self_cure = 6,
      rf_progression = 7,
      rf_stabilisation = 8,
      rf_reactivation = 9,
      rf_relapse = c(10, 11, 12),
      p_imm = 13,
      p_MDR_rec = 14,
      p_mort_Tx_TB = 15)
  }),
  
  tar_target(prm_ch_1, {
    r <- list(
      # Natural history
      progression0 = 0.0826, # Progression from Lf to I # Juan- fast progressors fixed, Menzies
      stabilisation0 = 0.872, # Stabilisation from Lf to Ls # Juan- annual transition prob from fast to slow, Menzies
      reactivation0 = 0.0006, # Reactivation from Ls to I # Juan- slow progressors fixed, Menzies
      self_cure0 = 1/6,
      mort_TB0 = 1/6,
      relapse0 = c(0.032, 0.14, 0.0015), # Juan fixed last two (had only two `Rs`)
      mort = 1/66,
      MDR_acqu = 0.01,
      Tx = 2,
      Tx2 = 0.5
    )
    
    p <- list(
      SL_trans = 0.8
    )
    
    # 1 row FL, 2 row SL- 1.Success 2.Failed 3.Died 4.LTFU
    Tx_outc <- tx_outcomes_all[["China"]][c("FL", "SL"), ]
    tmp1 <- Tx_outc / rowSums(Tx_outc) # Row-wise proportions of outcomes
    tmp2 <- cbind(rowSums(tmp1[, 1:2]), tmp1[,3], tmp1[,4]) # Collapse success + failure into "completed"
    p$cure <- tmp1[,1] / tmp2[,1] # Find proportion succeeding, out of those completing treatment
    
    rTx <- c(r[["Tx"]], r[["Tx2"]]) # Treatment rates
    r$muTx <- (tmp2[, 2] / tmp2[, 1]) * rTx # Mortality rates
    r$ltfu <- (tmp2[, 3] / tmp2[, 1]) * rTx # LTFU rates
    
    list(p = p, r = r)
  }),
  
  tar_target(prct_bac_conf_ch_1, prct_pulm_bac_conf_df %>%
               dplyr::filter(country == "China") %>%
               dplyr::pull(prct_pulm_bac_conf)),
  
  tar_target(bounds_ch_1, {
    bounds <- matrix(NA, ncol = 2, nrow = length(unlist(thetai_m3_1)))
    
    bounds[thetai_m3_1$r_beta, ] <- matrix(c(1e-6, 1e-6, 50, 50), nrow = 2)
    bounds[thetai_m3_1$r_Tx_init_new, ] <- c(1e-6, 12)
    bounds[thetai_m3_1$rf_Tx_init_old, ] <- c(1e-6, 1)
    bounds[thetai_m3_1$rf_mort_TB, ] <- c(1e-6, 1.5)
    bounds[thetai_m3_1$rf_self_cure, ] <- c(0.5, 1.5)
    bounds[thetai_m3_1$rf_progression, ] <- c(0.5, 1.5)
    bounds[thetai_m3_1$rf_stabilisation, ] <- c(0.5, 1.5)
    bounds[thetai_m3_1$rf_reactivation, ] <- c(0.5, 1.5)
    bounds[thetai_m3_1$rf_relapse, ] <- matrix(c(0.5, 0.5, 0.5, 1.5, 1.5, 1.5), nrow = 3)
    bounds[thetai_m3_1$p_imm, ] <- c(0.5, 0.9)
    bounds[thetai_m3_1$p_MDR_rec, ] <- c(1e-6, prct_bac_conf_ch_1)
    bounds[thetai_m3_1$p_mort_Tx_TB, ] <- c(0.1, 0.9)
    
    return(bounds)
  }),
  
  tar_target(extract_data_ch_1, extract_calib_data(country_name = "China",
                                                   old_inc_year = 2015,
                                                   prct_pulm_bac_conf = prct_bac_conf_ch_1)),
  
  tar_target(data_ch_1, {
    list(
      incd_new = extract_data_ch_1$incd_new,
      incd_old = extract_data_ch_1$incd_old,
      incd_rr = extract_data_ch_1$incd_rr,
      tb_mort = extract_data_ch_1$tb_mort,
      notifs = extract_data_ch_1$notifs_pu, # NOT calibrated; range between all notifs and bac confirmed notifs
      notifs_rr = extract_data_ch_1$notifs_rr_pu
    )
  }),
  
  tar_target(data_ch_notif_bounds_1,
             list(notif_l = extract_data_ch_1$notifs_pu[1],
                  notif_u = extract_data_ch_1$notifs_pu[3])),
  
  
  tar_target(lhd_incd_old_ch_1, get_distribution_fns(data_ch_1$incd_old, 'lognormal')[[1]]),
  tar_target(lhd_incd_new_ch_1, get_distribution_fns(data_ch_1$incd_new, 'lognormal')[[1]]),
  tar_target(lhd_incd_rr_ch_1, get_distribution_fns(data_ch_1$incd_rr, 'lognormal')[[1]]),
  tar_target(lhd_tb_mort_ch_1, get_distribution_fns(data_ch_1$tb_mort, 'lognormal')[[1]]),
  tar_target(lhd_notifs_ch_1, get_distribution_fns(data_ch_1$notifs, 'lognormal')[[1]]),
  tar_target(lhd_notifs_rr_ch_1, get_distribution_fns(data_ch_1$notifs_rr, 'lognormal')[[1]]),
  
  tar_target(lhd_list_ch_1, list(lhd_incd_old_ch_1,
                                 lhd_incd_new_ch_1 ,
                                 lhd_incd_rr_ch_1,
                                 lhd_tb_mort_ch_1,
                                 lhd_notifs_ch_1,
                                 lhd_notifs_rr_ch_1)),
  
  tar_target(theta_initial_screening_results_ch_1,
             {
               print("Finding better space coverage with LHS")
               # 1. Better space coverage with LHS
               library(lhs)
               n_params <- nrow(bounds_ch_1)
               lhs_design <- randomLHS(iters_init_screen, n_params)
               
               theta_matrix <- t(apply(lhs_design, 1, function(x) {
                 bounds_ch_1[,1] + x * (bounds_ch_1[,2] - bounds_ch_1[,1])
               }))
               
               # 2. Evaluate all candidates
               print("Evaluating all candidates")
               log_likelihoods <- numeric(nrow(theta_matrix))
               for (i in seq_len(nrow(theta_matrix))) {
                 if (i %% 500 == 0 || i == 1) {
                   cat("Evaluating theta", i, "of", nrow(theta_matrix), "\n")
                 }
                 theta <- theta_matrix[i, ]
                 log_likelihoods[i] <- get_objective_m3_1(
                   theta = theta,
                   allocate_parameters = allocate_parameters_m3_1,
                   goveqs = goveqs_m3_1,
                   goveqs_scaleup = goveqs_scaleup,
                   bounds = bounds_ch_1,
                   prm = prm_ch_1,
                   ref = ref_m3_1,
                   thetai = thetai_m3_1,
                   lhd_list = lhd_list_ch_1,
                   lhd_switch = lhd_switch_m2m3_1,
                   lhd_weight = lhd_weight_m2m3_1,
                   data_notif_bounds = data_ch_notif_bounds_1,
                   old_inc_year = 2015,
                   output_names = output_names_m2m3_1)[["lhd"]]
               }
               
               # 3. Take top 10 (can optimise these but takes time)
               print("Taking top 10")
               results <- data.frame(theta_matrix)
               results$loglikelihood <- log_likelihoods
               sorted_results <- results[order(-results$loglikelihood), ]
               
               # Make sure to keep only mortTB > mort on SL Tx
               sorted_results <- sorted_results %>% filter(.[[thetai_m3_1$rf_mort_TB]] > prm_ch_1$r$muTx[2]/prm_ch_1$r$mort_TB0)
               
               best_theta <- as.numeric(sorted_results[1, 1:n_params])
               
               list(
                 sorted_results = sorted_results,
                 best_theta = best_theta,
                 cov0 = cov(sorted_results[1:10, 1:n_params])
               )
             }
  ),
  
  tar_target(
    top_theta_list_ch_1,
    {
      # Take only the first theta and repeat it 4 times
      first_theta <- theta_initial_screening_results_ch_1$best_theta
      out <- rep(list(first_theta), 4)
      names(out) <- paste0("theta_", seq_along(out))
      out
    }
  ),
  
  tar_target(
    output_mcmc_ch_1,
    {
      # --- Stage 1 ---
      stage1 <- MCMC_adaptive_1(
        obj = get_objective_m3_1,
        allocate_parameters = allocate_parameters_m3_1,
        goveqs = goveqs_m3_1,
        goveqs_scaleup = goveqs_scaleup,
        theta0 = unlist(top_theta_list_ch_1),
        niter = iters_mcmc_stage1,
        mul_sigma = 1.5,
        bounds = bounds_ch_1,
        prm = prm_ch_1,
        ref = ref_m3_1,
        thetai = thetai_m3_1,
        cov0 = theta_initial_screening_results_ch_1$cov0,
        lhd_list = lhd_list_ch_1,
        lhd_switch = lhd_switch_m2m3_1,
        lhd_weight = lhd_weight_m2m3_1,
        data_notif_bounds = data_ch_notif_bounds_1,
        old_inc_year = 2015,
        output_names = output_names_m2m3_1,
        hiv_model = FALSE)
      
      # --- Process Stage 1 Output ---
      best_idx <- which.max(stage1$values)
      theta_best <- stage1$samples[, best_idx]
      pars_idx_for_cov_est <- (round(ncol(stage1$samples)/2):ncol(stage1$samples))
      samples_tail <- stage1$samples[, pars_idx_for_cov_est]
      cov_est <- if (ncol(samples_tail) > 1) cov(t(samples_tail)) else 1e-2 * diag(length(theta_best))
      
      # --- Stage 2 ---
      stage2 <- MCMC_adaptive_1(
        obj = get_objective_m3_1,
        allocate_parameters = allocate_parameters_m3_1,
        goveqs = goveqs_m3_1,
        goveqs_scaleup = goveqs_scaleup,
        theta0 = theta_best,
        niter = iters_mcmc_stage2,
        mul_sigma = 1.5,
        bounds = bounds_ch_1,
        prm = prm_ch_1,
        ref = ref_m3_1,
        thetai = thetai_m3_1,
        cov0 = cov_est,
        lhd_list = lhd_list_ch_1,
        lhd_switch = lhd_switch_m2m3_1,
        lhd_weight = lhd_weight_m2m3_1,
        data_notif_bounds = data_ch_notif_bounds_1,
        old_inc_year = 2015,
        output_names = output_names_m2m3_1,
        hiv_model = FALSE)
      
      return(stage2)
    },
    pattern = map(top_theta_list_ch_1),
    iteration = "list"
  ),
  # Use this to run chains in parallel
  # tar_make_future(output_mcmc_ch_1, workers = 4)
  
  tar_target(trace_plot_ch_1, plot_mcmc_traces_1(output_mcmc_ch_1, thetai_m3_1)),
  tar_target(save_trace_ch_1,
             ggsave("outputs/ch_1_trace.pdf", plot = trace_plot_ch_1,
                    width = 12, height = 8), format = "file"),
  # Posterior densities
  tar_target(trace_ll_plot_ch_1, plot_mcmc_ll_1(output_mcmc_ch_1)),
  tar_target(save_trace_ll_ch_1,
             ggsave("outputs/ch_1_trace_ll.pdf", plot = trace_ll_plot_ch_1,
                    width = 12, height = 8), format = "file"),
  
  # Summary statistics of parameters
  tar_target(posterior_summary_ch_1,
             summarise_posterior_with_rhat_1(output_mcmc_ch_1, thetai_m3_1, burnin_frac = 0.5)$summary),
  # Posterior samples across all chains
  tar_target(posterior_samples_ch_1,
             summarise_posterior_with_rhat_1(output_mcmc_ch_1, thetai_m3_1, burnin_frac = 0.5)$posterior_samples),
  
  # Distribution of posterior samples
  tar_target(posterior_dist_plot_ch_1, plot_posterior_samples_1(posterior_samples_ch_1, bounds = bounds_ch_1)),
  tar_target(save_posterior_dist_plot_ch_1,
             ggsave("outputs/ch_1_posterior_dist.pdf", plot = posterior_dist_plot_ch_1,
                    width = 12, height = 8), format = "file"),
  
  # Posterior predictive checks
  tar_target(ppc_unnorm_plot_ch_1, plot_ppc_1(mcmc_outputs = output_mcmc_ch_1, data = data_ch_1)$p_ppc_unnorm),
  tar_target(save_ppc_unnorm_plot_ch_1,
             ggsave("outputs/ch_1_ppc_unnorm.pdf", plot = ppc_unnorm_plot_ch_1,
                    width = 12, height = 8), format = "file"),
  tar_target(ppc_plot_ch_1, plot_ppc_1(mcmc_outputs = output_mcmc_ch_1, data = data_ch_1)$p_ppc),
  tar_target(save_ppc_plot_ch_1,
             ggsave("outputs/ch_1_ppc.pdf", plot = ppc_plot_ch_1,
                    width = 12, height = 8), format = "file"),
  
  # Intervention projections
  tar_target(pop_2024_ch, {df <- read.csv("data/TB_burden_countries_2025-11-27.csv") %>%
    dplyr::filter(country == "China") %>%
    dplyr::filter(year == 2024) %>%
    dplyr::pull(e_pop_num)
  }),
  
  tar_target(counterfactuals_ch_1,
             apply_counterfactuals_m3_1(posterior_samples = posterior_samples_ch_1,
                                        n_sams = intx_sams,
                                        allocate_parameters = allocate_parameters_m3_1,
                                        goveqs = goveqs_m3_1,
                                        goveqs_scaleup = goveqs_scaleup,
                                        bounds = bounds_ch_1,
                                        prm = prm_ch_1,
                                        ref = ref_m3_1,
                                        thetai = thetai_m3_1,
                                        pop_country = pop_2024_ch),
  ),
  
  tar_target(counterfactuals_summary_ch_1, counterfactuals_ch_1[[2]]),
  tar_target(counterfactuals_draws_ch_1, counterfactuals_ch_1[[1]]),
  
  tar_target(past_cf_plot_ch_1, plot_counterfactuals_1(counterfactuals_summary_ch_1)$p_past),
  tar_target(save_past_cf_plot_ch_1,
             ggsave("outputs/ch_1_past_cf.pdf", plot = past_cf_plot_ch_1,
                    width = 8, height = 4), format = "file"),
  tar_target(future_cf_plot_ch_1, plot_counterfactuals_1(counterfactuals_summary_ch_1)$p_future),
  tar_target(save_future_cf_plot_ch_1,
             ggsave("outputs/ch_1_future_cf.pdf", plot = future_cf_plot_ch_1,
                    width = 12, height = 6), format = "file"),
  
  
  
  # ----------------------------------------------------------------------------
  # South Africa
  # M4
  # ----------------------------------------------------------------------------
  
  tar_target(states0_m4_1, c("U")),
  tar_target(states1_m4_1, c("Lf", "Ls", "I", "Rlo", "Rhi", "R")),
  tar_target(states2_m4_1, c("Tx", "Tx2")),
  tar_target(gps_m4_1, list(strains = c("ds", "rr"), hiv = c("h0", "h1", "hart"))),
  
  # Step-by-step get_addresses() calls
  tar_target(res1_m4_1, get_addresses(list(states0_m4_1, gps_m4_1[["hiv"]]))),
  tar_target(res2_m4_1, get_addresses(list(states1_m4_1, gps_m4_1[["strains"]], gps_m4_1[["hiv"]]),
                                      res1_m4_1[["i"]], res1_m4_1[["s"]],
                                      res1_m4_1[["d"]], res1_m4_1[["lim"]])),
  tar_target(res3_m4_1, get_addresses(list(states2_m4_1, gps_m4_1[["strains"]], gps_m4_1[["hiv"]]),
                                      res2_m4_1[["i"]], res2_m4_1[["s"]],
                                      res2_m4_1[["d"]], res2_m4_1[["lim"]])),
  tar_target(ref_m4_1, {
    # i
    i <- res3_m4_1[["i"]]
    lim <- res3_m4_1[["lim"]]
    
    # Add auxiliaries
    names <- c("inc", "noti", "mort_tb", "mort_tb_ds", "mort_tb_rr", "noti_rr",
               "inc_h1hart", "mort_tb_h0", "mort_tb_h1hart")
    lgths <- c(2, 1, 1, 1, 1, 1,
               1, 1, 1)  # Number of indices for each auxiliary
    
    aux_list <- list()
    for (j in seq_along(names)) {
      inds <- (lim + 1):(lim + lgths[j])
      aux_list[[names[j]]] <- inds
      lim <- max(inds)
    }
    
    i[["aux"]] <- aux_list
    i[["nx"]] <- lim
    
    # s
    s <- res3_m4_1[["s"]]
    # If you want to create or update infectious/prevalent here:
    s[["infectious"]] <- unique(c(s[["I"]],
                                  intersect(s[["Tx"]], s[["rr"]])))
    s[["prevalent"]] <- unique(c(s[["infectious"]], s[["Tx"]], s[["Tx2"]]))
    
    # d
    d <- unlist(res3_m4_1[["d"]])
    
    list(i = i, s = s, d = d)
    
  }),
  
  tar_target(thetai_m4_1, {
    # Note that elements labelled 'rf' are multipliers on rates (useful when we
    # have a central estimate from literature, but want to incorporate uncertainty)
    list(
      r_beta = c(1, 2),
      r_Tx_init_new = 3,
      rf_Tx_init_old = 4,
      rf_mort_TB = c(5, 6),
      rf_self_cure = 7,
      rf_progression = 8,
      rf_stabilisation = 9,
      rf_reactivation = 10,
      rf_relapse = c(11, 12, 13),
      p_imm = 14,
      p_MDR_rec = 15,
      r_ART_init = 16,
      rf_HIV_mort = 17,
      p_HIV_relrate = 18)
  }),
  
  tar_target(hiv_sa_1, get_hiv_data("data/Estimates_2025_en.csv", "South Africa")),
  
  # NOTE - We need treatment outcomes stratified by HIV status also not just FL and SL
  tar_target(prm_sa_1, {
    r <- list(
      # Natural history
      progression0 = 0.0826, # Progression from Lf to I
      stabilisation0 = 0.872, # Stabilisation from Lf to Ls
      reactivation0 = 0.0006, # Reactivation from Ls to I 
      self_cure0 = 1/6,
      mort_TB0 = 1/6,
      relapse0 = c(0.032, 0.14, 0.0015),
      mort = 1/66,
      Dx = 52,
      MDR_acqu = 0.01,
      Tx = 2,
      Tx2 = 0.5
    )
    
    p <- list(
      SL_trans = c(0.8)
    )
    
    # 1 row FL, 2 row SL, 3 row HIV+ FL - 1.Success 2.Failed 3.Died 4.LTFU
    Tx_outc <- tx_outcomes_all[["South Africa"]][c("FL", "SL", "HIV"), ]
    tmp1 <- Tx_outc / rowSums(Tx_outc) # Row-wise proportions of outcomes
    tmp2 <- cbind(rowSums(tmp1[, 1:2]), tmp1[,3], tmp1[,4]) # Collapse success + failure into "completed"
    p$cure <- tmp1[,1] / tmp2[,1] # Find proportion succeeding, out of those completing treatment
    
    rTx <- c(r[["Tx"]], r[["Tx2"]], r[["Tx"]]) # Treatment rates
    r$muTx <- (tmp2[, 2] / tmp2[, 1]) * rTx # Mortality rates
    r$ltfu <- (tmp2[, 3] / tmp2[, 1]) * rTx # LTFU rates
    
    
    # HIV related
    p[["rHIV"]] <- approxfun(
      x = hiv_sa_1$HIV_years,
      y = hiv_sa_1$rHIV,
      method = "linear",
      rule = 2
    ) # HIV incidence per-person function from 1980-2035
    ART_start = hiv_sa_1$ART_start
    
    list(p = p, r = r, ART_start = ART_start)
  }),
  
  
  tar_target(prct_bac_conf_sa_1, prct_pulm_bac_conf_df %>%
               dplyr::filter(country == "South Africa") %>%
               dplyr::pull(prct_pulm_bac_conf)),
  
  tar_target(bounds_sa_1, {
    bounds <- matrix(NA, ncol = 2, nrow = length(unlist(thetai_m4_1)))
    
    bounds[thetai_m4_1$r_beta, ] <- matrix(c(1e-6, 1e-6, 35, 35), nrow = 2)
    bounds[thetai_m4_1$r_Tx_init_new, ] <- c(0, 12)
    bounds[thetai_m4_1$rf_Tx_init_old, ] <- c(0, 1)
    bounds[thetai_m4_1$rf_mort_TB, ] <- matrix(c(0, 0, 1.5, 20), nrow = 2)
    bounds[thetai_m4_1$rf_self_cure, ] <- c(0.5, 1.5)
    bounds[thetai_m4_1$rf_progression, ] <- c(0.5, 1.5)
    bounds[thetai_m4_1$rf_stabilisation, ] <- c(0.5, 1.5)
    bounds[thetai_m4_1$rf_reactivation, ] <- c(0.5, 1.5)
    bounds[thetai_m4_1$rf_relapse, ] <- matrix(c(0.5, 0.5, 0.5, 1.5, 1.5, 1.5), nrow = 3)
    bounds[thetai_m4_1$p_imm, ] <- c(0.5, 0.9)
    bounds[thetai_m4_1$p_MDR_rec, ] <- c(0, prct_bac_conf_sa_1)
    bounds[thetai_m4_1$r_ART_init, ] <- c(0, 10)
    bounds[thetai_m4_1$rf_HIV_mort, ] <- c(1, 10)
    bounds[thetai_m4_1$p_HIV_relrate, ] <- c(1, 100)
    
    return(bounds)
  }),
  
  tar_target(extract_data_sa_1, extract_calib_data(country_name = "South Africa", 
                                                   old_inc_year = 2016,
                                                   prct_pulm_bac_conf = prct_bac_conf_sa_1)),
  
  tar_target(extract_hiv_data_sa_1, extract_hiv_calib_data(country_name = "South Africa")),
  
  tar_target(data_sa_1, {
    list(
      incd_new = extract_data_sa_1$incd_new,
      incd_old = extract_data_sa_1$incd_old,
      incd_rr = extract_data_sa_1$incd_rr,
      tb_mort = extract_data_sa_1$tb_mort,
      notifs = extract_data_sa_1$notifs_pu, # NOT calibrated; range between all notifs and bac confirmed notifs
      notifs_rr = extract_data_sa_1$notifs_rr_pu,
      incd_h1hart = extract_hiv_data_sa_1$incd_h1hart, # total TB incidence in HIV positive people in 2024
      tb_mort_h0 = extract_hiv_data_sa_1$tb_mort_h0, # mortality in 2024 (var = e_mort_exc_tbhiv_100k) in HIV negative people in 2024
      tb_mort_h1hart = extract_hiv_data_sa_1$tb_mort_h1hart, # mortality in 2024 (var = e_mort_tbhiv_100k) in HIV positive people in 2024
      HIV_prev = hiv_sa_1$HIVprev_2024, # hiv prevalence in 2024
      ART_covg = hiv_sa_1$ARTcovg_2024 # art coverage in hiv positive people in 2024
    )
  }),
  
  tar_target(data_sa_notif_bounds_1,
             list(notif_l = extract_data_sa_1$notifs_pu[1],
                  notif_u = extract_data_sa_1$notifs_pu[3])),
  
  tar_target(lhd_incd_old_sa_1, get_distribution_fns(data_sa_1$incd_old, 'lognormal')[[1]]),
  tar_target(lhd_incd_new_sa_1, get_distribution_fns(data_sa_1$incd_new, 'lognormal')[[1]]),
  tar_target(lhd_incd_rr_sa_1, get_distribution_fns(data_sa_1$incd_rr, 'lognormal')[[1]]),
  tar_target(lhd_tb_mort_sa_1, get_distribution_fns(data_sa_1$tb_mort, 'lognormal')[[1]]),
  tar_target(lhd_notifs_sa_1, get_distribution_fns(data_sa_1$notifs, 'lognormal')[[1]]),
  tar_target(lhd_notifs_rr_sa_1, get_distribution_fns(data_sa_1$notifs_rr, 'lognormal')[[1]]),
  tar_target(lhd_incd_h1hart_sa_1, get_distribution_fns(data_sa_1$incd_h1hart, 'lognormal')[[1]]),
  tar_target(lhd_tb_mort_h0_sa_1, get_distribution_fns(data_sa_1$tb_mort_h0, 'lognormal')[[1]]),
  tar_target(lhd_tb_mort_h1hart_sa_1, get_distribution_fns(data_sa_1$tb_mort_h1hart, 'lognormal')[[1]]),
  tar_target(lhd_hiv_prev_sa_1, get_distribution_fns(data_sa_1$HIV_prev, 'lognormal')[[1]]),
  tar_target(lhd_ART_covg_sa_1, get_distribution_fns(data_sa_1$ART_covg, 'beta')[[1]]),
  
  
  tar_target(lhd_list_sa_1, list(lhd_incd_old_sa_1,
                                 lhd_incd_new_sa_1 ,
                                 lhd_incd_rr_sa_1, 
                                 lhd_tb_mort_sa_1,
                                 lhd_notifs_sa_1, 
                                 lhd_notifs_rr_sa_1,
                                 lhd_incd_h1hart_sa_1,
                                 lhd_tb_mort_h0_sa_1,
                                 lhd_tb_mort_h1hart_sa_1,
                                 lhd_hiv_prev_sa_1,
                                 lhd_ART_covg_sa_1)),
  
  tar_target(lhd_switch_m4_1, c(1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1)),
  tar_target(lhd_weight_m4_1, c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
  
  tar_target(output_names_m4_1, c("incd_old",
                                  "incd_new",
                                  "incd_rr",
                                  "tb_mort",
                                  "tb_mort_ds",
                                  "tb_mort_rr",
                                  "notifs",
                                  "prev",
                                  "notifs_rr",
                                  "incd_h1hart",
                                  "tb_mort_h0",
                                  "tb_mort_h1hart",
                                  "hiv_prev",
                                  "ART_covg")),
  
  tar_target(theta_initial_screening_results_sa_1, {

    print("Finding better space coverage with LHS")
    # 1. Better space coverage with LHS
    library(lhs)
    n_params <- nrow(bounds_sa_1)
    lhs_design <- randomLHS(iters_init_screen, n_params)
    
    theta_matrix <- t(apply(lhs_design, 1, function(x) {
      bounds_sa_1[,1] + x * (bounds_sa_1[,2] - bounds_sa_1[,1])
    }))
    
    # 2. Evaluate all candidates
    print("Evaluating all candidates")
    log_likelihoods <- numeric(nrow(theta_matrix))
    for (i in seq_len(nrow(theta_matrix))) {
      if (i %% 100 == 0 || i == 1) {
        cat("Evaluating theta", i, "of", nrow(theta_matrix), "\n")
      }
      theta <- theta_matrix[i, ]
      log_likelihoods[i] <- get_objective_m4_1(
        theta = theta,
        allocate_parameters = allocate_parameters_m4_1,
        goveqs = goveqs_m4_1,
        goveqs_scaleup = goveqs_scaleup_hiv,
        bounds = bounds_sa_1, 
        prm = prm_sa_1, 
        ref = ref_m4_1, 
        thetai = thetai_m4_1,
        lhd_list = lhd_list_sa_1,
        lhd_switch = lhd_switch_m4_1,
        lhd_weight = lhd_weight_m4_1,
        data_notif_bounds = data_sa_notif_bounds_1,
        old_inc_year = 2016,
        output_names = output_names_m4_1)[["lhd"]]
    }
    
    # 3. Take top 10 (can optimise these but takes time)
    print("Taking top 10")
    results <- data.frame(theta_matrix)
    results$loglikelihood <- log_likelihoods
    sorted_results <- results[order(-results$loglikelihood), ]
    
    
    # Make sure to keep only mortTB > mort on SL Tx
    sorted_results <- sorted_results %>% filter(.[[thetai_m4_1$rf_mort_TB[1]]] > prm_sa_1$r$muTx[2]/prm_sa_1$r$mort_TB0)
    
    best_theta <- as.numeric(sorted_results[1, 1:n_params])
    
    list(
      sorted_results = sorted_results,
      best_theta = best_theta,
      cov0 = cov(sorted_results[1:10, 1:n_params])
    )
  }),
  
  tar_target(
    top_theta_list_sa_1,
    {
      # Take only the first theta and repeat it 4 times
      first_theta <- theta_initial_screening_results_sa_1$best_theta
      # first_theta <- theta_initial_screening_results_sa_1$sorted_results[2, 1:nrow(bounds_sa_1)]
      out <- rep(list(first_theta), 4)
      names(out) <- paste0("theta_", seq_along(out))
      out
    }
  ),
  
  tar_target(
    output_mcmc_sa_1,
    {
      # --- Stage 1 ---
      stage1 <- MCMC_adaptive_1(
        obj = get_objective_m4_1,
        allocate_parameters = allocate_parameters_m4_1,
        goveqs = goveqs_m4_1,
        goveqs_scaleup = goveqs_scaleup_hiv,
        theta0 = unlist(top_theta_list_sa_1),
        niter = iters_mcmc_stage1,
        mul_sigma = 1.5,
        bounds = bounds_sa_1,
        prm = prm_sa_1,
        ref = ref_m4_1,
        thetai = thetai_m4_1,
        cov0 = theta_initial_screening_results_sa_1$cov0,
        lhd_list = lhd_list_sa_1,
        lhd_switch = lhd_switch_m4_1,
        lhd_weight = lhd_weight_m4_1,
        data_notif_bounds = data_sa_notif_bounds_1,
        old_inc_year = 2016,
        output_names = output_names_m4_1,
        hiv_model = TRUE)
      
      # --- Process Stage 1 Output ---
      best_idx <- which.max(stage1$values)
      theta_best <- stage1$samples[, best_idx]
      pars_idx_for_cov_est <- (round(ncol(stage1$samples)/2):ncol(stage1$samples))
      samples_tail <- stage1$samples[, pars_idx_for_cov_est]
      cov_est <- if (ncol(samples_tail) > 1) cov(t(samples_tail)) else 1e-2 * diag(length(theta_best))
      
      # --- Stage 2 ---
      stage2 <- MCMC_adaptive_1(
        obj = get_objective_m4_1,
        allocate_parameters = allocate_parameters_m4_1,
        goveqs = goveqs_m4_1,
        goveqs_scaleup = goveqs_scaleup_hiv,
        theta0 = theta_best,
        niter = iters_mcmc_stage2,
        mul_sigma = 1.5,
        bounds = bounds_sa_1,
        prm = prm_sa_1,
        ref = ref_m4_1,
        thetai = thetai_m4_1,
        cov0 = cov_est,
        lhd_list = lhd_list_sa_1,
        lhd_switch = lhd_switch_m4_1,
        lhd_weight = lhd_weight_m4_1,
        data_notif_bounds = data_sa_notif_bounds_1,
        old_inc_year = 2016,
        output_names = output_names_m4_1,
        hiv_model = TRUE)
      
      return(stage2)
    },
    pattern = map(top_theta_list_sa_1),
    iteration = "list"
  ),
  # Use this to run chains in parallel
  # tar_make_future(output_mcmc_sa_1, workers = 4)
  
  tar_target(trace_plot_sa_1, plot_mcmc_traces_1(output_mcmc_sa_1, thetai_m4_1)),
  tar_target(save_trace_sa_1,
             ggsave("outputs/sa_1_trace.pdf", plot = trace_plot_sa_1,
                    width = 12, height = 8), format = "file"),
  # Posterior densities
  tar_target(trace_ll_plot_sa_1, plot_mcmc_ll_1(output_mcmc_sa_1)),
  tar_target(save_trace_ll_sa_1,
             ggsave("outputs/sa_1_trace_ll.pdf", plot = trace_ll_plot_sa_1,
                    width = 12, height = 8), format = "file"),
  
  # Summary statistics of parameters
  tar_target(posterior_summary_sa_1,
             summarise_posterior_with_rhat_1(output_mcmc_sa_1, thetai_m4_1, burnin_frac = 0.5)$summary),
  # Posterior samples across all saains
  tar_target(posterior_samples_sa_1,
             summarise_posterior_with_rhat_1(output_mcmc_sa_1, thetai_m4_1, burnin_frac = 0.5)$posterior_samples),
  
  # Distribution of posterior samples
  tar_target(posterior_dist_plot_sa_1, plot_posterior_samples_1(posterior_samples_sa_1, bounds = bounds_sa_1)),
  tar_target(save_posterior_dist_plot_sa_1,
             ggsave("outputs/sa_1_posterior_dist.pdf", plot = posterior_dist_plot_sa_1,
                    width = 12, height = 8), format = "file"),
  
  # Posterior predictive checks
  tar_target(ppc_unnorm_plot_sa_1, plot_ppc_1(mcmc_outputs = output_mcmc_sa_1, data = data_sa_1)$p_ppc_unnorm),
  tar_target(save_ppc_unnorm_plot_sa_1,
             ggsave("outputs/sa_1_ppc_unnorm.pdf", plot = ppc_unnorm_plot_sa_1,
                    width = 12, height = 8), format = "file"),
  tar_target(ppc_plot_sa_1, plot_ppc_1(mcmc_outputs = output_mcmc_sa_1, data = data_sa_1)$p_ppc),
  tar_target(save_ppc_plot_sa_1,
             ggsave("outputs/sa_1_ppc.pdf", plot = ppc_plot_sa_1,
                    width = 12, height = 8), format = "file"),
  
  # Intervention projections
  tar_target(pop_2024_sa, {df <- read.csv("data/TB_burden_countries_2025-11-27.csv") %>%
    dplyr::filter(country == "South Africa") %>%
    dplyr::filter(year == 2024) %>% 
    dplyr::pull(e_pop_num)
  }),
  
  tar_target(counterfactuals_sa_1,
             apply_counterfactuals_m4_1(posterior_samples = posterior_samples_sa_1, 
                                        n_sams = intx_sams,
                                        allocate_parameters = allocate_parameters_m4_1, 
                                        goveqs = goveqs_m4_1, 
                                        goveqs_scaleup = goveqs_scaleup_hiv,
                                        bounds = bounds_sa_1,
                                        prm = prm_sa_1, 
                                        ref = ref_m4_1,
                                        thetai = thetai_m4_1,
                                        lhd_list = lhd_list_sa_1,
                                        lhd_switch = lhd_switch_m4_1,
                                        pop_country = pop_2024_sa),
  ),
  
  tar_target(counterfactuals_summary_sa_1, counterfactuals_sa_1$summary_dsol),
  tar_target(counterfactuals_draws_sa_1, counterfactuals_sa_1$draws),
  
  tar_target(past_cf_plot_sa_1, plot_counterfactuals_1(counterfactuals_summary_sa_1)$p_past),
  tar_target(save_past_cf_plot_sa_1,
             ggsave("outputs/sa_1_past_cf.pdf", plot = past_cf_plot_sa_1,
                    width = 8, height = 4), format = "file"),
  tar_target(future_cf_plot_sa_1, plot_counterfactuals_1(counterfactuals_summary_sa_1)$p_future),
  tar_target(save_future_cf_plot_sa_1,
             ggsave("outputs/sa_1_future_cf.pdf", plot = future_cf_plot_sa_1,
                    width = 12, height = 6), format = "file"),
  
  # ----------------------------------------------------------------------------
  # India
  # M1
  # ----------------------------------------------------------------------------
  
  tar_target(states0_m1_1, c("U")),
  tar_target(states1_m1_1, c("Lf", "Ls", "I", "Rlo", "Rhi", "R")),
  tar_target(states2_m1_1, c("Tx", "Tx2")),
  tar_target(gps_m1_1, list(strains = c("ds", "rr"), sectors = c("pu", "pr_en", "pr_nen"))),
  
  # Step-by-step get_addresses() calls
  tar_target(res1_m1_1, get_addresses(list(states0_m1_1))),
  tar_target(res2_m1_1, get_addresses(list(states1_m1_1, gps_m1_1[["strains"]]),
                                   res1_m1_1[["i"]], res1_m1_1[["s"]],
                                   res1_m1_1[["d"]], res1_m1_1[["lim"]])),
  tar_target(res3_m1_1, get_addresses(list(states2_m1_1, gps_m1_1[["strains"]], gps_m1_1[["sectors"]]),
                                   res2_m1_1[["i"]], res2_m1_1[["s"]],
                                   res2_m1_1[["d"]], res2_m1_1[["lim"]])),
  tar_target(ref_m1_1, {
    # i
    i <- res3_m1_1[["i"]]
    lim <- res3_m1_1[["lim"]]
    
    # Add auxiliaries
    names <- c("inc", "noti_pu", "mort_tb", "mort_tb_ds", "mort_tb_rr", "noti_pr_en", "noti_rr_pu", "noti_rr_pr_en")
    lgths <- c(2, 1, 1, 1, 1, 1, 1, 1)  # Number of indices for each auxiliary
    
    aux_list <- list()
    for (j in seq_along(names)) {
      inds <- (lim + 1):(lim + lgths[j])
      aux_list[[names[j]]] <- inds
      lim <- max(inds)
    }
    
    i[["aux"]] <- aux_list
    i[["nx"]] <- lim
    
    # s
    s <- res3_m1_1[["s"]]
    # If you want to create or update infectious/prevalent here:
    s[["infectious"]] <- unique(c(s[["I"]],
                                  intersect(s[["Tx"]], s[["rr"]])))
    s[["prevalent"]] <- unique(c(s[["infectious"]], s[["Tx"]], s[["Tx2"]]))
    
    # d
    d <- unlist(res3_m1_1[["d"]])
    
    list(i = i, s = s, d = d)
    
  }),
  
  tar_target(prm_in_1, {
    r <- list(
      # Natural history
      progression0 = 0.0826, # Progression from Lf to I
      stabilisation0 = 0.872, # Stabilisation from Lf to Ls
      reactivation0 = 0.0006, # Reactivation from Ls to I
      self_cure0 = 1/6,
      mort_TB0 = 1/6,
      relapse0 = c(0.032, 0.14, 0.0015),
      mort = 1/66,
      MDR_acqu = 0.01,
      Tx = 2,
      Tx2 = 0.5
    )
    
    p <- list(
      MDR_rec_pr_nen = c(0.1), # only fix for private non-engaged
      SL_trans = c(0.8, 0.5)
    )
    
    # 1 row FL, 2 row SL - 1.Success 2.Failed 3.Died 4.LTFU
    Tx_outc <- tx_outcomes_all[["India"]][c("FL", "SL"), ]
    tmp1 <- Tx_outc / rowSums(Tx_outc) # Row-wise proportions of outcomes
    tmp2 <- cbind(rowSums(tmp1[, 1:2]), tmp1[,3], tmp1[,4]) # Collapse success + failure into "completed"
    p$cure <- tmp1[,1] / tmp2[,1] # Find proportion succeeding, out of those completing treatment
    p$cure[3] <- p$cure[2] # for pr_en
    p$cure[4] <- p$cure[2] # for pr_nen
    
    rTx <- c(r[["Tx"]], r[["Tx2"]]) # Treatment rates
    r$muTx <- (tmp2[, 2] / tmp2[, 1]) * rTx # Mortality rates
    r$ltfu <- (tmp2[, 3] / tmp2[, 1]) * rTx # LTFU rates
    
    r$muTx[3] <- r$muTx[2] # for muTx2 of pr_en
    r$muTx[4] <- r$muTx[2] # for muTx2 of pr_nen
    
    r$ltfu[3] <- r$ltfu[2] # for ltfu2 of pr_en
    r$ltfu[4] <- r$ltfu[2] # for ltfu2 of pr_nen
    
    list(p = p, r = r)
  }),
  
  tar_target(thetai_m1_1, {
    # Note that elements labelled 'rf' are multipliers on rates (useful when we
    # have a central estimate from literature, but want to incorporate uncertainty)
    list(
      r_beta = c(1, 2),
      r_Tx_init_new = 3,
      rf_Tx_init_old = 4,
      p_pu = 5,
      p_pr_en = 6,
      rf_mort_TB = 7,
      rf_self_cure = 8,
      rf_progression = 9,
      rf_stabilisation = 10,
      rf_reactivation = 11,
      rf_relapse = c(12, 13, 14),
      p_imm = 15,
      p_MDR_rec = c(16, 17))
  }),
  
  tar_target(bounds_in_1, {
    bounds <- matrix(NA, ncol = 2, nrow = length(unlist(thetai_m1_1)))
    
    bounds[thetai_m1_1$r_beta, ] <- matrix(c(1e-6, 1e-6, 50, 50), nrow = 2)
    bounds[thetai_m1_1$r_Tx_init_new, ] <- c(0, 12)
    bounds[thetai_m1_1$rf_Tx_init_old, ] <- c(0, 1)
    bounds[thetai_m1_1$p_pu, ] <- c(0, 1)
    bounds[thetai_m1_1$p_pr_en, ] <- c(0, 1)
    bounds[thetai_m1_1$rf_mort_TB, ] <- c(0, 1.5)
    bounds[thetai_m1_1$rf_self_cure, ] <- c(0.5, 1.5)
    bounds[thetai_m1_1$rf_progression, ] <- c(0.5, 1.5)
    bounds[thetai_m1_1$rf_stabilisation, ] <- c(0.5, 1.5)
    bounds[thetai_m1_1$rf_reactivation, ] <- c(0.5, 1.5)
    bounds[thetai_m1_1$rf_relapse, ] <- matrix(c(0.5, 0.5, 0.5, 1.5, 1.5, 1.5), nrow = 3)
    bounds[thetai_m1_1$p_imm, ] <- c(0.5, 0.9)
    bounds[thetai_m1_1$p_MDR_rec, ] <- matrix(c(0, 0, 0.57, 0.28), nrow = 2)
    
    return(bounds)
  }),
  
  tar_target(data_in_1, {
    list(
      incd_new = c(160, 187, 218), # total TB incidence in 2024
      incd_old = c(166, 192, 223), # total TB incidence in 2023
      incd_rr = c(7.3, 8.9, 10.5), # incidence of RR-TB in 2024
      tb_mort = c(16.2, 21.3, 26.5), # mortality in 2024
      notifs_pu = c(74.43, 102.105, 129.78), # (ds+rr) pu notifications in 2023 1853274/1428000000*100000 
      notifs_pr_en = c(12.75, 29.355, 45.96), # (ds+rr) pr_en notification in 2023 656262/1428000000*100000
      notifs_rr_pu = 4.1*c(0.9, 1, 1.1), # 59486/1428000000 * 1e5
      notifs_rr_pr_en = 0.3*c(0.9, 1, 1.1) # 4443/1428000000 * 1e5
    )
  }),
  
  # From India TB report 2024, page 231 for 2023
  tar_target(data_in_notif_bounds_1, 
             list(notif_l_pu    = 74.43, # 1062870/1428000000*100000
                  notif_u_pu    = 129.78, # 1853274/1428000000*100000 
                  notif_l_pr_en = 12.75, # 182063/1428000000*100000
                  notif_u_pr_en = 45.96)), # 656262/1428000000*100000
  
  tar_target(lhd_incd_old, get_distribution_fns(data_in_1$incd_old, 'lognormal')[[1]]),
  tar_target(lhd_incd_new, get_distribution_fns(data_in_1$incd_old, 'lognormal')[[1]]),
  tar_target(lhd_incd_rr, get_distribution_fns(data_in_1$incd_rr, 'lognormal')[[1]]),
  tar_target(lhd_tb_mort, get_distribution_fns(data_in_1$tb_mort, 'lognormal')[[1]]),
  tar_target(lhd_notifs_pu, get_distribution_fns(data_in_1$notifs_pu, 'lognormal')[[1]]),
  tar_target(lhd_notifs_pr_en, get_distribution_fns(data_in_1$notifs_pr_en, 'lognormal')[[1]]),
  tar_target(lhd_notifs_rr_pu, get_distribution_fns(data_in_1$notifs_rr_pu, 'lognormal')[[1]]),
  tar_target(lhd_notifs_rr_pr_en, get_distribution_fns(data_in_1$notifs_rr_pr_en, 'lognormal')[[1]]),
  
  tar_target(lhd_list_in_1, list(lhd_incd_old,
                                 lhd_incd_new,
                                 lhd_incd_rr, 
                                 lhd_tb_mort,
                                 lhd_notifs_pu, 
                                 lhd_notifs_pr_en,
                                 lhd_notifs_rr_pu, 
                                 lhd_notifs_rr_pr_en)),
  
  tar_target(lhd_switch_m1_1, c(1, 1, 1, 1, 0, 0, 1, 1)),
  tar_target(lhd_weight_m1_1, c(1, 1, 1, 1, 1, 1, 1, 1)),
  
  tar_target(output_names_m1_1, c("incd_old",
                                  "incd_new",
                                  "incd_rr",
                                  "tb_mort",
                                  "tb_mort_ds",
                                  "tb_mort_rr",
                                  "notifs_pu",
                                  "prev",
                                  "notifs_pr_en",
                                  "notifs_rr_pu",
                                  "notifs_rr_pr_en")),
  
  tar_target(theta_initial_screening_results_in_1, {

    print("Finding better space coverage with LHS")
    # 1. Better space coverage with LHS
    library(lhs)
    n_params <- nrow(bounds_in_1)
    lhs_design <- randomLHS(iters_init_screen, n_params)
    
    theta_matrix <- t(apply(lhs_design, 1, function(x) {
      bounds_in_1[,1] + x * (bounds_in_1[,2] - bounds_in_1[,1])
    }))
    
    # 2. Evaluate all candidates
    print("Evaluating all candidates")
    log_likelihoods <- numeric(nrow(theta_matrix))
    for (i in seq_len(nrow(theta_matrix))) {
      if (i %% 100 == 0 || i == 1) {
        cat("Evaluating theta", i, "of", nrow(theta_matrix), "\n")
      }
      theta <- theta_matrix[i, ]
      log_likelihoods[i] <- get_objective_m1_1(
        theta = theta, 
        allocate_parameters = allocate_parameters_m1_1,
        goveqs = goveqs_m1_1,
        goveqs_scaleup = goveqs_scaleup,
        bounds = bounds_in_1, 
        prm = prm_in_1, 
        ref = ref_m1_1, 
        thetai = thetai_m1_1,
        lhd_list = lhd_list_in_1,
        lhd_switch = lhd_switch_m1_1,
        lhd_weight = lhd_weight_m1_1,
        data_notif_bounds = data_in_notif_bounds_1,
        old_inc_year = 2023,
        output_names = output_names_m1_1)[["lhd"]]
    }
    
    # 3. Take top 10 (can optimise these but takes time)
    print("Taking top 10")
    results <- data.frame(theta_matrix)
    results$loglikelihood <- log_likelihoods
    sorted_results <- results[order(-results$loglikelihood), ]
    
    # Make sure to keep only mortTB > mort on SL Tx
    sorted_results <- sorted_results %>% filter(.[[thetai_m1_1$rf_mort_TB]] > prm_in_1$r$muTx[2]/prm_in_1$r$mort_TB0)
    
    best_theta <- as.numeric(sorted_results[1, 1:n_params])
    
    list(
      sorted_results = sorted_results,
      best_theta = best_theta,
      cov0 = cov(sorted_results[1:10, 1:n_params])
    )
  }),
  
  tar_target(
    top_theta_list_in_1,
    {
      # Take only the first theta and repeat it 4 times
      first_theta <- theta_initial_screening_results_in_1$best_theta
      out <- rep(list(first_theta), 4)
      names(out) <- paste0("theta_", seq_along(out))
      out
    }
  ),
  
  tar_target(
    output_mcmc_in_1,
    {
      # --- Stage 1 ---
      stage1 <- MCMC_adaptive_1(
        obj = get_objective_m1_1,
        allocate_parameters = allocate_parameters_m1_1,
        goveqs = goveqs_m1_1,
        goveqs_scaleup = goveqs_scaleup,
        theta0 = unlist(top_theta_list_in_1),
        niter = iters_mcmc_stage1,
        mul_sigma = 1.5,
        bounds = bounds_in_1,
        prm = prm_in_1,
        ref = ref_m1_1,
        thetai = thetai_m1_1,
        cov0 = theta_initial_screening_results_in_1$cov0,
        lhd_list = lhd_list_in_1,
        lhd_switch = lhd_switch_m1_1,
        lhd_weight = lhd_weight_m1_1,
        data_notif_bounds = data_in_notif_bounds_1,
        old_inc_year = 2023,
        output_names = output_names_m1_1,
        hiv_model = FALSE)
      
      # --- Process Stage 1 Output ---
      best_idx <- which.max(stage1$values)
      theta_best <- stage1$samples[, best_idx]
      pars_idx_for_cov_est <- (round(ncol(stage1$samples)/2):ncol(stage1$samples))
      samples_tail <- stage1$samples[, pars_idx_for_cov_est]
      cov_est <- if (ncol(samples_tail) > 1) cov(t(samples_tail)) else 1e-2 * diag(length(theta_best))
      
      # --- Stage 2 ---
      stage2 <- MCMC_adaptive_1(
        obj = get_objective_m1_1,
        allocate_parameters = allocate_parameters_m1_1,
        goveqs = goveqs_m1_1,
        goveqs_scaleup = goveqs_scaleup,
        theta0 = theta_best,
        niter = iters_mcmc_stage2,
        mul_sigma = 1.5,
        bounds = bounds_in_1,
        prm = prm_in_1,
        ref = ref_m1_1,
        thetai = thetai_m1_1,
        cov0 = cov_est,
        lhd_list = lhd_list_in_1,
        lhd_switch = lhd_switch_m1_1,
        lhd_weight = lhd_weight_m1_1,
        data_notif_bounds = data_in_notif_bounds_1,
        old_inc_year = 2023,
        output_names = output_names_m1_1,
        hiv_model = FALSE)
      
      return(stage2)
    },
    pattern = map(top_theta_list_in_1),
    iteration = "list"
  ),
  # Use this to run chains in parallel
  # tar_make_future(output_mcmc_in_1, workers = 4)
  
  
  tar_target(trace_plot_in_1, plot_mcmc_traces_1(output_mcmc_in_1, thetai_m1_1)),
  tar_target(save_trace_in_1,
             ggsave("outputs/in_1_trace.pdf", plot = trace_plot_in_1,
                    width = 12, height = 8), format = "file"),
  # Posterior densities
  tar_target(trace_ll_plot_in_1, plot_mcmc_ll_1(output_mcmc_in_1)),
  tar_target(save_trace_ll_in_1,
             ggsave("outputs/in_1_trace_ll.pdf", plot = trace_ll_plot_in_1,
                    width = 12, height = 8), format = "file"),
  
  # Summary statistics of parameters
  tar_target(posterior_summary_in_1,
             summarise_posterior_with_rhat_1(output_mcmc_in_1, thetai_m1_1, burnin_frac = 0.5)$summary),
  # Posterior samples across all saains
  tar_target(posterior_samples_in_1,
             summarise_posterior_with_rhat_1(output_mcmc_in_1, thetai_m1_1, burnin_frac = 0.5)$posterior_samples),
  
  # Distribution of posterior samples
  tar_target(posterior_dist_plot_in_1, plot_posterior_samples_1(posterior_samples_in_1, bounds = bounds_in_1)),
  tar_target(save_posterior_dist_plot_in_1,
             ggsave("outputs/in_1_posterior_dist.pdf", plot = posterior_dist_plot_in_1,
                    width = 12, height = 8), format = "file"),
  
  # Posterior predictive checks
  tar_target(ppc_unnorm_plot_in_1, plot_ppc_1(mcmc_outputs = output_mcmc_in_1, data = data_in_1)$p_ppc_unnorm),
  tar_target(save_ppc_unnorm_plot_in_1,
             ggsave("outputs/in_1_ppc_unnorm.pdf", plot = ppc_unnorm_plot_in_1,
                    width = 12, height = 8), format = "file"),
  tar_target(ppc_plot_in_1, plot_ppc_1(mcmc_outputs = output_mcmc_in_1, data = data_in_1)$p_ppc),
  tar_target(save_ppc_plot_in_1,
             ggsave("outputs/in_1_ppc.pdf", plot = ppc_plot_in_1,
                    width = 12, height = 8), format = "file"),
  
  tar_target(pop_2024_in, {df <- read.csv("data/TB_burden_countries_2025-11-27.csv") %>%
    dplyr::filter(country == "India") %>%
    dplyr::filter(year == 2024) %>%
    dplyr::pull(e_pop_num)
  }),
  
  tar_target(counterfactuals_in_1,
             apply_counterfactuals_m1_1(posterior_samples = posterior_samples_in_1, 
                                            n_sams = intx_sams,
                                            allocate_parameters = allocate_parameters_m1_1,
                                            goveqs = goveqs_m1_1,
                                            goveqs_scaleup = goveqs_scaleup,
                                            bounds = bounds_in_1,
                                            prm = prm_in_1, 
                                            ref = ref_m1_1,
                                            thetai = thetai_m1_1,
                                            lhd_list = lhd_list_in_1,
                                            lhd_switch = lhd_switch_m1_1,
                                            pop_country = pop_2024_in)
  ),
  
  tar_target(counterfactuals_summary_in_1, counterfactuals_in_1$summary_dsol),
  tar_target(counterfactuals_draws_in_1, counterfactuals_in_1$draws),
  
  tar_target(past_cf_plot_in_1, plot_counterfactuals_1(counterfactuals_summary_in_1)$p_past),
  tar_target(save_past_cf_plot_in_1,
             ggsave("outputs/in_1_past_cf.pdf", plot = past_cf_plot_in_1,
                    width = 8, height = 4), format = "file"),
  
  tar_target(future_cf_plot_in_1, plot_counterfactuals_1(counterfactuals_summary_in_1)$p_future),
  tar_target(save_future_cf_plot_in_1,
             ggsave("outputs/in_1_future_cf.pdf", plot = future_cf_plot_in_1,
                    width = 12, height = 6), format = "file"),
  
  # India sensitivity of not increasing engagement
  
  tar_target(counterfactuals_in_1_sens,
             apply_counterfactuals_m1_1_sens(posterior_samples = posterior_samples_in_1, 
                                        n_sams = intx_sams,
                                        allocate_parameters = allocate_parameters_m1_1,
                                        goveqs = goveqs_m1_1,
                                        goveqs_scaleup = goveqs_scaleup,
                                        bounds = bounds_in_1,
                                        prm = prm_in_1, 
                                        ref = ref_m1_1,
                                        thetai = thetai_m1_1,
                                        lhd_list = lhd_list_in_1,
                                        lhd_switch = lhd_switch_m1_1,
                                        pop_country = pop_2024_in)
  ),
  
  tar_target(counterfactuals_summary_in_1_sens, counterfactuals_in_1_sens$summary_dsol),
  tar_target(counterfactuals_draws_in_1_sens, counterfactuals_in_1_sens$draws),
  
  tar_target(future_cf_plot_in_1_sens, plot_counterfactuals_1(counterfactuals_summary_in_1_sens)$p_future),
  tar_target(save_future_cf_plot_in_1_sens,
             ggsave("outputs/in_1_sens_future_cf.pdf", plot = future_cf_plot_in_1_sens,
                    width = 12, height = 6), format = "file"),
  tar_target(out_india_sensitivity, 
             calculate_india_sensitivity(primary_draws = counterfactuals_in_1$draws,
                                         sensitivity_draws = counterfactuals_in_1_sens$draws)),
  
  tar_target(future_cf_plot_in_sensitivity_combined,
             patchwork::wrap_plots(future_cf_plot_in_1,
                                   future_cf_plot_in_1_sens,
                                   ncol = 1) +
               patchwork::plot_annotation(tag_levels = "A")),
  
  tar_target(save_future_cf_plot_in_sensitivity_combined,
             ggsave("outputs/future_cf_plot_in_sensitivity_combined.pdf", plot = future_cf_plot_in_sensitivity_combined,
                    width = 8, height = 8), format = "file")
  
  
)