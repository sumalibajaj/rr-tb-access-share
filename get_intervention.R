apply_counterfactuals_m1_1 <- function(posterior_samples, 
                                           n_sams,
                                           allocate_parameters,
                                           goveqs,
                                           goveqs_scaleup,
                                           bounds,
                                           prm, 
                                           ref, 
                                           thetai, 
                                           lhd_list,
                                           lhd_switch,
                                           pop_country){
  # Sample evenly spaced 250 parameters vectors from the posteriors after burnin and thinning
  total_samples <- nrow(posterior_samples)
  dx <- round(total_samples/n_sams)
  indices <- seq(1, total_samples, by = dx)
  indices <- indices[seq_len(min(n_sams, length(indices)))] # Ensure exactly 250 samples
  
  xs <- posterior_samples[indices, 1:nrow(bounds)]
  n_samples <- nrow(xs)
  dsol_list <- vector("list", length = n_samples)
  
  # Main loop over all selected posterior samples
  for (ii in seq_len(n_samples)) {
    # Print to see progress of intervention iterations
    if(ii %% 50 == 0) 
      print(paste("Intervention iteration: ", ii)) 
    
    theta <- as.numeric(xs[ii, ])
    
    # Allocate parameters based on input theta (parameter vector)
    alloc <- allocate_parameters(theta, thetai, prm$p, prm$r)
    r <- alloc$r
    p <- alloc$p
    i <- ref$i
    d <- ref$d
    
    # 1) Actual past trajectory (with RNTCP)
    init <- rep(0, i[["nx"]])
    names(init) <- c(
      "U", "Lf_ds", "Lf_rr", "Ls_ds", "Ls_rr", "I_ds", "I_rr",
      "Rlo_ds", "Rlo_rr", "Rhi_ds", "Rhi_rr", "R_ds", "R_rr",
      "Tx_ds_pu", "Tx_ds_pr_en", "Tx_ds_pr_nen", "Tx_rr_pu", "Tx_rr_pr_en", "Tx_rr_pr_nen",
      "Tx2_ds_pu", "Tx2_ds_pr_en", "Tx2_ds_pr_nen", "Tx2_rr_pu", "Tx2_rr_pr_en", "Tx2_rr_pr_nen",
      "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort", "cu_mort_ds", "cu_mort_rr")
    seed <- 1e-6
    init[i[["U"]]] <- 1 - seed; init[i[["I"]][["ds"]]] <- seed
    
    # equilibrium (no MDR, no RNTCP)
    # In absence of MDR and RNTCP
    p0 <- p; r0 <- r 
    p0[["pu"]] <- 0
    r0[["beta"]][2] <- 0
    r0[["MDR_acqu"]] <- 0
    r0[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters0 <- list(r = r0, p = p0)
    soln0 <- deSolve::ode(y = init, times = seq(0, 2000, by = 1000), func = goveqs, 
                          parms = parameters0,
                          method = "bdf")
    init1 <- tail(soln0[,-1],1)
    
    # MDR intro (no RNTCP)
    p1 <- p; r1 <- r 
    p1[["pu"]] <- 0
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters1 <- list(r = r1, p = p1)
    soln1 <- deSolve::ode(y = init1, times = seq(1970, 1997, by = 1), func = goveqs, 
                          parms = parameters1,
                          method = "bdf")
    init2 <- tail(soln1[,-1],1)
    
    # RNTCP scale-up
    p2 <- p; r2 <- r
    r2[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters2 <- list(r_start = r2, r_end = r2, p_start = p1, p_end = p2, times = c(1997, 2007))
    parameters2$goveqs <- goveqs
    soln2 <- deSolve::ode(y = init2, times = seq(1997, 2015, by = 1), func = goveqs_scaleup, 
                          parms = parameters2,
                          method = "bdf")
    init3 <- tail(soln2[,-1],1)
    
    # Treatment initiation scaleup
    r3_1 <- r
    r3_1[["Tx_init"]] <- r[["Tx_init_old"]]
    r3_2 <- r
    r3_2[["Tx_init"]] <- r[["Tx_init_new"]]
    p3 <- p
    parameters3 <- list(r_start = r3_1, r_end = r3_2, p_start = p3, p_end = p3, times = c(2015, 2025))
    parameters3$goveqs <- goveqs
    soln3 <- deSolve::ode(y = init3, times = seq(2015, 2025, by = 1), func = goveqs_scaleup, 
                          parms = parameters3,
                          method = "bdf")
    
    solution_df0 <- as.data.frame(soln0)
    colnames(solution_df0) <- c("time", d[1:i$nstates], 
                                "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df1 <- as.data.frame(soln1)
    colnames(solution_df1) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df1 <- solution_df1 |> dplyr::filter(time > 1970)
    solution_df2 <- as.data.frame(soln2)
    colnames(solution_df2) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df2 <- solution_df2 |> dplyr::filter(time > 1997)
    solution_df3 <- as.data.frame(soln3)
    colnames(solution_df3) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df3 <- solution_df3 |> dplyr::filter(time > 2015)
    
    sol_actual <- rbind(solution_df1, solution_df2, solution_df3)
    
    dsol_actual <- sol_actual %>%
      dplyr::select(time, cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr) %>%
      mutate(across(-time, ~ c(.x[1], diff(.x))),
             across(-time, ~ .x*1e5),
             scenario="Past actual") %>%
      rename(tb_mort = cu_mort_tb,
             inc_rr = cu_inc_rr,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    init_forward <- tail(soln3[,-1],1)  # forward intervention scenarios use this
    
    # 2) Counterfactual past (no p_MDR_rec, everything else same)
    seed <- 1e-6
    init[i[["U"]]] <- 1 - seed; init[i[["I"]][["ds"]]] <- seed
    
    # equilibrium (no MDR, no RNTCP)
    # In absence of MDR and RNTCP
    p0 <- p
    r0 <- r 
    p0[["pu"]] <- 0
    r0[["beta"]][2] <- 0
    r0[["MDR_acqu"]] <- 0
    r0[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters0 <- list(r = r0, p = p0)
    soln0 <- deSolve::ode(y = init, times = seq(0, 2000, by = 1000), func = goveqs, 
                          parms = parameters0,
                          method = "bdf")
    init1 <- tail(soln0[,-1],1)
    
    # MDR intro (no RNTCP)
    p1 <- p; r1 <- r 
    p1[["pu"]] <- 0
    p1[["MDR_rec"]][1] <- 0
    p1[["MDR_rec"]][2] <- 0
    p1[["MDR_rec_pr_nen"]] <- 0
    r1[["muTx"]][2] <- r[["mort_TB"]]
    r1[["muTx"]][3] <- r[["mort_TB"]]
    r1[["muTx"]][4] <- r[["mort_TB"]]
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters1 <- list(r = r1, p = p1)
    soln1 <- deSolve::ode(y = init1, times = seq(1970, 1997, by = 1), func = goveqs, 
                          parms = parameters1,
                          method = "bdf")
    init2 <- tail(soln1[,-1],1)
    
    # RNTCP scale-up
    p2 <- p; r2 <- r
    p2[["MDR_rec"]][1] <- 0
    p2[["MDR_rec"]][2] <- 0
    p2[["MDR_rec_pr_nen"]] <- 0
    r2[["muTx"]][2] <- r[["mort_TB"]]
    r2[["muTx"]][3] <- r[["mort_TB"]]
    r2[["muTx"]][4] <- r[["mort_TB"]]
    r2[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters2 <- list(r_start = r2, r_end = r2, p_start = p1, p_end = p2, times = c(1997, 2007))
    parameters2$goveqs <- goveqs
    soln2 <- deSolve::ode(y = init2, times = seq(1997, 2015, by = 1), func = goveqs_scaleup, 
                          parms = parameters2,
                          method = "bdf")
    init3 <- tail(soln2[,-1],1)
    
    # Treatment initiation scaleup
    r3_1 <- r
    r3_1[["Tx_init"]] <- r[["Tx_init_old"]]
    r3_2 <- r
    r3_2[["Tx_init"]] <- r[["Tx_init_new"]]
    p3 <- p
    p3[["MDR_rec"]][1] <- 0
    p3[["MDR_rec"]][2] <- 0
    p3[["MDR_rec_pr_nen"]] <- 0
    r3_1[["muTx"]][2] <- r[["mort_TB"]]
    r3_2[["muTx"]][2] <- r[["mort_TB"]]
    r3_1[["muTx"]][3] <- r[["mort_TB"]]
    r3_2[["muTx"]][3] <- r[["mort_TB"]]
    r3_1[["muTx"]][4] <- r[["mort_TB"]]
    r3_2[["muTx"]][4] <- r[["mort_TB"]]
    parameters3 <- list(r_start = r3_1, r_end = r3_2, p_start = p3, p_end = p3, times = c(2015, 2025))
    parameters3$goveqs <- goveqs
    soln3 <- deSolve::ode(y = init3, times = seq(2015, 2025, by = 1), func = goveqs_scaleup, 
                          parms = parameters3,
                          method = "bdf")
    
    solution_df0 <- as.data.frame(soln0)
    colnames(solution_df0) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df1 <- as.data.frame(soln1)
    colnames(solution_df1) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df1 <- solution_df1 |> dplyr::filter(time > 1970)
    solution_df2 <- as.data.frame(soln2)
    colnames(solution_df2) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df2 <- solution_df2 |> dplyr::filter(time > 1997)
    solution_df3 <- as.data.frame(soln3)
    colnames(solution_df3) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df3 <- solution_df3 |> dplyr::filter(time > 2015)
    
    sol_cf <- rbind(solution_df1, solution_df2, solution_df3)
    
    dsol_cf <- sol_cf %>%
      as.data.frame() %>%
      dplyr::select(time, cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr) %>%
      mutate(across(-time, ~ c(.x[1], diff(.x))),
             across(-time, ~ .x*1e5),
             scenario="Past no p_MDR_rec")%>%
      rename(tb_mort = cu_mort_tb,
             inc_rr = cu_inc_rr,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    # 3) Forward projections using actual trajectory till 2025
    # add 0s for capturing treatment starters
    init <- c(init_forward)
    
    # scenario 0 - base model
    r0 <- r
    p0 <- p
    # p0[["MDR_rec_pr_nen"]] <- p0[["MDR_rec"]][2] # I think this was done to maintain that nen is not better than en? very similar
    r0[["Tx_init"]] <- r[["Tx_init_new"]]
    parameters0 <- list(r = r0, p = p0)
    soln0_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs, 
                            parms = parameters0,
                            method = "bdf")
    
    # scenario 1 - improve DST as final condition in scaleup (and private sector engagement)
    p1 <- p0
    r1 <- r0
    r1[["Tx_init"]] <- r[["Tx_init_new"]]
    # Improve public-sector DST recognition
    p1$MDR_rec[1] <- max(p1$MDR_rec[1], 0.75)
    
    # Expand combined public + engaged-private coverage to at least 90%
    target_coverage <- 0.90
    
    p1$pr_en <- if (p1$pu >= target_coverage) {
      p1$pr_en
    } else {
      max(p1$pr_en, (target_coverage - p1$pu) / (1 - p1$pu))
    }
    
    # Improve engaged-private DST recognition to at least public-sector level
    p1$MDR_rec[2] <- max(p1$MDR_rec[2], p1$MDR_rec[1])
    
    parameters1 <- list(r_start = r0, r_end = r1, p_start = p0, p_end = p1, times = c(2026, 2028))
    parameters1$goveqs <- goveqs
    soln1_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs_scaleup, 
                            parms = parameters1,
                            method = "bdf")
    
    # scenario 2 - improve SL treatment as final condition in public sector only (and private sector engagement)
    # a) shift 50% of bad outcomes (LTFU and mortality) into success
    # b) halve LTFU when on Tx2
    # c) halve mortality when on Tx2
    p2 <- p1; r2 <- r1
    vec <- c(r2$Tx2, r2$ltfu[2], r2$muTx[2])
    propvec <- vec/sum(vec) # current proportions of Tx2 outcomes
    newvec <- c(propvec[1] + 0.5 * (1 - propvec[1]), 0.5 * propvec[2], 0.5 * propvec[3]) # new props of Tx2 outcomes
    r2$ltfu[2] <- r2$Tx2 * newvec[2] / newvec[1] # scale lower LTFU rate with new props
    r2$muTx[2] <- r2$Tx2 * newvec[3] / newvec[1] # scale lower muTx2 rate with new props
    p2$cure[2] <- (p2$cure[2] + 1)*0.5 # reducing failure by 50% i.e increasing success. r_Tx2 remains the same
    
    # improving SL outcomes in pr-en as done in public above
    r2$ltfu[3] <- r2$ltfu[2]
    r2$muTx[3] <- r2$muTx[2]
    p2$cure[3] <- p2$cure[2]
    
    parameters2 <- list(r_start = r0, r_end = r2, p_start = p0, p_end = p2, times = c(2026, 2028))
    parameters2$goveqs <- goveqs
    soln2_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs_scaleup, 
                            parms = parameters2,
                            method = "bdf")
    
    # 4) Combine outputs
    soln0_f <- as.data.frame(soln0_f) %>% mutate(scenario = "Baseline")
    soln1_f <- as.data.frame(soln1_f) %>% mutate(scenario = "Improved DST")
    soln2_f <- as.data.frame(soln2_f) %>% mutate(scenario = "Improved DST + SL outcomes")
    
    soln_f <- bind_rows(soln0_f, soln1_f, soln2_f)
    colnames(soln_f) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", 
                          "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr",
                          "cu_noti_pr_en",
                          "cu_noti_rr_pu",
                          "cu_noti_rr_pr_en",
                          "scenario")
    
    dsol_f <- soln_f %>%
      dplyr::select("time", "scenario", "cu_inc_rr", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr") %>%
      group_by(scenario) %>%
      mutate(across(c(cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr),
                    ~ c(.x[1], diff(.x))),
             across(-time, ~ .x * 1e5)) %>% # to get per 100k per year
      rename(inc_rr = cu_inc_rr,
             tb_mort = cu_mort_tb,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    dsol <- bind_rows(dsol_actual, dsol_cf, dsol_f)
    dsol$iter <- ii
    dsol_list[[ii]] <- dsol
  }
  
  combined_dsol_all <- bind_rows(dsol_list)
  
  # Summarise inc_rr, tb_mort, tb_mort_ds, tb_mort_rr in the past and future scenarios
  # pivot long
  long_dsol <- combined_dsol_all %>%
    pivot_longer(cols=c(inc_rr, tb_mort, tb_mort_ds, tb_mort_rr),
                 names_to="indicator", values_to="value")
  
  # summarise the raw values
  summary_dsol <- long_dsol %>%
    group_by(time, scenario, indicator) %>%
    summarise(
      median = median(value),
      lower  = quantile(value, 0.025),
      upper  = quantile(value, 0.975),
      .groups = "drop"
    )

  
  return(list(draws = combined_dsol_all,
              summary_dsol = summary_dsol))
  
}




apply_counterfactuals_m2_1 <- function(posterior_samples, 
                                              n_sams,
                                              allocate_parameters,
                                              goveqs,
                                              goveqs_scaleup,
                                              bounds,
                                              prm, 
                                              ref, 
                                              thetai, 
                                              lhd_list,
                                              lhd_switch,
                                              pop_country){
  # Sample evenly spaced 250 parameters vectors from the posteriors after burnin and thinning
  total_samples <- nrow(posterior_samples)
  dx <- round(total_samples/n_sams)
  indices <- seq(1, total_samples, by = dx)
  indices <- indices[seq_len(min(n_sams, length(indices)))] # Ensure exactly 250 samples
  
  xs <- posterior_samples[indices, 1:nrow(bounds)]
  n_samples <- nrow(xs)
  dsol_list <- vector("list", length = n_samples)
  
  # Main loop over all selected posterior samples
  for (ii in seq_len(n_samples)) {
    # Print to see progress of intervention iterations
    if(ii %% 50 == 0) 
      print(paste("Intervention iteration: ", ii)) 
    
    theta <- as.numeric(xs[ii, ])
    
    # Allocate parameters based on input theta (parameter vector)
    alloc <- allocate_parameters(theta, thetai, prm$p, prm$r)
    r <- alloc$r
    p <- alloc$p
    i <- ref$i
    d <- ref$d
    
    # 1) Actual past trajectory (with RNTCP)
    init <- rep(0, i[["nx"]])
    names(init) <- c(
      "U", "Lf_ds", "Lf_rr", "Ls_ds", "Ls_rr", "I_ds", "I_rr",
      "Rlo_ds", "Rlo_rr", "Rhi_ds", "Rhi_rr", "R_ds", "R_rr",
      "Tx_ds_pu", "Tx_ds_pr_nen", "Tx_rr_pu", "Tx_rr_pr_nen",
      "Tx2_ds_pu", "Tx2_ds_pr_nen", "Tx2_rr_pu", "Tx2_rr_pr_nen",
      "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort", "cu_mort_ds", "cu_mort_rr", 
      "cu_noti_rr_pu")
    seed <- 1e-6
    init[i[["U"]]] <- 1 - seed; init[i[["I"]][["ds"]]] <- seed
    
    # equilibrium (no MDR, no RNTCP)
    # In absence of MDR and RNTCP
    p0 <- p; r0 <- r 
    p0[["pu"]] <- 0
    r0[["beta"]][2] <- 0
    r0[["MDR_acqu"]] <- 0
    r0[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters0 <- list(r = r0, p = p0)
    soln0 <- deSolve::ode(y = init, times = seq(0, 2000, by = 1000), func = goveqs, 
                          parms = parameters0,
                          method = "bdf")
    init1 <- tail(soln0[,-1],1)
    
    # MDR intro (no RNTCP)
    p1 <- p; r1 <- r 
    p1[["pu"]] <- 0
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters1 <- list(r = r1, p = p1)
    soln1 <- deSolve::ode(y = init1, times = seq(1970, 1997, by = 1), func = goveqs, 
                          parms = parameters1,
                          method = "bdf")
    init2 <- tail(soln1[,-1],1)
    
    # RNTCP scale-up
    p2 <- p; r2 <- r
    r2[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters2 <- list(r_start = r2, r_end = r2, p_start = p1, p_end = p2, times = c(1997, 2007))
    parameters2$goveqs <- goveqs
    soln2 <- deSolve::ode(y = init2, times = seq(1997, 2015, by = 1), func = goveqs_scaleup, 
                          parms = parameters2,
                          method = "bdf")
    init3 <- tail(soln2[,-1],1)
    
    # Treatment initiation scaleup
    r3_1 <- r
    r3_1[["Tx_init"]] <- r[["Tx_init_old"]]
    r3_2 <- r
    r3_2[["Tx_init"]] <- r[["Tx_init_new"]]
    p3 <- p
    parameters3 <- list(r_start = r3_1, r_end = r3_2, p_start = p3, p_end = p3, times = c(2015, 2025))
    parameters3$goveqs <- goveqs
    soln3 <- deSolve::ode(y = init3, times = seq(2015, 2025, by = 1), func = goveqs_scaleup, 
                          parms = parameters3,
                          method = "bdf")
    
    solution_df0 <- as.data.frame(soln0)
    colnames(solution_df0) <- c("time", d[1:i$nstates], 
                                "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr_pu")
    solution_df1 <- as.data.frame(soln1)
    colnames(solution_df1) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr_pu")
    solution_df1 <- solution_df1 |> dplyr::filter(time > 1970)
    solution_df2 <- as.data.frame(soln2)
    colnames(solution_df2) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr_pu")
    solution_df2 <- solution_df2 |> dplyr::filter(time > 1997)
    solution_df3 <- as.data.frame(soln3)
    colnames(solution_df3) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr_pu")
    solution_df3 <- solution_df3 |> dplyr::filter(time > 2015)
    
    sol_actual <- rbind(solution_df1, solution_df2, solution_df3)
    
    dsol_actual <- sol_actual %>%
      dplyr::select(time, cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr) %>%
      mutate(across(-time, ~ c(.x[1], diff(.x))),
             across(-time, ~ .x*1e5),
             scenario="Past actual") %>%
      rename(tb_mort = cu_mort_tb,
             inc_rr = cu_inc_rr,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    init_forward <- tail(soln3[,-1],1)  # forward intervention scenarios use this
    
    # 2) Counterfactual past (no p_MDR_rec, everything else same)
    seed <- 1e-6
    init[i[["U"]]] <- 1 - seed; init[i[["I"]][["ds"]]] <- seed
    
    # equilibrium (no MDR, no RNTCP)
    # In absence of MDR and RNTCP
    p0 <- p
    r0 <- r 
    p0[["pu"]] <- 0
    r0[["beta"]][2] <- 0
    r0[["MDR_acqu"]] <- 0
    r0[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters0 <- list(r = r0, p = p0)
    soln0 <- deSolve::ode(y = init, times = seq(0, 2000, by = 1000), func = goveqs, 
                          parms = parameters0,
                          method = "bdf")
    init1 <- tail(soln0[,-1],1)
    
    # MDR intro (no RNTCP)
    p1 <- p; r1 <- r 
    p1[["pu"]] <- 0
    p1[["MDR_rec"]][1] <- 0
    p1[["MDR_rec_pr_nen"]] <- 0
    r1[["muTx"]][2] <- r[["mort_TB"]]
    r1[["muTx"]][3] <- r[["mort_TB"]]
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters1 <- list(r = r1, p = p1)
    soln1 <- deSolve::ode(y = init1, times = seq(1970, 1997, by = 1), func = goveqs, 
                          parms = parameters1,
                          method = "bdf")
    init2 <- tail(soln1[,-1],1)
    
    # RNTCP scale-up
    p2 <- p; r2 <- r
    p2[["MDR_rec"]][1] <- 0
    p2[["MDR_rec_pr_nen"]] <- 0
    r2[["muTx"]][2] <- r[["mort_TB"]]
    r2[["muTx"]][3] <- r[["mort_TB"]]
    r2[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters2 <- list(r_start = r2, r_end = r2, p_start = p1, p_end = p2, times = c(1997, 2007))
    parameters2$goveqs <- goveqs
    soln2 <- deSolve::ode(y = init2, times = seq(1997, 2015, by = 1), func = goveqs_scaleup, 
                          parms = parameters2,
                          method = "bdf")
    init3 <- tail(soln2[,-1],1)
    
    # Treatment initiation scaleup
    r3_1 <- r
    r3_1[["Tx_init"]] <- r[["Tx_init_old"]]
    r3_2 <- r
    r3_2[["Tx_init"]] <- r[["Tx_init_new"]]
    p3 <- p
    p3[["MDR_rec"]][1] <- 0
    p3[["MDR_rec_pr_nen"]] <- 0
    r3_1[["muTx"]][2] <- r[["mort_TB"]]
    r3_2[["muTx"]][2] <- r[["mort_TB"]]
    r3_1[["muTx"]][3] <- r[["mort_TB"]]
    r3_2[["muTx"]][3] <- r[["mort_TB"]]
    parameters3 <- list(r_start = r3_1, r_end = r3_2, p_start = p3, p_end = p3, times = c(2015, 2025))
    parameters3$goveqs <- goveqs
    soln3 <- deSolve::ode(y = init3, times = seq(2015, 2025, by = 1), func = goveqs_scaleup, 
                          parms = parameters3,
                          method = "bdf")
    
    solution_df0 <- as.data.frame(soln0)
    colnames(solution_df0) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr_pu")
    solution_df1 <- as.data.frame(soln1)
    colnames(solution_df1) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr_pu")
    solution_df1 <- solution_df1 |> dplyr::filter(time > 1970)
    solution_df2 <- as.data.frame(soln2)
    colnames(solution_df2) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr_pu")
    solution_df2 <- solution_df2 |> dplyr::filter(time > 1997)
    solution_df3 <- as.data.frame(soln3)
    colnames(solution_df3) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr_pu")
    solution_df3 <- solution_df3 |> dplyr::filter(time > 2015)
    
    sol_cf <- rbind(solution_df1, solution_df2, solution_df3)
    
    dsol_cf <- sol_cf %>%
      as.data.frame() %>%
      dplyr::select(time, cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr) %>%
      mutate(across(-time, ~ c(.x[1], diff(.x))),
             across(-time, ~ .x*1e5),
             scenario="Past no p_MDR_rec")%>%
      rename(tb_mort = cu_mort_tb,
             inc_rr = cu_inc_rr,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    # 3) Forward projections using actual trajectory till 2025
    # add 0s for capturing treatment starters
    init <- c(init_forward)
    
    # scenario 0 - base model
    r0 <- r
    p0 <- p
    r0[["Tx_init"]] <- r[["Tx_init_new"]]
    parameters0 <- list(r = r0, p = p0)
    soln0_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs, 
                            parms = parameters0,
                            method = "bdf")
    
    # scenario 1 - improve DST as final condition in scaleup in public sector (and increasing private sector engagement)
    p1 <- p0
    r1 <- r0
    r1[["Tx_init"]] <- r[["Tx_init_new"]]
    # Increasing p_pu to 0.9 for DST scaleup in public and increased private sector engagement
    p1[["pu"]] <- max(p1[["pu"]], 0.9)
    p1$MDR_rec <- max(p1$MDR_rec, 0.75) 
    parameters1 <- list(r_start = r0, r_end = r1, p_start = p0, p_end = p1, times = c(2026, 2028))
    parameters1$goveqs <- goveqs
    soln1_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs_scaleup, 
                            parms = parameters1,
                            method = "bdf")
    
    
    # scenario 2 - improve SL treatment as final condition in public sector (and increasing private sector engagement)
    # a) shift 50% of bad outcomes (LTFU and mortality) into success
    # b) halve LTFU when on Tx2
    # c) halve mortality when on Tx2
    p2 <- p1; r2 <- r1
    vec <- c(r2$Tx2, r2$ltfu[2], r2$muTx[2])
    propvec <- vec/sum(vec) # current proportions of Tx2 outcomes
    newvec <- c(propvec[1] + 0.5 * (1 - propvec[1]), 0.5 * propvec[2], 0.5 * propvec[3]) # new props of Tx2 outcomes
    r2$ltfu[2] <- r2$Tx2 * newvec[2] / newvec[1] # scale lower LTFU rate with new props
    r2$muTx[2] <- r2$Tx2 * newvec[3] / newvec[1] # scale lower muTx2 rate with new props
    p2$cure[2] <- (p2$cure[2] + 1)*0.5 # reducing failure by 50% i.e increasing success. r_Tx2 remains the same
    
    parameters2 <- list(r_start = r0, r_end = r2, p_start = p0, p_end = p2, times = c(2026, 2028))
    parameters2$goveqs <- goveqs
    soln2_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs_scaleup, 
                            parms = parameters2,
                            method = "bdf")
    
    # 4) Combine outputs
    soln0_f <- as.data.frame(soln0_f) %>% mutate(scenario = "Baseline")
    soln1_f <- as.data.frame(soln1_f) %>% mutate(scenario = "Improved DST")
    soln2_f <- as.data.frame(soln2_f) %>% mutate(scenario = "Improved DST + SL outcomes")
    
    soln_f <- bind_rows(soln0_f, soln1_f, soln2_f)
    colnames(soln_f) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", 
                          "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr", 
                          "cu_noti_rr_pu",
                          "scenario")
    
    dsol_f <- soln_f %>%
      dplyr::select("time", "scenario", "cu_inc_rr", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr") %>%
      group_by(scenario) %>%
      mutate(across(c(cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr),
                    ~ c(.x[1], diff(.x))),
             across(-time, ~ .x * 1e5)) %>% # to get per 100k per year
      rename(inc_rr = cu_inc_rr,
             tb_mort = cu_mort_tb,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    dsol <- bind_rows(dsol_actual, dsol_cf, dsol_f)
    dsol$iter <- ii
    dsol_list[[ii]] <- dsol
  }
  
  combined_dsol_all <- bind_rows(dsol_list)
  
  # Summarise inc_rr, tb_mort, tb_mort_ds, tb_mort_rr in the past and future scenarios
  # pivot long
  long_dsol <- combined_dsol_all %>%
    pivot_longer(cols=c(inc_rr, tb_mort, tb_mort_ds, tb_mort_rr),
                 names_to="indicator", values_to="value")
  
  # summarise the raw values
  summary_dsol <- long_dsol %>%
    group_by(time, scenario, indicator) %>%
    summarise(
      median = median(value),
      lower  = quantile(value, 0.025),
      upper  = quantile(value, 0.975),
      .groups = "drop"
    )

  
  return(list(draws = combined_dsol_all,
              summary_dsol = summary_dsol))
  
}


apply_counterfactuals_m3_1 <- function(posterior_samples, 
                                       n_sams,
                                       allocate_parameters, 
                                       goveqs, 
                                       goveqs_scaleup,
                                       bounds, 
                                       prm, 
                                       ref, 
                                       thetai, 
                                       pop_country) {
  
  # Select evenly spaced posterior samples
  total_samples <- nrow(posterior_samples)
  dx <- max(1, round(total_samples/n_sams))
  
  indices <- seq(1, total_samples, by = dx)
  indices <- indices[seq_len(min(n_sams, length(indices)))]
  
  xs <- posterior_samples[indices, 1:nrow(bounds), drop = FALSE]
  # xs$sim_rr_inc24 <- NA
  n_samples <- nrow(xs)
  dsol_list <- vector("list", length = n_samples)
  
  # Loop over posterior samples
  for (ii in seq_len(n_samples)) {
    
    if (ii %% 50 == 0) {
      print(paste("Intervention iteration:", ii))
    }
    
    theta <- as.numeric(xs[ii, ])
    
    # Allocate parameters
    alloc <- allocate_parameters(theta, thetai, prm$p, prm$r)
    r <- alloc$r
    p <- alloc$p
    i <- ref$i
    d <- ref$d
    
    # 1. ACTUAL HISTORICAL TRAJECTORY
    init <- rep(0, i[["nx"]])
    names(init) <- c("U",
                     "Lf_ds", "Lf_rr",
                     "Ls_ds", "Ls_rr",
                     "I_ds", "I_rr",
                     "Rlo_ds", "Rlo_rr",
                     "Rhi_ds", "Rhi_rr",
                     "R_ds", "R_rr",
                     "Tx_ds", "Tx_rr",
                     "Tx2_ds", "Tx2_rr",
                     "cu_inc_ds", "cu_inc_rr",
                     "cu_noti",
                     "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr",
                     "cu_noti_rr")
    
    seed <- 1e-6
    init[i[["U"]]] <- 1 - seed
    init[i[["I"]][["ds"]]] <- seed
    
    # Equilibrium: no MDR
    p0 <- p
    r0 <- r
    
    r0[["beta"]][2] <- 0
    r0[["MDR_acqu"]] <- 0
    r0[["Tx_init"]] <- r[["Tx_init_old"]]
    
    parameters0 <- list(r = r0, 
                        p = p0) 
    
    soln0 <- deSolve::ode(y = init,
                          times = seq(0, 2000, by = 1000),
                          func = goveqs,
                          parms = parameters0,
                          method = "bdf")
    
    init1 <- tail(soln0[, -1], 1)
    
    # MDR introduced, but no MDR recognition yet
    # 1970-1997
    p1 <- p
    r1 <- r
    
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    
    parameters1 <- list(r = r1, 
                        p = p1)
    
    soln1 <- deSolve::ode(y = init1,
                          times = seq(1970, 1997, by = 1),
                          func = goveqs,
                          parms = parameters1,
                          method = "bdf")
    
    init2 <- tail(soln1[, -1], 1)
    
    # MDR recognition/DST scale-up
    # 1997-2007
    p2 <- p
    r2 <- r
    
    r2[["Tx_init"]] <- r[["Tx_init_old"]]
    
    parameters2 <- list(r_start = r2,
                        r_end   = r2,
                        p_start = p1,
                        p_end   = p2,
                        times   = c(1997, 2007))
    parameters2$goveqs <- goveqs
    
    soln2 <- deSolve::ode(y = init2,
                          times = seq(1997, 2015, by = 1),
                          func = goveqs_scaleup,
                          parms = parameters2,
                          method = "bdf")
    
    init3 <- tail(soln2[, -1], 1)
    
    # Treatment initiation scale-up
    # 2015-2025
    r3_start <- r
    r3_start[["Tx_init"]] <- r[["Tx_init_old"]]
    r3_end <- r
    r3_end[["Tx_init"]] <- r[["Tx_init_new"]]
    p3 <- p
    
    parameters3 <- list(r_start = r3_start,
                        r_end   = r3_end,
                        p_start = p3,
                        p_end   = p3,
                        times   = c(2015, 2025))
    parameters3$goveqs <- goveqs
    
    soln3 <- deSolve::ode(y = init3,
                          times = seq(2015, 2025, by = 1),
                          func = goveqs_scaleup,
                          parms = parameters3,
                          method = "bdf")
    
    # Convert historical solutions
    solution_names <- c("time",
                        d[1:i$nstates],
                        "cu_inc_ds", "cu_inc_rr",
                        "cu_noti",
                        "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr",
                        "cu_noti_rr")
    
    solution_df1 <- as.data.frame(soln1)
    colnames(solution_df1) <- solution_names
    solution_df1 <- solution_df1 %>%
      dplyr::filter(time > 1970)
    
    solution_df2 <- as.data.frame(soln2)
    colnames(solution_df2) <- solution_names
    solution_df2 <- solution_df2 %>%
      dplyr::filter(time > 1997)
    
    solution_df3 <- as.data.frame(soln3)
    colnames(solution_df3) <- solution_names
    solution_df3 <- solution_df3 %>%
      dplyr::filter(time > 2015)
    
    sol_actual <- dplyr::bind_rows(solution_df1,
                                   solution_df2,
                                   solution_df3)
    
    dsol_actual <- sol_actual %>%
      dplyr::select(time,
                    cu_inc_rr,
                    cu_mort_tb,
                    cu_mort_tb_ds,
                    cu_mort_tb_rr) %>%
      dplyr::mutate(dplyr::across(-time, ~ c(.x[1], diff(.x))),
                    dplyr::across(-time, ~ .x * 1e5), scenario = "Past actual") %>%
      dplyr::rename(inc_rr     = cu_inc_rr,
                    tb_mort    = cu_mort_tb,
                    tb_mort_ds = cu_mort_tb_ds,
                    tb_mort_rr = cu_mort_tb_rr)
    
    
    # Actual 2025 state used for ALL future interventions
    init_forward <- tail(soln3[, -1], 1)
    
    # 2. HISTORICAL COUNTERFACTUAL
    # MDR exists, but MDR_rec remains zero
    init_cf <- init
    
    # Equilibrium
    p0_cf <- p
    r0_cf <- r
    
    r0_cf[["beta"]][2] <- 0
    r0_cf[["MDR_acqu"]] <- 0
    r0_cf[["Tx_init"]] <- r[["Tx_init_old"]]
    
    parameters0_cf <- list(r = r0_cf,
                           p = p0_cf)
    
    soln0_cf <- deSolve::ode(y = init_cf,
                             times = seq(0, 2000, by = 1000),
                             func = goveqs,
                             parms = parameters0_cf,
                             method = "bdf")
    
    init1_cf <- tail(soln0_cf[, -1], 1)
    
    # MDR introduced but NEVER recognised in the public sector
    # 1970-1997
    p1_cf <- p
    r1_cf <- r
    
    p1_cf[["MDR_rec"]] <- 0
    r1_cf[["Tx_init"]] <- r[["Tx_init_old"]]
    r1_cf[["muTx"]][2] <- r[["mort_TB"]] # RR patients don't receive SL treatment mortality benefit
    
    parameters1_cf <- list(r = r1_cf,
                           p = p1_cf)
    
    soln1_cf <- deSolve::ode(y = init1_cf,
                             times = seq(1970, 1997, by = 1),
                             func = goveqs,
                             parms = parameters1_cf,
                             method = "bdf")
    
    init2_cf <- tail(soln1_cf[, -1], 1)
    
    # 1997-2015: q scales up, MDR recognition remains zero
    p2_cf <- p
    p2_cf[["MDR_rec"]]<- 0
    
    r2_cf <- r
    r2_cf[["Tx_init"]] <- r[["Tx_init_old"]]
    r2_cf[["muTx"]][2] <- r[["mort_TB"]] # RR patients don't receive SL treatment mortality benefit
    
    parameters2_cf <- list(r_start = r2_cf,
                           r_end   = r2_cf,
                           p_start = p1_cf,
                           p_end   = p2_cf,
                           times   = c(1997, 2007),
                           goveqs  = goveqs)
    parameters2_cf$goveqs <- goveqs
    
    soln2_cf <- deSolve::ode(y = init2_cf,
                             times = seq(1997, 2015, by = 1),
                             func = goveqs_scaleup,
                             parms = parameters2_cf,
                             method = "bdf")
    
    init3_cf <- tail(soln2_cf[, -1], 1)
    
    # 2015-2025: treatment initiation scales up
    # MDR recognition remains zero
    r3_start_cf <- r
    r3_start_cf[["Tx_init"]] <- r[["Tx_init_old"]]
    r3_start_cf[["muTx"]][2] <- r[["mort_TB"]]
    
    r3_end_cf <- r
    r3_end_cf[["Tx_init"]] <- r[["Tx_init_new"]]
    r3_end_cf[["muTx"]][2] <- r[["mort_TB"]]
    
    p3_cf <- p
    p3_cf[["MDR_rec"]] <- 0
    
    parameters3_cf <- list(r_start = r3_start_cf,
                           r_end   = r3_end_cf,
                           p_start = p3_cf,
                           p_end   = p3_cf,
                           times   = c(2015, 2025),
                           goveqs  = goveqs)
    parameters3_cf$goveqs <- goveqs
    
    soln3_cf <- deSolve::ode(y = init3_cf,
                             times = seq(2015, 2025, by = 1),
                             func = goveqs_scaleup,
                             parms = parameters3_cf,
                             method = "bdf")
    
    # Convert counterfactual outputs
    
    solution_df1_cf <- as.data.frame(soln1_cf)
    colnames(solution_df1_cf) <- solution_names
    solution_df1_cf <- solution_df1_cf %>%
      dplyr::filter(time > 1970)
    
    solution_df2_cf <- as.data.frame(soln2_cf)
    colnames(solution_df2_cf) <- solution_names
    solution_df2_cf <- solution_df2_cf %>%
      dplyr::filter(time > 1997)
    
    solution_df3_cf <- as.data.frame(soln3_cf)
    colnames(solution_df3_cf) <- solution_names
    solution_df3_cf <- solution_df3_cf %>%
      dplyr::filter(time > 2015)
    
    sol_cf <- dplyr::bind_rows(solution_df1_cf, solution_df2_cf, solution_df3_cf)
    
    dsol_cf <- sol_cf %>%
      dplyr::select(time,
                    cu_inc_rr,
                    cu_mort_tb,
                    cu_mort_tb_ds,
                    cu_mort_tb_rr) %>%
      dplyr::mutate(dplyr::across(-time, ~ c(.x[1], diff(.x))),
                    dplyr::across(-time, ~ .x * 1e5), scenario = "Past no MDR recognition") %>%
      dplyr::rename(inc_rr     = cu_inc_rr,
                    tb_mort    = cu_mort_tb,
                    tb_mort_ds = cu_mort_tb_ds,
                    tb_mort_rr = cu_mort_tb_rr)
    
    # 3. FUTURE INTERVENTIONS
    #    All start from ACTUAL 2025 state
    init_f <- init_forward
    
    # Scenario 0: Baseline
    p0_f <- p
    r0_f <- r
    
    r0_f[["Tx_init"]] <- r[["Tx_init_new"]]
    
    parameters0_f <- list(r = r0_f,
                          p = p0_f)
    
    soln0_f <- deSolve::ode(y = as.numeric(init_f),
                            times = seq(2025, 2036, by = 1),
                            func = goveqs,
                            parms = parameters0_f,
                            method = "bdf")
    
    # Scenario 1: Improved DST/MDR recognition in public sector
    p1_f <- p0_f
    r1_f <- r0_f
    
    p1_f[["MDR_rec"]] <- max(p1_f[["MDR_rec"]], 0.75)
    
    parameters1_f <- list(r_start = r0_f,
                          r_end   = r1_f,
                          p_start = p0_f,
                          p_end   = p1_f,
                          times   = c(2026, 2028))
    
    parameters1_f$goveqs <- goveqs
    
    soln1_f <- deSolve::ode(y = as.numeric(init_f),
                            times = seq(2025, 2036, by = 1),
                            func = goveqs_scaleup,
                            parms = parameters1_f,
                            method = "bdf")
    
    # Scenario 2: Improved DST + SL treatment outcomes
    p2_f <- p1_f
    r2_f <- r1_f
    
    vec <- c(r2_f$Tx2, r2_f$ltfu[2], r2_f$muTx[2])
    propvec <- vec / sum(vec)
    newvec <- c(propvec[1] + 0.5 * (1 - propvec[1]),
                0.5 * propvec[2],
                0.5 * propvec[3])
    r2_f$ltfu[2] <- r2_f$Tx2*newvec[2]/newvec[1]
    r2_f$muTx[2] <- r2_f$Tx2*newvec[3]/newvec[1]
    p2_f$cure[2] <- (p2_f$cure[2] + 1)*0.5
    
    parameters2_f <- list(r_start = r0_f,
                          r_end   = r2_f,
                          p_start = p0_f,
                          p_end   = p2_f,
                          times   = c(2026, 2028))
    
    parameters2_f$goveqs <- goveqs
    
    soln2_f <- deSolve::ode(y = as.numeric(init_f),
                            times = seq(2025, 2036, by = 1),
                            func = goveqs_scaleup,
                            parms = parameters2_f,
                            method = "bdf")
    
    # 4. Process future outputs
    soln0_f <- as.data.frame(soln0_f) %>%
      dplyr::mutate(scenario = "Baseline")
    
    soln1_f <- as.data.frame(soln1_f) %>%
      dplyr::mutate(scenario = "Improved DST")
    
    soln2_f <- as.data.frame(soln2_f) %>%
      dplyr::mutate(scenario = "Improved DST + SL outcomes")
    
    soln_f <- dplyr::bind_rows(soln0_f, soln1_f, soln2_f)
    colnames(soln_f) <- c(solution_names, "scenario")
    
    
    dsol_f <- soln_f %>%
      dplyr::select(time,
                    scenario,
                    cu_inc_rr,
                    cu_mort_tb,
                    cu_mort_tb_ds,
                    cu_mort_tb_rr) %>%
      dplyr::group_by(scenario) %>%
      mutate(across(c(cu_inc_rr, 
                      cu_mort_tb, 
                      cu_mort_tb_ds, 
                      cu_mort_tb_rr),
                    ~ c(.x[1], diff(.x))),
             across(-time, ~ .x * 1e5)) %>% # to get per 100k per year
      rename(inc_rr       = cu_inc_rr,
             tb_mort      = cu_mort_tb,
             tb_mort_ds   = cu_mort_tb_ds,
             tb_mort_rr   = cu_mort_tb_rr)
    
    # Combine
    dsol <- dplyr::bind_rows(dsol_actual, dsol_cf, dsol_f)
    dsol$iter <- ii
    dsol_list[[ii]] <- dsol
  }
  
  # Combine posterior draws
  combined_dsol_all <- dplyr::bind_rows(dsol_list)
  
  # Incidence / mortality summaries
  long_dsol <- combined_dsol_all %>%
    tidyr::pivot_longer(cols = c(inc_rr,
                                 tb_mort,
                                 tb_mort_ds,
                                 tb_mort_rr),
                        names_to = "indicator",
                        values_to = "value")
  
  
  summary_dsol <- long_dsol %>%
    dplyr::group_by(time, scenario, indicator) %>%
    dplyr::summarise(median = median(value, na.rm = TRUE),
                     lower  = quantile(value, 0.025, na.rm = TRUE),
                     upper  = quantile(value, 0.975, na.rm = TRUE),
                     .groups = "drop")
  
  
  return(list(draws = combined_dsol_all,
              summary_dsol = summary_dsol))
}

apply_counterfactuals_m4_1 <- function(posterior_samples, 
                                       n_sams,
                                       allocate_parameters, 
                                       goveqs, 
                                       goveqs_scaleup,
                                       bounds, 
                                       prm, 
                                       ref, 
                                       thetai, 
                                       lhd_list,
                                       lhd_switch,
                                       pop_country){
  # Sample evenly spaced 250 parameters vectors from the posteriors after burnin and thinning
  total_samples <- nrow(posterior_samples)
  dx <- round(total_samples/n_sams)
  indices <- seq(1, total_samples, by = dx)
  indices <- indices[seq_len(min(n_sams, length(indices)))] # Ensure exactly 250 samples
  
  xs <- posterior_samples[indices, 1:nrow(bounds)]
  xs$sim_rr_inc24 <- NA
  
  # # Pre-allocate storage for results
  n_samples <- nrow(xs)
  dsol_list <- vector("list", length = n_samples)
  
  # Main loop over all selected posterior samples
  for (ii in seq_len(n_samples)) {
    # Print to see progress of intervention iterations
    if(ii %% 50 == 0) 
      print(paste("Intervention iteration: ", ii)) 
    
    theta <- as.numeric(xs[ii, ])
    
    # Allocate parameters based on input theta (parameter vector)
    alloc <- allocate_parameters(theta, thetai, prm$p, prm$r)
    r <- alloc$r
    p <- alloc$p
    i <- ref$i
    s <- ref$s
    d <- ref$d
    
    # 1) Actual past trajectory
    init <- rep(0, i[["nx"]])
    names(init) <- c(
      d[1:i[["nstates"]]],
      "cu_inc_ds",
      "cu_inc_rr",
      "cu_noti",
      "cu_mort_tb",
      "cu_mort_tb_ds",
      "cu_mort_tb_rr",
      "cu_noti_rr",
      "cu_inc_h1hart",
      "cu_mort_tb_h0",
      "cu_mort_tb_h1hart"
    )
    seed <- 1e-6
    init[i[["U"]][["h0"]]] <- 1 - seed; init[i[["I"]][["ds"]][["h0"]]] <- seed
    
    # In absence of MDR and no ART
    p0 <- p; r0 <- r
    r0[["beta"]][2] <- 0
    r0[["MDR_acqu"]] <- 0
    r0[["Tx_init"]] <- r[["Tx_init_old"]]
    r0[["ART_init"]] <- 0
    parameters0 <- list(r = r0, p = p0, i = i)
    soln0 <- deSolve::ode(y = init, times = seq(0, 2000, by = 1000), func = goveqs, 
                          parms = parameters0,
                          method = "bdf")  
    soln0[, -1] <- pmax(soln0[, -1], 0)
    
    # 2. Introduce MDR-TB and no ART
    init1 <- tail(soln0[, -1], 1)
    if (is.null(init1)) return(-1e10)
    p1 <- p; r1 <- r
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    r1[["ART_init"]] <- 0
    parameters1 <- list(r = r1, p = p1, i = i)
    soln1 <- deSolve::ode(y = init1, times = seq(1970, 2000, by = 1), func = goveqs, 
                          parms = parameters1,
                          method = "bdf")
    soln1[, -1] <- pmax(soln1[, -1], 0)
    
    # 3. Treatment initiation scale-up and ART scaleup
    init2 <- tail(soln1[, -1], 1)
    if (is.null(init2)) return(-1e10)
    r2_1 <- r
    r2_1[["Tx_init"]] <- r[["Tx_init_old"]]
    r2_2 <- r
    r2_2[["Tx_init"]] <- r[["Tx_init_new"]]
    p2 <- p
    r2_1[["ART_init"]] <- 0
    r2_2[["ART_init"]] <- r[["ART_init"]]
    parameters2 <- list(r_start = r2_1, r_end = r2_2, p_start = p2, p_end = p2, times = c(2000, 2025), i = i)
    parameters2$goveqs <- goveqs
    soln2 <- deSolve::ode(y = init2, times = seq(2000, 2025, by = 1), func = goveqs_scaleup, 
                          parms = parameters2,
                          method = "bdf")
    soln2[, -1] <- pmax(soln2[, -1], 0)
    
    # Calculate outputs
    solution_df0 <- as.data.frame(soln0)
    colnames(solution_df0) <- c("time", d[1:i$nstates], 
                                "cu_inc_ds", "cu_inc_rr", 
                                "cu_noti", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr",
                                "cu_inc_h1hart", "cu_mort_tb_h0", "cu_mort_tb_h1hart")
    solution_df1 <- as.data.frame(soln1)
    colnames(solution_df1) <- c("time", d[1:i$nstates], 
                                "cu_inc_ds", "cu_inc_rr", 
                                "cu_noti", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr",
                                "cu_inc_h1hart", "cu_mort_tb_h0", "cu_mort_tb_h1hart")
    solution_df1 <- solution_df1 |> dplyr::filter(time > 1970)
    solution_df2 <- as.data.frame(soln2)
    colnames(solution_df2) <- c("time", d[1:i$nstates], 
                                "cu_inc_ds", "cu_inc_rr", 
                                "cu_noti", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr",
                                "cu_inc_h1hart", "cu_mort_tb_h0", "cu_mort_tb_h1hart")
    solution_df2 <- solution_df2 |> dplyr::filter(time > 2000)
    
    sol_actual <- rbind(solution_df1, solution_df2)
    dsol_actual <- sol_actual %>%
      dplyr::select(time, cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr) %>%
      mutate(across(-time, ~ c(.x[1], diff(.x))),
             across(-time, ~ .x*1e5),
             scenario="Past actual") %>%
      rename(tb_mort = cu_mort_tb,
             inc_rr = cu_inc_rr,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    init_forward <- tail(soln2[,-1],1)  # forward intervention scenarios use this
    
    # 2) Counterfactual past (no p_MDR_rec, everything else same)
    seed <- 1e-6
    init[i[["U"]][["h0"]]] <- 1 - seed; init[i[["I"]][["ds"]][["h0"]]] <- seed
    
    # In absence of MDR and no ART
    p0 <- p; r0 <- r
    r0[["beta"]][2] <- 0
    r0[["MDR_acqu"]] <- 0
    r0[["Tx_init"]] <- r[["Tx_init_old"]]
    r0[["ART_init"]] <- 0
    parameters0 <- list(r = r0, p = p0, i = i)
    soln0 <- deSolve::ode(y = init, times = seq(0, 2000, by = 1000), func = goveqs, 
                          parms = parameters0,
                          method = "bdf")  
    soln0[, -1] <- pmax(soln0[, -1], 0)
    
    # 2. Introduce MDR-TB and no ART, but not RR services
    init1 <- tail(soln0[, -1], 1)
    if (is.null(init1)) return(-1e10)
    p1 <- p
    r1 <- r
    p1[["MDR_rec"]] <- 0
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    r1[["ART_init"]] <- 0
    r1["muTx"][2] <- r[["mort_TB"]] # SL for HIV- and HIV+ and ART people
    parameters1 <- list(r = r1, p = p1, i = i)
    soln1 <- deSolve::ode(y = init1, times = seq(1970, 2000, by = 1), func = goveqs, 
                          parms = parameters1,
                          method = "bdf")
    soln1[, -1] <- pmax(soln1[, -1], 0)
    
    
    # 3. Treatment initiation scale-up and ART scaleup, but not RR services
    init2 <- tail(soln1[, -1], 1)
    r2_1 <- r
    r2_1[["Tx_init"]] <- r[["Tx_init_old"]]
    r2_2 <- r
    r2_2[["Tx_init"]] <- r[["Tx_init_new"]]
    p2 <- p
    p2[["MDR_rec"]] <- 0
    r2_1[["ART_init"]] <- 0
    r2_2[["ART_init"]] <- r[["ART_init"]]
    r2_1["muTx"][2] <- r[["mort_TB"]]
    r2_2["muTx"][2] <- r[["mort_TB"]]
    parameters2 <- list(r_start = r2_1, r_end = r2_2, p_start = p2, p_end = p2, times = c(2000, 2025), i = i)
    parameters2$goveqs <- goveqs
    soln2 <- deSolve::ode(y = init2, times = seq(2000, 2025, by = 1), func = goveqs_scaleup, 
                          parms = parameters2,
                          method = "bdf")
    soln2[, -1] <- pmax(soln2[, -1], 0)
    
    # Calculate outputs
    solution_df0 <- as.data.frame(soln0)
    colnames(solution_df0) <- c("time", d[1:i$nstates], 
                                "cu_inc_ds", "cu_inc_rr", 
                                "cu_noti", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr",
                                "cu_inc_h1hart", "cu_mort_tb_h0", "cu_mort_tb_h1hart")
    solution_df1 <- as.data.frame(soln1)
    colnames(solution_df1) <- c("time", d[1:i$nstates], 
                                "cu_inc_ds", "cu_inc_rr", 
                                "cu_noti", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr",
                                "cu_inc_h1hart", "cu_mort_tb_h0", "cu_mort_tb_h1hart")
    solution_df1 <- solution_df1 |> dplyr::filter(time > 1970)
    solution_df2 <- as.data.frame(soln2)
    colnames(solution_df2) <- c("time", d[1:i$nstates], 
                                "cu_inc_ds", "cu_inc_rr", 
                                "cu_noti", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr",
                                "cu_inc_h1hart", "cu_mort_tb_h0", "cu_mort_tb_h1hart")
    solution_df2 <- solution_df2 |> dplyr::filter(time > 2000)
    
    sol_cf <- rbind(solution_df1, solution_df2)
    dsol_cf <- sol_cf %>%
      as.data.frame() %>%
      dplyr::select(time, cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr) %>%
      mutate(across(-time, ~ c(.x[1], diff(.x))),
             across(-time, ~ .x*1e5),
             scenario="Past no p_MDR_rec")%>%
      rename(tb_mort = cu_mort_tb,
             inc_rr = cu_inc_rr,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    xs$sim_rr_inc24[ii] <- dsol_cf$inc_rr[nrow(dsol_cf)]
    
    # 3) Forward projections using actual trajectory till 2025
    init <- init_forward
    
    # scenario 0 - base model
    r0 <- r
    p0 <- p
    r0[["Tx_init"]] <- r0[["Tx_init_new"]]
    parameters0 <- list(r = r0, p = p0, i = i)
    soln0_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs, 
                            parms = parameters0,
                            method = "bdf")
    soln0_f[, -1] <- pmax(soln0_f[, -1], 0)
    
    # scenario 1 - improve DST as final condition in scaleup
    p1 <- p0
    r1 <- r0
    r1[["Tx_init"]] <- r1[["Tx_init_new"]]
    p1[["MDR_rec"]] <- max(p1[["MDR_rec"]], 0.75)
    parameters1 <- list(r_start = r0, r_end = r1, p_start = p0, p_end = p1, times = c(2026, 2028), i = i)
    parameters1$goveqs <- goveqs
    soln1_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs_scaleup, 
                            parms = parameters1,
                            method = "bdf")
    soln1_f[, -1] <- pmax(soln1_f[, -1], 0)
    
    # scenario 2 - improve SL treatment as final condition in scaleup
    # a) shift 50% of bad outcomes (LTFU and mortality) into success
    # b) halve LTFU when on Tx2
    # c) halve mortality when on Tx2
    p2 <- p1; r2 <- r1
    vec <- c(r2$Tx2, r2$ltfu[2], r2$muTx[2])
    propvec <- vec/sum(vec) # current proportions of Tx2 outcomes
    newvec <- c(propvec[1] + 0.5 * (1 - propvec[1]), 0.5 * propvec[2], 0.5 * propvec[3]) # new props of Tx2 outcomes
    # newvec <- c(propvec[1] + 0.99 * (1 - propvec[1]), 0.01 * propvec[2], 0.01 * propvec[3]) # new props of Tx2 outcomes
    r2$ltfu[2] <- r2$Tx2 * newvec[2] / newvec[1] # scale lower LTFU rate with new props
    r2$muTx[2] <- r2$Tx2 * newvec[3] / newvec[1] # scale lower muTx2 rate with new props
    # p2$cure[2] <- (p2$cure[2] + 1)*0.99 # reducing failure by 50% i.e increasing success. r_Tx2 remains the same
    p2$cure[2] <- (p2$cure[2] + 1)*0.5 # reducing failure by 50% i.e increasing success. r_Tx2 remains the same
    parameters2 <- list(r_start = r0, r_end = r2, p_start = p0, p_end = p2, times = c(2026, 2028), i = i)
    parameters2$goveqs <- goveqs
    soln2_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs_scaleup, 
                            parms = parameters2,
                            method = "bdf")
    soln2_f[, -1] <- pmax(soln2_f[, -1], 0)
    
    # 4) Combine outputs
    soln0_f <- as.data.frame(soln0_f) %>% mutate(scenario = "Baseline")
    soln1_f <- as.data.frame(soln1_f) %>% mutate(scenario = "Improved DST")
    soln2_f <- as.data.frame(soln2_f) %>% mutate(scenario = "Improved DST + SL outcomes")
    
    soln_f <- bind_rows(soln0_f, soln1_f, soln2_f)
    colnames(soln_f) <- c("time", d[1:i$nstates], 
                          "cu_inc_ds", "cu_inc_rr", 
                          "cu_noti", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_rr",
                          "cu_inc_h1hart", "cu_mort_tb_h0", "cu_mort_tb_h1hart",
                          "scenario")
    dsol_f <- soln_f %>%
      dplyr::select("time", "scenario", 
                    "cu_inc_rr", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr") %>%
      group_by(scenario) %>%
      mutate(across(c(cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr),
                    ~ c(.x[1], diff(.x))),
             across(-time, ~ .x * 1e5)) %>% # to get per 100k per year
      rename(inc_rr = cu_inc_rr,
             tb_mort = cu_mort_tb,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    dsol <- bind_rows(dsol_actual, dsol_cf, dsol_f)
    dsol$iter <- ii
    dsol_list[[ii]] <- dsol
  }
  
  combined_dsol_all <- bind_rows(dsol_list)
  
  # pivot long
  long_dsol <- combined_dsol_all %>%
    pivot_longer(cols=c(inc_rr,tb_mort, tb_mort_ds, tb_mort_rr),
                 names_to="indicator", values_to="value")
  
  # summarise the raw values
  summary_dsol <- long_dsol %>%
    group_by(time, scenario, indicator) %>%
    summarise(
      median = median(value),
      lower  = quantile(value, 0.025),
      upper  = quantile(value, 0.975),
      .groups = "drop"
    )
  
  return(list(draws = combined_dsol_all,
              summary_dsol = summary_dsol))
  
}




plot_counterfactuals_1 <- function(counterfactuals_summary) {
  # Past data
  past_data <- counterfactuals_summary %>%
    dplyr::filter(scenario %in% c("Past actual",
                                  "Past no MDR recognition")) %>%
    dplyr::filter(time > 1980) %>%
    dplyr::mutate(scenario = dplyr::recode(scenario,
                                           "Past actual" = "Actual",
                                           "Past no MDR recognition" = "Counterfactual: no RR services"))
  
  # Future data
  future_data <- counterfactuals_summary %>%
    dplyr::filter(scenario %in% c("Baseline",
                                  "Improved DST",
                                  "Improved DST + SL outcomes")) %>%
    dplyr::filter(time > 2025) %>%
    dplyr::mutate(time_year = time - 1)
  
  # Past counterfactual plot
  p_past <- ggplot(past_data %>%
                     dplyr::filter(indicator %in% c("inc_rr", "tb_mort_rr")),
                   aes(x = time, y = median, color = scenario, fill = scenario)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.2, color = NA) +
    facet_wrap(~indicator, scales = "free_y", 
               labeller = labeller(indicator = c("inc_rr" = "RR-TB incidence",
                                                 "tb_mort_rr" = "RR-TB mortality"))) +
    expand_limits(y = 0) +
    scale_color_manual(values = c("Actual" = "#DDA448",
                                  "Counterfactual: no RR services" = "#ef476f")) +
    scale_fill_manual(values = c("Actual" = "#DDA448", 
                                 "Counterfactual: no RR services" = "#ef476f")) +
    labs(title = "Historical counterfactual",
         y = "Rate per 100k population",
         x = "Year",
         fill = "Scenario",
         color = "Scenario") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          axis.text.x = element_text(angle = 25, hjust = 0.7))
  
  
  # Future intervention plot
  p_future <- ggplot(future_data %>%
                       dplyr::filter(indicator %in% c("inc_rr", "tb_mort_rr")),
                     aes(x = time_year,
                         y = median,
                         color = scenario,
                         fill = scenario)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.2, color = NA) +
    facet_wrap(~indicator, 
               scales = "free_y",
               labeller = ggplot2::labeller(indicator = c("inc_rr" = "RR-TB incidence",
                                                          "tb_mort_rr" = "RR-TB mortality"))) +
    expand_limits(y = 0) +
    scale_x_continuous(breaks = seq(min(future_data$time_year),
                                    max(future_data$time_year),
                                    by = 2)) +
    scale_color_manual(values = c("Baseline" = "#DDA448",
                                  "Improved DST" = "#785589",
                                  "Improved DST + SL outcomes" = "#14591D")) +
    scale_fill_manual(values = c("Baseline" = "#DDA448",
                                 "Improved DST" = "#785589",
                                 "Improved DST + SL outcomes" = "#14591D")) +
    labs(title = "Model projections over time",
         y = "Rate per 100k population",
         x = "Year",
         fill = "Scenario",
         color = "Scenario") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          axis.text.x = element_text(angle = 25, hjust = 0.7)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE),
           color = guide_legend(nrow = 2, byrow = TRUE))
  
  
  return(
    list(
      p_past = p_past,
      p_future = p_future
    )
  )
}



# For India sensitivity of not increasing engagement
apply_counterfactuals_m1_1_sens <- function(posterior_samples, 
                                        n_sams,
                                        obj,
                                        allocate_parameters,
                                        goveqs,
                                        goveqs_scaleup,
                                        bounds,
                                        prm, 
                                        ref, 
                                        thetai, 
                                        lhd_list,
                                        lhd_switch,
                                        pop_country){
  # Sample evenly spaced 250 parameters vectors from the posteriors after burnin and thinning
  total_samples <- nrow(posterior_samples)
  dx <- round(total_samples/n_sams)
  indices <- seq(1, total_samples, by = dx)
  indices <- indices[seq_len(min(n_sams, length(indices)))] # Ensure exactly 250 samples
  
  xs <- posterior_samples[indices, 1:nrow(bounds)]
  n_samples <- nrow(xs)
  dsol_list <- vector("list", length = n_samples)
  
  # Main loop over all selected posterior samples
  for (ii in seq_len(n_samples)) {
    # Print to see progress of intervention iterations
    if(ii %% 50 == 0) 
      print(paste("Intervention iteration: ", ii)) 
    
    theta <- as.numeric(xs[ii, ])
    
    # Allocate parameters based on input theta (parameter vector)
    alloc <- allocate_parameters(theta, thetai, prm$p, prm$r)
    r <- alloc$r
    p <- alloc$p
    i <- ref$i
    d <- ref$d
    
    # 1) Actual past trajectory (with RNTCP)
    init <- rep(0, i[["nx"]])
    names(init) <- c(
      "U", "Lf_ds", "Lf_rr", "Ls_ds", "Ls_rr", "I_ds", "I_rr",
      "Rlo_ds", "Rlo_rr", "Rhi_ds", "Rhi_rr", "R_ds", "R_rr",
      "Tx_ds_pu", "Tx_ds_pr_en", "Tx_ds_pr_nen", "Tx_rr_pu", "Tx_rr_pr_en", "Tx_rr_pr_nen",
      "Tx2_ds_pu", "Tx2_ds_pr_en", "Tx2_ds_pr_nen", "Tx2_rr_pu", "Tx2_rr_pr_en", "Tx2_rr_pr_nen",
      "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort", "cu_mort_ds", "cu_mort_rr", 
      "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    seed <- 1e-6
    init[i[["U"]]] <- 1 - seed; init[i[["I"]][["ds"]]] <- seed
    
    # equilibrium (no MDR, no RNTCP)
    # In absence of MDR and RNTCP
    p0 <- p; r0 <- r 
    p0[["pu"]] <- 0
    r0[["beta"]][2] <- 0
    r0[["MDR_acqu"]] <- 0
    r0[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters0 <- list(r = r0, p = p0)
    soln0 <- deSolve::ode(y = init, times = seq(0, 2000, by = 1000), func = goveqs, 
                          parms = parameters0,
                          method = "bdf")
    init1 <- tail(soln0[,-1],1)
    
    # MDR intro (no RNTCP)
    p1 <- p; r1 <- r 
    p1[["pu"]] <- 0
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters1 <- list(r = r1, p = p1)
    soln1 <- deSolve::ode(y = init1, times = seq(1970, 1997, by = 1), func = goveqs, 
                          parms = parameters1,
                          method = "bdf")
    init2 <- tail(soln1[,-1],1)
    
    # RNTCP scale-up
    p2 <- p; r2 <- r
    r2[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters2 <- list(r_start = r2, r_end = r2, p_start = p1, p_end = p2, times = c(1997, 2007))
    parameters2$goveqs <- goveqs
    soln2 <- deSolve::ode(y = init2, times = seq(1997, 2015, by = 1), func = goveqs_scaleup, 
                          parms = parameters2,
                          method = "bdf")
    init3 <- tail(soln2[,-1],1)
    
    # Treatment initiation scaleup
    r3_1 <- r
    r3_1[["Tx_init"]] <- r[["Tx_init_old"]]
    r3_2 <- r
    r3_2[["Tx_init"]] <- r[["Tx_init_new"]]
    p3 <- p
    parameters3 <- list(r_start = r3_1, r_end = r3_2, p_start = p3, p_end = p3, times = c(2015, 2025))
    parameters3$goveqs <- goveqs
    soln3 <- deSolve::ode(y = init3, times = seq(2015, 2025, by = 1), func = goveqs_scaleup, 
                          parms = parameters3,
                          method = "bdf")
    
    solution_df0 <- as.data.frame(soln0)
    colnames(solution_df0) <- c("time", d[1:i$nstates], 
                                "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df1 <- as.data.frame(soln1)
    colnames(solution_df1) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df1 <- solution_df1 |> dplyr::filter(time > 1970)
    solution_df2 <- as.data.frame(soln2)
    colnames(solution_df2) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df2 <- solution_df2 |> dplyr::filter(time > 1997)
    solution_df3 <- as.data.frame(soln3)
    colnames(solution_df3) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df3 <- solution_df3 |> dplyr::filter(time > 2015)
    
    sol_actual <- rbind(solution_df1, solution_df2, solution_df3)
    
    dsol_actual <- sol_actual %>%
      dplyr::select(time, cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr) %>%
      mutate(across(-time, ~ c(.x[1], diff(.x))),
             across(-time, ~ .x*1e5),
             scenario="Past actual") %>%
      rename(tb_mort = cu_mort_tb,
             inc_rr = cu_inc_rr,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    init_forward <- tail(soln3[,-1],1)  # forward intervention scenarios use this
    
    # 2) Counterfactual past (no p_MDR_rec, everything else same)
    seed <- 1e-6
    init[i[["U"]]] <- 1 - seed; init[i[["I"]][["ds"]]] <- seed
    
    # equilibrium (no MDR, no RNTCP)
    # In absence of MDR and RNTCP
    p0 <- p
    r0 <- r 
    p0[["pu"]] <- 0
    r0[["beta"]][2] <- 0
    r0[["MDR_acqu"]] <- 0
    r0[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters0 <- list(r = r0, p = p0)
    soln0 <- deSolve::ode(y = init, times = seq(0, 2000, by = 1000), func = goveqs, 
                          parms = parameters0,
                          method = "bdf")
    init1 <- tail(soln0[,-1],1)
    
    # MDR intro (no RNTCP)
    p1 <- p; r1 <- r 
    p1[["pu"]] <- 0
    p1[["MDR_rec"]][1] <- 0
    p1[["MDR_rec"]][2] <- 0
    p1[["MDR_rec_pr_nen"]] <- 0
    r1[["muTx"]][2] <- r[["mort_TB"]]
    r1[["muTx"]][3] <- r[["mort_TB"]]
    r1[["muTx"]][4] <- r[["mort_TB"]]
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters1 <- list(r = r1, p = p1)
    soln1 <- deSolve::ode(y = init1, times = seq(1970, 1997, by = 1), func = goveqs, 
                          parms = parameters1,
                          method = "bdf")
    init2 <- tail(soln1[,-1],1)
    
    # RNTCP scale-up
    p2 <- p; r2 <- r
    p2[["MDR_rec"]][1] <- 0
    p2[["MDR_rec"]][2] <- 0
    p2[["MDR_rec_pr_nen"]] <- 0
    r2[["muTx"]][2] <- r[["mort_TB"]]
    r2[["muTx"]][3] <- r[["mort_TB"]]
    r2[["muTx"]][4] <- r[["mort_TB"]]
    r2[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters2 <- list(r_start = r2, r_end = r2, p_start = p1, p_end = p2, times = c(1997, 2007))
    parameters2$goveqs <- goveqs
    soln2 <- deSolve::ode(y = init2, times = seq(1997, 2015, by = 1), func = goveqs_scaleup, 
                          parms = parameters2,
                          method = "bdf")
    init3 <- tail(soln2[,-1],1)
    
    # Treatment initiation scaleup
    r3_1 <- r
    r3_1[["Tx_init"]] <- r[["Tx_init_old"]]
    r3_2 <- r
    r3_2[["Tx_init"]] <- r[["Tx_init_new"]]
    p3 <- p
    p3[["MDR_rec"]][1] <- 0
    p3[["MDR_rec"]][2] <- 0
    p3[["MDR_rec_pr_nen"]] <- 0
    r3_1[["muTx"]][2] <- r[["mort_TB"]]
    r3_2[["muTx"]][2] <- r[["mort_TB"]]
    r3_1[["muTx"]][3] <- r[["mort_TB"]]
    r3_2[["muTx"]][3] <- r[["mort_TB"]]
    r3_1[["muTx"]][4] <- r[["mort_TB"]]
    r3_2[["muTx"]][4] <- r[["mort_TB"]]
    parameters3 <- list(r_start = r3_1, r_end = r3_2, p_start = p3, p_end = p3, times = c(2015, 2025))
    parameters3$goveqs <- goveqs
    soln3 <- deSolve::ode(y = init3, times = seq(2015, 2025, by = 1), func = goveqs_scaleup, 
                          parms = parameters3,
                          method = "bdf")
    
    solution_df0 <- as.data.frame(soln0)
    colnames(solution_df0) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df1 <- as.data.frame(soln1)
    colnames(solution_df1) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df1 <- solution_df1 |> dplyr::filter(time > 1970)
    solution_df2 <- as.data.frame(soln2)
    colnames(solution_df2) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df2 <- solution_df2 |> dplyr::filter(time > 1997)
    solution_df3 <- as.data.frame(soln3)
    colnames(solution_df3) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", 
                                "cu_mort_tb_ds", "cu_mort_tb_rr", "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    solution_df3 <- solution_df3 |> dplyr::filter(time > 2015)
    
    sol_cf <- rbind(solution_df1, solution_df2, solution_df3)
    
    dsol_cf <- sol_cf %>%
      as.data.frame() %>%
      dplyr::select(time, cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr) %>%
      mutate(across(-time, ~ c(.x[1], diff(.x))),
             across(-time, ~ .x*1e5),
             scenario="Past no p_MDR_rec")%>%
      rename(tb_mort = cu_mort_tb,
             inc_rr = cu_inc_rr,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    # 3) Forward projections using actual trajectory till 2025
    # add 0s for capturing treatment starters
    init <- c(init_forward)
    
    # scenario 0 - base model
    r0 <- r
    p0 <- p
    # p0[["MDR_rec_pr_nen"]] <- p0[["MDR_rec"]][2] # I think this was done to maintain that nen is not better than en? very similar
    r0[["Tx_init"]] <- r[["Tx_init_new"]]
    parameters0 <- list(r = r0, p = p0)
    soln0_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs, 
                            parms = parameters0,
                            method = "bdf")
    
    # scenario 1 - improve DST as final condition in scaleup
    p1 <- p0
    r1 <- r0
    r1[["Tx_init"]] <- r[["Tx_init_new"]]
    # p1$MDR_rec[1] <- max(p1$MDR_rec[1], 0.70) 
    p1$MDR_rec[1] <- max(p1$MDR_rec[1], 0.75) # estimating ~0.68 to impact is less
    parameters1 <- list(r_start = r0, r_end = r1, p_start = p0, p_end = p1, times = c(2026, 2028))
    parameters1$goveqs <- goveqs
    soln1_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs_scaleup, 
                            parms = parameters1,
                            method = "bdf")
    
    
    # scenario 2 - improve SL treatment as final condition in PUBLIC sector only
    # a) shift 50% of bad outcomes (LTFU and mortality) into success
    # b) halve LTFU when on Tx2
    # c) halve mortality when on Tx2
    p2 <- p1; r2 <- r1
    # r2$relapse[2] <- r1$relapse[2]*2 # try this to see impact on inc_rr?
    vec <- c(r2$Tx2, r2$ltfu[2], r2$muTx[2])
    propvec <- vec/sum(vec) # current proportions of Tx2 outcomes
    newvec <- c(propvec[1] + 0.5 * (1 - propvec[1]), 0.5 * propvec[2], 0.5 * propvec[3]) # new props of Tx2 outcomes
    r2$ltfu[2] <- r2$Tx2 * newvec[2] / newvec[1] # scale lower LTFU rate with new props
    r2$muTx[2] <- r2$Tx2 * newvec[3] / newvec[1] # scale lower muTx2 rate with new props
    p2$cure[2] <- (p2$cure[2] + 1)*0.5 # reducing failure by 50% i.e increasing success. r_Tx2 remains the same
    
    parameters2 <- list(r_start = r0, r_end = r2, p_start = p0, p_end = p2, times = c(2026, 2028))
    parameters2$goveqs <- goveqs
    soln2_f <- deSolve::ode(y = as.numeric(init), times = seq(2025, 2036, by = 1), func = goveqs_scaleup, 
                            parms = parameters2,
                            method = "bdf")
    
    
    # 4) Combine outputs
    soln0_f <- as.data.frame(soln0_f) %>% mutate(scenario = "Baseline")
    soln1_f <- as.data.frame(soln1_f) %>% mutate(scenario = "Improved DST (no PPM expansion)")
    soln2_f <- as.data.frame(soln2_f) %>% mutate(scenario = "Improved DST + SL outcomes (no PPM expansion)")
    
    soln_f <- bind_rows(soln0_f, soln1_f, soln2_f)
    colnames(soln_f) <- c("time", d[1:i$nstates], "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", 
                          "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr",
                          "cu_noti_pr_en",
                          "cu_noti_rr_pu",
                          "cu_noti_rr_pr_en",
                          "scenario")
    
    dsol_f <- soln_f %>%
      dplyr::select("time", "scenario", "cu_inc_rr", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr") %>%
      group_by(scenario) %>%
      mutate(across(c(cu_inc_rr, cu_mort_tb, cu_mort_tb_ds, cu_mort_tb_rr),
                    ~ c(.x[1], diff(.x))),
             across(-time, ~ .x * 1e5)) %>% # to get per 100k per year
      rename(inc_rr = cu_inc_rr,
             tb_mort = cu_mort_tb,
             tb_mort_ds = cu_mort_tb_ds,
             tb_mort_rr = cu_mort_tb_rr)
    
    dsol <- bind_rows(dsol_actual, dsol_cf, dsol_f)
    dsol$iter <- ii
    dsol_list[[ii]] <- dsol
  }
  
  combined_dsol_all <- bind_rows(dsol_list)
  
  # Summarise inc_rr, tb_mort, tb_mort_ds, tb_mort_rr in the past and future scenarios
  # pivot long
  long_dsol <- combined_dsol_all %>%
    pivot_longer(cols=c(inc_rr, tb_mort, tb_mort_ds, tb_mort_rr),
                 names_to="indicator", values_to="value")
  
  # summarise the raw values
  summary_dsol <- long_dsol %>%
    group_by(time, scenario, indicator) %>%
    summarise(
      median = median(value),
      lower  = quantile(value, 0.025),
      upper  = quantile(value, 0.975),
      .groups = "drop"
    )

  
  return(list(draws = combined_dsol_all,
              summary_dsol = summary_dsol))
  
}

calculate_india_sensitivity <- function(primary_draws, sensitivity_draws) {
  
  # Population data (same as make_agg_global)
  population_df <- dplyr::bind_rows(
    readxl::read_excel(
      "data/WPP2024_GEN_F01_DEMOGRAPHIC_INDICATORS_COMPACT.xlsx",
      sheet = "Estimates", skip = 16
    ),
    readxl::read_excel(
      "data/WPP2024_GEN_F01_DEMOGRAPHIC_INDICATORS_COMPACT.xlsx",
      sheet = "Medium variant", skip = 16
    )
  ) %>%
    dplyr::filter(
      `Type` == "Country/Area",
      `Region, subregion, country or area *` == "India"
    ) %>%
    dplyr::rename(
      time = Year,
      population_k = `Total Population, as of 1 July (thousands)`
    ) %>%
    dplyr::distinct(time, .keep_all = TRUE) %>%
    dplyr::transmute(
      time = as.numeric(time),
      population = as.numeric(population_k) * 1000
    )
  
  # Select combined intervention from each model
  primary <- primary_draws %>%
    dplyr::filter(
      scenario == "Improved DST + SL outcomes",
      time > 2025
    ) %>%
    dplyr::mutate(analysis = "With PPM")
  
  sensitivity <- sensitivity_draws %>%
    dplyr::filter(
      scenario == "Improved DST + SL outcomes (no PPM expansion)",
      time > 2025
    ) %>%
    dplyr::mutate(analysis = "Without PPM")
  
  combined <- dplyr::bind_rows(primary, sensitivity) %>%
    dplyr::left_join(population_df, by = "time") %>%
    tidyr::pivot_longer(
      cols = c(inc_rr, tb_mort_rr),
      names_to = "indicator",
      values_to = "rate"
    ) %>%
    dplyr::mutate(count = rate * population / 1e5)
  
  if (anyNA(combined$count)) {
    stop("Missing population or incidence/mortality values")
  }
  
  # Cumulative outcomes for each posterior draw
  cumulative <- combined %>%
    dplyr::group_by(iter, analysis, indicator) %>%
    dplyr::summarise(
      cumulative = sum(count),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = analysis,
      values_from = cumulative
    ) %>%
    dplyr::mutate(
      additional_averted = `Without PPM` - `With PPM`,
      additional_reduction_pct =
        100 * additional_averted / `Without PPM`
    )
  
  # Posterior median and 95% credible intervals
  cumulative %>%
    dplyr::group_by(indicator) %>%
    dplyr::summarise(
      additional_averted_median = median(additional_averted),
      additional_averted_lower = quantile(additional_averted, 0.025),
      additional_averted_upper = quantile(additional_averted, 0.975),
      reduction_pct_median = median(additional_reduction_pct),
      reduction_pct_lower = quantile(additional_reduction_pct, 0.025),
      reduction_pct_upper = quantile(additional_reduction_pct, 0.975),
      .groups = "drop"
    )
}