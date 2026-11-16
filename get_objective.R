get_objective_m1_1 <- function(theta, 
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
                               output_names) {
  failed_result <- function() {
    c(list(lhd = -1e10),
      setNames(as.list(rep(NA_real_, length(output_names))), output_names),
      list(sfin = NA))
  }
  # Check if the current values of parameters are within bounds (brute force priors)
  mat <- rbind(bounds[, 2] - theta, theta - bounds[, 1])
  
  # Check if any value is out of bounds
  safe_min <- suppressWarnings(min(mat, na.rm = TRUE))
  if (is.na(safe_min) || safe_min < 0) {
    return(failed_result())
    print("Bounds violated")
  }
  else {
    # Allocate parameters based on input theta (parameter vector)
    param_list <- allocate_parameters(theta, thetai, prm$p, prm$r)
    r <- param_list$r
    p <- param_list$p
    i <- ref$i
    s <- ref$s
    d <- ref$d
    
    # Solve the models
    # 1. Equilibrium model
    init <- rep(0, i[["nx"]])
    names(init) <- c(
      "U", "Lf_ds", "Lf_rr", "Ls_ds", "Ls_rr", "I_ds", "I_rr",
      "Rlo_ds", "Rlo_rr", "Rhi_ds", "Rhi_rr", "R_ds", "R_rr",
      "Tx_ds_pu", "Tx_ds_pr_en", "Tx_ds_pr_nen", "Tx_rr_pu", "Tx_rr_pr_en", "Tx_rr_pr_nen",
      "Tx2_ds_pu", "Tx2_ds_pr_en", "Tx2_ds_pr_nen", "Tx2_rr_pu", "Tx2_rr_pr_en", "Tx2_rr_pr_nen",
      "cu_inc_ds", "cu_inc_rr", "cu_noti_pu", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr",
      "cu_noti_pr_en", "cu_noti_rr_pu", "cu_noti_rr_pr_en")
    seed <- 1e-6
    init[i[["U"]]] <- 1 - seed; init[i[["I"]][["ds"]]] <- seed
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
    soln0[, -1] <- pmax(soln0[, -1], 0)
    init1 <- tail(soln0[, -1], 1)
    if (is.null(init1)) return(failed_result())
    
    # 2. Introduce MDR-TB
    # In absence of RNTCP
    p1 <- p; r1 <- r 
    p1[["pu"]] <- 0
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters1 <- list(r = r1, p = p1)
    soln1 <- deSolve::ode(y = init1, times = seq(1970, 1997, by = 1), func = goveqs, 
                          parms = parameters1,
                          method = "bdf")
    soln1[, -1] <- pmax(soln1[, -1], 0)
    init2 <- tail(soln1[, -1], 1)
    if (is.null(init2)) return(failed_result())
    
    # 3) RNTCP scale-up
    p2 <- p; r2 <- r
    r2[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters2 <- list(r_start = r2, r_end = r2, p_start = p1, p_end = p2, times = c(1997, 2007))
    parameters2$goveqs <- goveqs
    soln2 <- deSolve::ode(y = init2, times = seq(1997, 2015, by = 1), func = goveqs_scaleup, 
                          parms = parameters2,
                          method = "bdf")
    soln2[, -1] <- pmax(soln2[, -1], 0)
    init3 <- tail(soln2[, -1], 1)
    if (is.null(init3)) return(failed_result())
    
    # 4) Treatment initiation scale-up
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
    soln3[, -1] <- pmax(soln3[, -1], 0)
    init3 <- tail(soln3[, -1], 1)
    if (is.null(init3)) return(failed_result())
    
    
    # Convert historical solutions
    # Not keeping the equilibrium solution
    solution_names <- c("time",
                        d[1:i$nstates],
                        "cu_inc_ds", 
                        "cu_inc_rr", 
                        "cu_noti_pu", 
                        "cu_mort_tb", 
                        "cu_mort_tb_ds", 
                        "cu_mort_tb_rr",
                        "cu_noti_pr_en", 
                        "cu_noti_rr_pu", 
                        "cu_noti_rr_pr_en")
    
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
    
    # NOTE sfin does not have time column, but sol and dsol have- so indexing will increase by 1
    sol <- dplyr::bind_rows(solution_df1,
                            solution_df2,
                            solution_df3)
    sfin <- tail(sol[, -1], 1) # Extract final state 
    dsol <- sol |> dplyr::mutate(across(-time, ~ c(.x[1], diff(.x)))) # Keep the time column as is
    
    # Prevalence
    sim_prev <- sum(sfin[ , s[["prevalent"]]]) * 1e5
    
    # Get objectives
    # Incidence for 2024
    last_row_dsol <- tail(dsol, 1)
    tmp <- last_row_dsol[, i[["aux"]][["inc"]] + 1] * 1e5
    sim_incd_new_all <- sum(tmp)
    sim_incd_rr <- as.numeric(tmp[2])
    
    # Incidence for older time point
    row_old <- dsol |> dplyr::filter(time == old_inc_year + 1)
    tmp_old <- row_old[, i[["aux"]][["inc"]] + 1] * 1e5
    sim_incd_old_all <- sum(tmp_old)
    sim_incd_old_rr  <- as.numeric(tmp_old[2])
    
    # TB mortality
    sim_tb_mort <- last_row_dsol[, i[["aux"]][["mort_tb"]] + 1] * 1e5
    
    # TB mortality - DS
    sim_tb_mort_ds <- last_row_dsol[, i[["aux"]][["mort_tb_ds"]] + 1] * 1e5
    
    # TB mortality - RR
    sim_tb_mort_rr <- last_row_dsol[, i[["aux"]][["mort_tb_rr"]] + 1] * 1e5
    
    # Notifications from public sector in 2023
    row_2023 <- dsol |> dplyr::filter(time == 2024)
    sim_notifs_pu <- row_2023[, i[["aux"]][["noti_pu"]] + 1] * 1e5
    
    # Notifications from engaged private sector in 2023
    sim_notifs_pr_en <- row_2023[, i[["aux"]][["noti_pr_en"]] + 1] * 1e5
    
    # RR notifications from public sector in 2023
    sim_notifs_rr_pu <- row_2023[, i[["aux"]][["noti_rr_pu"]] + 1] * 1e5
    
    # RR notifications from engaged private sector in 2023
    sim_notifs_rr_pr_en <- row_2023[, i[["aux"]][["noti_rr_pr_en"]] + 1] * 1e5
    
    # Combine all model simulated values in the correct order
    model_outputs <- list(sim_incd_old_all,
                          sim_incd_new_all,
                          sim_incd_rr,
                          sim_tb_mort,
                          sim_notifs_pu, 
                          sim_notifs_pr_en, 
                          sim_notifs_rr_pu, 
                          sim_notifs_rr_pr_en)
    
    all_outputs <- unlist(model_outputs)
    if (any(!is.finite(all_outputs)) || any(all_outputs < 0)) {
      print("Invalid ODE outputs detected. Penalizing.")
      return(failed_result())
    }
    
    # Default switch- include everything (1s) if not provided
    if (is.null(lhd_switch)) lhd_switch <- rep(1, length(model_outputs))
    if (is.null(lhd_weight)) lhd_weight <- rep(1, length(model_outputs))
    
    # Build safe wrappers for each provided likelihood function
    safe_funcs <- lapply(lhd_list, function(f) {
      if (is.null(f)) {
        # if no function provided, return a function giving 0 contribution
        function(x) 0
      } else if (!is.function(f)) {
        # if something else provided accidentally, treat as zero and warn
        warning("One entry in lhd_list is not a function; treating it as 0.")
        function(x) 0
      } else {
        function(x) {
          if (!is.finite(x) || x < 0) return(-1e10)
          return(f(x))
        }
      }
    })
    
    # Compute contributions and total log-likelihood
    contributions <- numeric(length(model_outputs))
    
    for (i in seq_along(model_outputs)) {
      sw <- as.numeric(lhd_switch[i])
      if (sw == 0) {
        contributions[i] <- 0
      } else {
        val <- safe_funcs[[i]](model_outputs[[i]])
        if (!is.finite(val) || length(val) != 1) val <- -1e10
        val <- val * as.numeric(lhd_weight[i])
        contributions[i] <- val
      }
    }
    
    # Make sure notifs are bounded by bac confirmed notifs and all notifs for all TB in public and private engaged sectors
    if (!is.finite(sim_notifs_pu)    || sim_notifs_pu    < data_notif_bounds$notif_l_pu    || sim_notifs_pu    > data_notif_bounds$notif_u_pu    ||
        !is.finite(sim_notifs_pr_en) || sim_notifs_pr_en < data_notif_bounds$notif_l_pr_en || sim_notifs_pr_en > data_notif_bounds$notif_u_pr_en) {
      lhd <- -1e10
    } else {
      lhd <- sum(contributions)
    }
    
    return(list(lhd = lhd, incd_old = sim_incd_old_all, incd_new = sim_incd_new_all,
                incd_rr = sim_incd_rr, tb_mort = sim_tb_mort,
                tb_mort_ds = sim_tb_mort_ds, tb_mort_rr = sim_tb_mort_rr,
                notifs_pu = sim_notifs_pu, prev = sim_prev, notifs_pr_en = sim_notifs_pr_en, 
                notifs_rr_pu = sim_notifs_rr_pu, notifs_rr_pr_en = sim_notifs_rr_pr_en, 
                sfin = sfin))
    
  }
}

get_objective_m2_1 <- function(theta, 
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
                               output_names) {
  
  failed_result <- function() {
    c(list(lhd = -1e10),
      setNames(as.list(rep(NA_real_, length(output_names))), output_names),
      list(sfin = NA))
  }
  
  # Check if the current values of parameters are within bounds (brute force priors)
  mat <- rbind(bounds[, 2] - theta, theta - bounds[, 1])
  
  # Check if any value is out of bounds
  safe_min <- suppressWarnings(min(mat, na.rm = TRUE))
  if (is.na(safe_min) || safe_min < 0) {
    return(failed_result())
    print("Bounds violated")
  }
  else {
    # Allocate parameters based on input theta (parameter vector)
    param_list <- allocate_parameters(theta, thetai, prm$p, prm$r)
    r <- param_list$r
    p <- param_list$p
    i <- ref$i
    s <- ref$s
    d <- ref$d
    
    # Solve the models
    # 1. Equilibrium model
    init <- rep(0, i[["nx"]])
    names(init) <- c(
      "U", "Lf_ds", "Lf_rr", "Ls_ds", "Ls_rr", "I_ds", "I_rr",
      "Rlo_ds", "Rlo_rr", "Rhi_ds", "Rhi_rr", "R_ds", "R_rr",
      "Tx_ds_pu", "Tx_ds_pr_nen", "Tx_rr_pu", "Tx_rr_pr_nen",
      "Tx2_ds_pu", "Tx2_ds_pr_nen", "Tx2_rr_pu", "Tx2_rr_pr_nen",
      "cu_inc_ds", "cu_inc_rr", "cu_noti", "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr",
      "cu_noti_rr")
    seed <- 1e-6
    init[i[["U"]]] <- 1 - seed; init[i[["I"]][["ds"]]] <- seed
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
    soln0[, -1] <- pmax(soln0[, -1], 0)
    init1 <- tail(soln0[, -1], 1)
    if (is.null(init1)) return(failed_result())
    
    # 2. Introduce MDR-TB
    # In absence of RNTCP
    p1 <- p; r1 <- r 
    p1[["pu"]] <- 0
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters1 <- list(r = r1, p = p1)
    soln1 <- deSolve::ode(y = init1, times = seq(1970, 1997, by = 1), func = goveqs, 
                          parms = parameters1,
                          method = "bdf")
    soln1[, -1] <- pmax(soln1[, -1], 0)
    init2 <- tail(soln1[, -1], 1)
    if (is.null(init2)) return(failed_result())
    
    # 3) RNTCP scale-up
    p2 <- p; r2 <- r
    r2[["Tx_init"]] <- r[["Tx_init_old"]]
    parameters2 <- list(r_start = r2, r_end = r2, p_start = p1, p_end = p2, times = c(1997, 2007))
    parameters2$goveqs <- goveqs
    soln2 <- deSolve::ode(y = init2, times = seq(1997, 2015, by = 1), func = goveqs_scaleup, 
                          parms = parameters2,
                          method = "bdf")
    soln2[, -1] <- pmax(soln2[, -1], 0)
    init3 <- tail(soln2[, -1], 1)
    if (is.null(init3)) return(failed_result())
    
    # 4) Treatment initiation scale-up
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
    soln3[, -1] <- pmax(soln3[, -1], 0)
    init3 <- tail(soln3[, -1], 1)
    if (is.null(init3)) return(failed_result())
    
    
    # Convert historical solutions
    # Not keeping the equilibrium solution
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
    
    # NOTE sfin does not have time column, but sol and dsol have- so indexing will increase by 1
    sol <- dplyr::bind_rows(solution_df1,
                            solution_df2,
                            solution_df3)
    sfin <- tail(sol[, -1], 1) # Extract final state 
    dsol <- sol |> dplyr::mutate(across(-time, ~ c(.x[1], diff(.x)))) # Keep the time column as is
    
    # Prevalence
    sim_prev <- sum(sfin[ , s[["prevalent"]]]) * 1e5
    
    # Get objectives
    # Incidence for 2024
    last_row_dsol <- tail(dsol, 1)
    tmp <- last_row_dsol[, i[["aux"]][["inc"]] + 1] * 1e5
    sim_incd_new_all <- sum(tmp)
    sim_incd_rr <- as.numeric(tmp[2])
    
    # Incidence for older time point
    row_old <- dsol |> dplyr::filter(time == old_inc_year + 1)
    tmp_old <- row_old[, i[["aux"]][["inc"]] + 1] * 1e5
    sim_incd_old_all <- sum(tmp_old)
    sim_incd_old_rr  <- as.numeric(tmp_old[2])
    
    # TB mortality
    sim_tb_mort <- last_row_dsol[, i[["aux"]][["mort_tb"]] + 1] * 1e5
    
    # TB mortality - DS
    sim_tb_mort_ds <- last_row_dsol[, i[["aux"]][["mort_tb_ds"]] + 1] * 1e5
    
    # TB mortality - RR
    sim_tb_mort_rr <- last_row_dsol[, i[["aux"]][["mort_tb_rr"]] + 1] * 1e5
    
    # Notifications from public sector in 2024
    row_2024 <- dsol |> dplyr::filter(time == 2025)
    sim_notifs <- row_2024[, i[["aux"]][["noti"]] + 1] * 1e5
    
    # RR notifications from public sector in 2024
    sim_notifs_rr <- row_2024[, i[["aux"]][["noti_rr"]] + 1] * 1e5
    
    # Combine all model simulated values in the correct order
    model_outputs <- list(sim_incd_old_all,
                          sim_incd_new_all,
                          sim_incd_rr,
                          sim_tb_mort,
                          sim_notifs, 
                          sim_notifs_rr)
    
    all_outputs <- unlist(model_outputs)
    if (any(!is.finite(all_outputs)) || any(all_outputs < 0)) {
      print("Invalid ODE outputs detected. Penalizing.")
      return(failed_result())
    }
    
    # Default switch- include everything (1s) if not provided
    if (is.null(lhd_switch)) lhd_switch <- rep(1, length(model_outputs))
    if (is.null(lhd_weight)) lhd_weight <- rep(1, length(model_outputs))
    
    # Build safe wrappers for each provided likelihood function
    safe_funcs <- lapply(lhd_list, function(f) {
      if (is.null(f)) {
        # if no function provided, return a function giving 0 contribution
        function(x) 0
      } else if (!is.function(f)) {
        # if something else provided accidentally, treat as zero and warn
        warning("One entry in lhd_list is not a function; treating it as 0.")
        function(x) 0
      } else {
        function(x) {
          if (!is.finite(x) || x < 0) return(-1e10)
          return(f(x))
        }
      }
    })
    
    # Compute contributions and total log-likelihood
    contributions <- numeric(length(model_outputs))
    
    for (i in seq_along(model_outputs)) {
      sw <- as.numeric(lhd_switch[i])
      if (sw == 0) {
        contributions[i] <- 0
      } else {
        val <- safe_funcs[[i]](model_outputs[[i]])
        if (!is.finite(val) || length(val) != 1) val <- -1e10
        val <- val * as.numeric(lhd_weight[i])
        contributions[i] <- val
      }
    }
    
    # Make sure notifs are bounded by bac confirmed notifs and all notifs for all TB
    if (!is.finite(sim_notifs) || sim_notifs < data_notif_bounds$notif_l || sim_notifs > data_notif_bounds$notif_u) {
      lhd <- -1e10
    } else {
      lhd <- sum(contributions)
    }
    
    return(list(lhd = lhd, incd_old = sim_incd_old_all, incd_new = sim_incd_new_all,
                incd_rr = sim_incd_rr, tb_mort = sim_tb_mort,
                tb_mort_ds = sim_tb_mort_ds, tb_mort_rr = sim_tb_mort_rr,
                notifs = sim_notifs, prev = sim_prev, 
                notifs_rr = sim_notifs_rr,
                sfin = sfin))
    
  }
}

get_objective_m3_1 <- function(theta, 
                               allocate_parameters,
                               goveqs,
                               goveqs_scaleup,
                               bounds, prm, ref, 
                               thetai,
                               lhd_list,
                               lhd_switch,
                               lhd_weight,
                               data_notif_bounds,
                               old_inc_year,
                               output_names) {
  
  failed_result <- function() {
    c(list(lhd = -1e10),
      setNames(as.list(rep(NA_real_, length(output_names))), output_names),
      list(sfin = NA))
  }
  
  # Check if the current values of parameters are within bounds (brute force priors)
  mat <- rbind(bounds[, 2] - theta, theta - bounds[, 1])
  
  # Check if any value is out of bounds
  safe_min <- suppressWarnings(min(mat, na.rm = TRUE))
  if (is.na(safe_min) || safe_min < 0) {
    return(failed_result())
    print("Bounds violated")
  }
  else {
    # Allocate parameters based on input theta (parameter vector)
    param_list <- allocate_parameters(theta, thetai, prm$p, prm$r)
    r <- param_list$r
    p <- param_list$p
    i <- ref$i
    s <- ref$s
    d <- ref$d
    
    # Solve the models
    # 1. Equilibrium model
    init <- rep(0, i[["nx"]])
    names(init) <- c(
      "U", "Lf_ds", "Lf_rr", "Ls_ds", "Ls_rr", "I_ds", "I_rr",
      "Rlo_ds", "Rlo_rr", "Rhi_ds", "Rhi_rr", "R_ds", "R_rr",
      "Tx_ds", "Tx_rr", "Tx2_ds", "Tx2_rr",
      "cu_inc_ds", "cu_inc_rr",
      "cu_noti",
      "cu_mort_tb", "cu_mort_tb_ds", "cu_mort_tb_rr",
      "cu_noti_rr")
    seed <- 1e-6
    init[i[["U"]]] <- 1 - seed; init[i[["I"]][["ds"]]] <- seed
    
    # In absence of MDR
    p0 <- p
    r0 <- r
    
    r0[["beta"]][2] <- 0
    r0[["MDR_acqu"]] <- 0
    r0[["Tx_init"]] <- r[["Tx_init_old"]]
    
    parameters0 <- list(r = r0, p = p0)
    
    soln0 <- deSolve::ode(y = init, 
                          times = seq(0, 2000, by = 1000), 
                          func = goveqs, 
                          parms = parameters0,
                          method = "bdf")
    
    soln0[, -1] <- pmax(soln0[, -1], 0)
    init1 <- tail(soln0[, -1], 1)
    if (is.null(init1)) return(failed_result())
    
    # MDR introduced, constant MDR recognition
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
    soln1[, -1] <- pmax(soln1[, -1], 0)
    
    init2 <- tail(soln1[, -1], 1)
    if (is.null(init2)) return(failed_result())
    
    # constant MDR recognition
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
    soln2[, -1] <- pmax(soln2[, -1], 0)
    
    init3 <- tail(soln2[, -1], 1)
    if (is.null(init3)) return(failed_result())
    
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
    soln3[, -1] <- pmax(soln3[, -1], 0)
    
    init3 <- tail(soln3[, -1], 1)
    if (is.null(init3)) return(failed_result())
    
    # Convert historical solutions
    # Not keeping the equilibrium solution
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
    
    # NOTE sfin does not have time column, but sol and dsol have- so indexing will increase by 1
    sol <- dplyr::bind_rows(solution_df1,
                            solution_df2,
                            solution_df3)
    sfin <- tail(sol[, -1], 1) # Extract final state 
    dsol <- sol |> dplyr::mutate(across(-time, ~ c(.x[1], diff(.x)))) # Keep the time column as is
    
    # Get prevalence (per 100k)
    sim_prev <- sum(sfin[ , s[["prevalent"]]]) * 1e5
    
    # Incidence for 2024
    last_row_dsol <- tail(dsol, 1)
    tmp <- last_row_dsol[, i[["aux"]][["inc"]] + 1] * 1e5
    sim_incd_new_all <- sum(tmp)
    sim_incd_rr <- as.numeric(tmp[2])
    
    # Incidence for older time point
    row_old <- dsol |> dplyr::filter(time == old_inc_year + 1)
    tmp_old <- row_old[, i[["aux"]][["inc"]] + 1] * 1e5
    sim_incd_old_all <- sum(tmp_old)
    sim_incd_old_rr  <- as.numeric(tmp_old[2])
    
    # TB mortality
    sim_tb_mort <- last_row_dsol[, i[["aux"]][["mort_tb"]] + 1] * 1e5
    
    # TB mortality - DS
    sim_tb_mort_ds <- last_row_dsol[, i[["aux"]][["mort_tb_ds"]] + 1] * 1e5
    
    # TB mortality - RR
    sim_tb_mort_rr <- last_row_dsol[, i[["aux"]][["mort_tb_rr"]] + 1] * 1e5
    
    # Notifications from 2024
    row_2024 <- dsol |> dplyr::filter(time == 2025) # same as last_row_dsol
    sim_notifs <- row_2024[, i[["aux"]][["noti"]] + 1] * 1e5
    
    # RR notifications in 2024
    sim_notifs_rr <- row_2024[, i[["aux"]][["noti_rr"]] + 1] * 1e5
    
    # Combine all model simulated values in the correct order
    model_outputs <- list(sim_incd_old_all,
                          sim_incd_new_all,
                          sim_incd_rr,
                          sim_tb_mort,
                          sim_notifs,
                          sim_notifs_rr)
    
    all_outputs <- unlist(model_outputs)
    if (any(!is.finite(all_outputs)) || any(all_outputs < 0)) {
      print("Invalid ODE outputs detected. Penalizing.")
      return(failed_result())
    }
    
    # Default switch- include everything (1s) if not provided
    if (is.null(lhd_switch)) lhd_switch <- rep(1, length(model_outputs))
    if (is.null(lhd_weight)) lhd_weight <- rep(1, length(model_outputs))
    
    # Build safe wrappers for each provided likelihood function
    safe_funcs <- lapply(lhd_list, function(f) {
      if (is.null(f)) {
        # if no function provided, return a function giving 0 contribution
        function(x) 0
      } else if (!is.function(f)) {
        # if something else provided accidentally, treat as zero and warn
        warning("One entry in lhd_list is not a function; treating it as 0.")
        function(x) 0
      } else {
        function(x) {
          if (!is.finite(x) || x < 0) return(-1e10)
          return(f(x))
        }
      }
    })
    
    # Compute contributions and total log-likelihood
    contributions <- numeric(length(model_outputs))
    
    for (i in seq_along(model_outputs)) {
      sw <- as.numeric(lhd_switch[i])
      if (sw == 0) {
        contributions[i] <- 0
      } else {
        val <- safe_funcs[[i]](model_outputs[[i]])
        if (!is.finite(val) || length(val) != 1) val <- -1e10
        val <- val * as.numeric(lhd_weight[i])
        contributions[i] <- val
      }
    }
    
    # Make sure notifs are bounded by bac confirmed notifs and all notifs for all TB
    if (!is.finite(sim_notifs) || sim_notifs < data_notif_bounds$notif_l || sim_notifs > data_notif_bounds$notif_u) {
      lhd <- -1e10
    } else {
      lhd <- sum(contributions)
    }
    
    return(list(lhd = lhd, incd_old = sim_incd_old_all, incd_new = sim_incd_new_all,
                incd_rr = sim_incd_rr, tb_mort = sim_tb_mort,
                tb_mort_ds = sim_tb_mort_ds, tb_mort_rr = sim_tb_mort_rr,
                notifs = sim_notifs, prev = sim_prev, 
                notifs_rr = sim_notifs_rr,
                sfin = sfin))
    
  }
}


get_objective_m4_1 <- function(theta, 
                                 allocate_parameters,
                                 goveqs,
                                 goveqs_scaleup,
                                 bounds, prm, ref, 
                                 thetai,
                                 lhd_list,
                                 lhd_switch,
                                 lhd_weight,
                                 data_notif_bounds,
                                 old_inc_year,
                               output_names) {
  
  failed_result <- function() {
    c(list(lhd = -1e10),
      setNames(as.list(rep(NA_real_, length(output_names))), output_names),
      list(sfin = NA))
  }
  
  # Check if the current values of parameters are within bounds (brute force priors)
  mat <- rbind(bounds[, 2] - theta, theta - bounds[, 1])
  
  # Check if any value is out of bounds
  safe_min <- suppressWarnings(min(mat, na.rm = TRUE))
  if (is.na(safe_min) || safe_min < 0) {
    return(failed_result())
    print("Bounds violated")
  }
  else {
    # Allocate parameters based on input theta (parameter vector)
    param_list <- allocate_parameters(theta, thetai, prm$p, prm$r)
    r <- param_list$r
    p <- param_list$p
    i <- ref$i
    s <- ref$s
    d <- ref$d
    
    # Solve the models
    # 1. Equilibrium model
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
    seed <- 1e-3
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
    init1 <- tail(soln0[, -1], 1)
    if (is.null(init1)) return(failed_result())
    
    # 2. Introduce MDR-TB and no ART
    p1 <- p; r1 <- r
    r1[["Tx_init"]] <- r[["Tx_init_old"]]
    r1[["ART_init"]] <- 0
    parameters1 <- list(r = r1, p = p1, i = i)
    soln1 <- deSolve::ode(y = init1, times = seq(1970, 2000, by = 1), func = goveqs, 
                          parms = parameters1,
                          method = "bdf")
    soln1[, -1] <- pmax(soln1[, -1], 0)
    init2 <- tail(soln1[, -1], 1)
    if (is.null(init2)) return(failed_result())
    
    # 3. Treatment initiation scale-up and ART scaleup
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
    init3 <- tail(soln2[, -1], 1)
    if (is.null(init3)) return(failed_result())
    
    # Convert historical solutions
    # Not keeping the equilibrium solution
    solution_names <- c("time",
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
                        "cu_mort_tb_h1hart")
    
    solution_df1 <- as.data.frame(soln1)
    colnames(solution_df1) <- solution_names
    solution_df1 <- solution_df1 %>%
      dplyr::filter(time > 1970)
    
    solution_df2 <- as.data.frame(soln2)
    colnames(solution_df2) <- solution_names
    solution_df2 <- solution_df2 %>%
      dplyr::filter(time > 2000)
    
    
    # NOTE sfin does not have time column, but sol and dsol have- so indexing will increase by 1
    sol <- dplyr::bind_rows(solution_df1,
                            solution_df2)
    sfin <- tail(sol[, -1], 1) # Extract final state 
    dsol <- sol |> dplyr::mutate(across(-time, ~ c(.x[1], diff(.x)))) # Keep the time column as is
    
    # Get prevalence (per 100k)
    sim_prev <- sum(sfin[ , s[["prevalent"]]]) * 1e5
    
    # Incidence for 2024
    last_row_dsol <- tail(dsol, 1)
    tmp <- last_row_dsol[, i[["aux"]][["inc"]] + 1] * 1e5
    sim_incd_new_all <- sum(tmp)
    sim_incd_rr <- as.numeric(tmp[2])
    
    # Incidence for 2010
    row_old <- dsol |> dplyr::filter(time == old_inc_year + 1)
    tmp_old <- row_old[, i[["aux"]][["inc"]] + 1] * 1e5
    sim_incd_old_all <- sum(tmp_old)
    sim_incd_old_rr  <- as.numeric(tmp_old[2])
    
    # TB mortality
    sim_tb_mort <- last_row_dsol[, i[["aux"]][["mort_tb"]] + 1] * 1e5
    
    # TB mortality - DS
    sim_tb_mort_ds <- last_row_dsol[, i[["aux"]][["mort_tb_ds"]] + 1] * 1e5
    
    # TB mortality - RR
    sim_tb_mort_rr <- last_row_dsol[, i[["aux"]][["mort_tb_rr"]] + 1] * 1e5
    
    # Notifications from 2024
    row_2024 <- dsol |> dplyr::filter(time == 2025) # same as last_row_dsol
    sim_notifs <- row_2024[, i[["aux"]][["noti"]] + 1] * 1e5
    
    # RR notifications from public sector in 2024
    sim_notifs_rr <- row_2024[, i[["aux"]][["noti_rr"]] + 1] * 1e5
    
    # HIV prevalence and ART coverage from final state
    sim_HIV_prev <- sum(sfin[c(s[["h1"]], s[["hart"]])]) / sum(sfin[1:i[["nstates"]]])
    sim_ART_covg <- sum(sfin[s[["hart"]]]) / sum(sfin[c(s[["h1"]], s[["hart"]])])
    
    # These three should come from HIV-stratified cumulative outputs
    # in goveqs_sa_1() if you want them exactly.
    sim_incd_h1hart    <- last_row_dsol[, i[["aux"]][["inc_h1hart"]] + 1] * 1e5
    sim_tb_mort_h0     <- last_row_dsol[, i[["aux"]][["mort_tb_h0"]] + 1] * 1e5
    sim_tb_mort_h1hart <- last_row_dsol[, i[["aux"]][["mort_tb_h1hart"]] + 1] * 1e5
    
    # Combine all model simulated values in the correct order
    model_outputs <- list(sim_incd_old_all,
                          sim_incd_new_all,
                          sim_incd_rr,
                          sim_tb_mort,
                          sim_notifs,
                          sim_notifs_rr,
                          sim_incd_h1hart,
                          sim_tb_mort_h0,
                          sim_tb_mort_h1hart,
                          sim_HIV_prev,
                          sim_ART_covg)
    
    all_outputs <- unlist(model_outputs)
    if (any(!is.finite(all_outputs)) || any(all_outputs < 0)) {
      print("Invalid ODE outputs detected. Penalizing.")
      return(failed_result())
    }
    
    # Default switch- include everything (1s) if not provided
    if (is.null(lhd_switch)) lhd_switch <- rep(1, length(model_outputs))
    if (is.null(lhd_weight)) lhd_weight <- rep(1, length(model_outputs))
    
    # Build safe wrappers for each provided likelihood function
    safe_funcs <- lapply(lhd_list, function(f) {
      if (is.null(f)) {
        # if no function provided, return a function giving 0 contribution
        function(x) 0
      } else if (!is.function(f)) {
        # if something else provided accidentally, treat as zero and warn
        warning("One entry in lhd_list is not a function; treating it as 0.")
        function(x) 0
      } else {
        function(x) {
          if (!is.finite(x) || x < 0) return(-1e10)
          return(f(x))
        }
      }
    })
    
    # Compute contributions and total log-likelihood
    contributions <- numeric(length(model_outputs))
    
    for (i in seq_along(model_outputs)) {
      sw <- as.numeric(lhd_switch[i])
      if (sw == 0) {
        contributions[i] <- 0
      } else {
        val <- safe_funcs[[i]](model_outputs[[i]])
        if (!is.finite(val) || length(val) != 1) val <- -1e10
        val <- val * as.numeric(lhd_weight[i])
        contributions[i] <- val
      }
    }
    
    # Make sure notifs are bounded by bac confirmed notifs and all notifs for all TB
    if (!is.finite(sim_notifs) || sim_notifs < data_notif_bounds$notif_l || sim_notifs > data_notif_bounds$notif_u) {
      lhd <- -1e10
    } else {
      lhd <- sum(contributions)
    }
    
    return(list(lhd = lhd, incd_old = sim_incd_old_all, incd_new = sim_incd_new_all,
                incd_rr = sim_incd_rr, tb_mort = sim_tb_mort,
                tb_mort_ds = sim_tb_mort_ds, tb_mort_rr = sim_tb_mort_rr,
                notifs = sim_notifs, prev = sim_prev, 
                notifs_rr = sim_notifs_rr, 
                incd_h1hart = sim_incd_h1hart, tb_mort_h0 = sim_tb_mort_h0, tb_mort_h1hart = sim_tb_mort_h1hart,
                hiv_prev = sim_HIV_prev, ART_covg = sim_ART_covg,
                sfin = sfin))
    
  }
}


