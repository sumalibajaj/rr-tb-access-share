# Extract country data relations
extract_country_data_relations <- function(
    country_names,
    prct_pulm_bac_conf_df,
    filename_burden = "data/TB_burden_countries_2025-11-27.csv",
    filename_rr = "data/MDR_RR_TB_burden_estimates_2025-11-27.csv",
    filename_notifs = "data/TB_notifications_2026-04-09.csv"
) {
  
  results <- lapply(country_names, function(ctry) {
    
    # Get country-specific proportion
    prct <- prct_pulm_bac_conf_df |>
      dplyr::filter(country == ctry) |>
      dplyr::pull(prct_pulm_bac_conf)
    
    if (length(prct) != 1 || is.na(prct) ||
        prct < 0 || prct > 1) {
      stop("Invalid or missing prct_pulm_bac_conf for: ", ctry)
    }
    
    # Extract calibration data
    dat <- extract_calib_data(
      filename_burden = filename_burden,
      filename_rr = filename_rr,
      filename_notifs = filename_notifs,
      country_name = ctry,
      old_inc_year = 2023, # can ignore old_inc data
      prct_pulm_bac_conf = prct
    )
    
    # Central estimates
    inc_new <- dat$incd_new[2]
    inc_rr  <- dat$incd_rr[2]
    mort    <- dat$tb_mort[2]
    noti    <- dat$notifs_pu[2]
    noti_rr <- dat$notifs_rr_pu[2]
    
    data.frame(
      country = ctry,
      prct_pulm_bac_conf = prct,
      
      # Original calibration data
      inc_total = inc_new,
      inc_rr = inc_rr,
      tb_mort = mort,
      notifs = noti,
      notifs_rr = noti_rr,
      
      # Derived indicators
      rr_inc_fraction = inc_rr / inc_new,
      rr_noti_inc_ratio = noti_rr / inc_rr,
      total_noti_inc_ratio = noti / inc_new,
      mort_inc_ratio = mort / inc_new
    )
  })
  
  dplyr::bind_rows(results)
}


# Extract DST uptake additional percentage
extract_dst_limit <- function(country_names,
                              filename_notifs = "data/TB_notifications_2026-04-09.csv") {
  
  df_notifs <- read.csv(filename_notifs)
  
  results <- df_notifs %>%
    dplyr::filter(country %in% country_names, year == 2024) %>%
    dplyr::mutate(pulm_tb_notif = c_newinc - new_ep,
                  clindx = new_clindx + ret_rel_clindx,
                  clindx_by_pulm_tb_notif = clindx/pulm_tb_notif,
                  newinc_pulm_clindx_rdx_by_clindx = newinc_pulm_clindx_rdx/clindx,
                  dst_limit = clindx_by_pulm_tb_notif* (1 - newinc_pulm_clindx_rdx_by_clindx)) %>%
    dplyr::select(country, pulm_tb_notif, clindx, newinc_pulm_clindx_rdx, 
                  clindx_by_pulm_tb_notif, newinc_pulm_clindx_rdx_by_clindx, 
                  dst_limit) %>%
    dplyr::rename(
      pulmonary_tb_notifications = pulm_tb_notif,
      clinically_diagnosed_tb = clindx,
      clinically_diagnosed_with_test = newinc_pulm_clindx_rdx,
      prop_pulmonary_tb_clinically_diagnosed = clindx_by_pulm_tb_notif,
      prop_clinically_diagnosed_with_test = newinc_pulm_clindx_rdx_by_clindx,
      dst_limit = dst_limit
    )
  
  return(results)
}


# Function to extract treatment outcomes
extract_tx_outcomes <- function(filename = "data/TB_outcomes_2026-02-20.csv", 
                                country_names) {
  
  df <- read.csv(filename)
  
  outcome_cols <- list(FL  = c("newrel_succ", "newrel_fail", "newrel_died", "newrel_lost"),
                       SL  = c("mdr_succ", "mdr_fail", "mdr_died", "mdr_lost"),
                       HIV = c("tbhiv_succ", "tbhiv_fail", "tbhiv_died", "tbhiv_lost"))
  
  # Check columns exist
  required <- c("country", "year", unlist(outcome_cols))
  missing <- setdiff(required, names(df))
  
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }
  
  results <- lapply(country_names, function(ctry) {
    
    df_country <- df %>%
      dplyr::filter(country == ctry,
                    year <= 2024)
    
    Tx_outc <- lapply(names(outcome_cols), function(type) {
      
      cols <- outcome_cols[[type]]
      
      # Allow missing HIV columns
      if (!all(cols %in% names(df_country))) {
        if (type == "HIV") return(rep(NA_real_, 4))
        stop("Missing ", type, " columns for ", ctry)
      }
      
      # Select latest year with all four outcomes available
      valid <- df_country %>%
        dplyr::filter(dplyr::if_all(dplyr::all_of(cols), ~ !is.na(.x))) %>%
        dplyr::arrange(dplyr::desc(year))
      
      if (nrow(valid) == 0) {
        if (type == "HIV") return(rep(NA_real_, 4))
        stop("No ", type, " outcomes for ", ctry)
      }
      
      as.numeric(valid[1, cols])
    })
    
    Tx_outc <- do.call(rbind, Tx_outc)
    
    rownames(Tx_outc) <- c("FL", "SL", "HIV")
    colnames(Tx_outc) <- c("Success", "Failed", "Died", "LTFU")
    
    Tx_outc
  })
  
  names(results) <- country_names
  
  return(results)
}


# Extract calibration data for a country
extract_calib_data <- function(filename_burden = "data/TB_burden_countries_2025-11-27.csv", 
                               filename_rr = "data/MDR_RR_TB_burden_estimates_2025-11-27.csv", 
                               filename_notifs = "data/TB_notifications_2026-04-09.csv", 
                               country_name, 
                               old_inc_year, 
                               prct_pulm_bac_conf){
  
  # All TB incidence data
  df_burden <- read.csv(filename_burden) %>%
    dplyr::filter(country == country_name)
  
  incd_old <- df_burden %>% 
    dplyr::filter(year == old_inc_year) %>%
    dplyr::select(e_inc_100k_lo, e_inc_100k, e_inc_100k_hi) %>% 
    as.numeric()
  incd_new <- df_burden %>% 
    dplyr::filter(year == 2024) %>%
    dplyr::select(e_inc_100k_lo, e_inc_100k, e_inc_100k_hi) %>% 
    as.numeric()
  
  # All TB mortality data
  tb_mort <- df_burden %>%
    dplyr::filter(year == 2024) %>%
    dplyr::select(e_mort_100k_lo, e_mort_100k, e_mort_100k_hi) %>% 
    as.numeric()
  
  pop_2024 = df_burden %>% 
    dplyr::filter(year == 2024) %>% 
    dplyr::pull(e_pop_num)
  
  # RR-TB incidence data
  df_rr <- read.csv(filename_rr) %>%
    dplyr::filter(country == country_name)
  
  incdRR_num <- df_rr %>% 
    dplyr::filter(year == 2024) %>%
    dplyr::select(e_inc_rr_num_lo, e_inc_rr_num, e_inc_rr_num_hi) %>% 
    as.numeric()
  incd_rr <- (incdRR_num/pop_2024)*1e5
  
  # RR-TB notifications public sector (public + engaged private)
  df_notifs <- read.csv(filename_notifs) %>%
    dplyr::filter(country == country_name)
  
  notifs_rr_pu_num <-  df_notifs %>% 
    dplyr::filter(year == 2024) %>%
    dplyr::select(conf_rr_nfqr) %>% 
    as.numeric()
  notifs_rr_pu <- (notifs_rr_pu_num/pop_2024)*1e5*c(0.9, 1, 1.1)
  
  # All TB notifications public sector (public + engaged private)
  notifs_pu_num_hi <- df_notifs %>% 
    dplyr::filter(year == 2024) %>%
    dplyr::select(c_newinc) %>% 
    as.numeric()
  notifs_pu_hi <- (notifs_pu_num_hi/pop_2024)*1e5
  notifs_pu_lo <- prct_pulm_bac_conf*notifs_pu_hi
  notifs_pu <- c(notifs_pu_lo, (notifs_pu_lo + notifs_pu_hi)/2, notifs_pu_hi)
  
  return(list(incd_new = incd_new,
              incd_old = incd_old,
              incd_rr = incd_rr,
              tb_mort = tb_mort,
              notifs_pu = notifs_pu,
              notifs_rr_pu = notifs_rr_pu))
  
}

# Extract HIV calibration data for a country
extract_hiv_calib_data <- function(filename_burden = "data/TB_burden_countries_2025-11-27.csv", 
                                   country_name){
  
  # All TB incidence data in HIV positive people
  df_burden <- read.csv(filename_burden) %>%
    dplyr::filter(country == country_name)
  
  incd_h1hart <- df_burden %>% 
    dplyr::filter(year == 2024) %>%
    dplyr::select(e_inc_tbhiv_100k_lo, e_inc_tbhiv_100k, e_inc_tbhiv_100k_hi) %>% 
    as.numeric()
  
  # All TB mortality data in HIV negative people
  tb_mort_h0 <- df_burden %>%
    dplyr::filter(year == 2024) %>%
    dplyr::select(e_mort_exc_tbhiv_100k_lo, e_mort_exc_tbhiv_100k, e_mort_exc_tbhiv_100k_hi) %>% 
    as.numeric()
  
  # All TB mortality data in HIV positive people (on and not on ART)
  tb_mort_h1hart <- df_burden %>%
    dplyr::filter(year == 2024) %>%
    dplyr::select(e_mort_tbhiv_100k_lo, e_mort_tbhiv_100k, e_mort_tbhiv_100k_hi) %>% 
    as.numeric()
  
  return(list(incd_h1hart = incd_h1hart,
              tb_mort_h0 = tb_mort_h0,
              tb_mort_h1hart = tb_mort_h1hart))
  
}
