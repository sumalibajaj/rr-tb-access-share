# Plot rr burden, and SL outcomes improvements

plot_combined_improvements <- function(
    tx_outcomes_all,
    filename = "data/TB_burden_countries_2025-11-27.csv",
    filename_rr = "data/MDR_RR_TB_burden_estimates_2025-11-27.csv",
    r_Tx2_fixed = 0.5,
    alpha = 0.5
) {
  
  # ==========================================
  # 1. RR-TB BURDEN: TOP 10 COUNTRIES
  # ==========================================
  
  pop <- read.csv(filename)
  rr_inc <- read.csv(filename_rr)
  
  pop_24 <- pop %>%
    dplyr::filter(year == 2024) %>%
    dplyr::select(country, e_pop_num, e_inc_100k)
  
  rr_inc_24 <- rr_inc %>%
    dplyr::filter(year == 2024) %>%
    dplyr::select(country, e_inc_rr_num)
  
  # Calculate global RR-TB burden before selecting top 10
  rr_all <- rr_inc_24 %>%
    dplyr::filter(!is.na(e_inc_rr_num)) %>%
    dplyr::arrange(dplyr::desc(e_inc_rr_num)) %>%
    dplyr::mutate(
      total_rr_inc = sum(e_inc_rr_num),
      rr_inc_prop = e_inc_rr_num / total_rr_inc,
      cumulative_prop = cumsum(e_inc_rr_num) / total_rr_inc
    )
  
  # Select top 10 countries
  rr_top10 <- rr_all %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::left_join(pop_24, by = "country") %>%
    dplyr::mutate(
      all_tb_cases = e_inc_100k * e_pop_num / 1e5,
      rr_tb_prop = e_inc_rr_num / all_tb_cases
    )
  
  # Country order based on RR-TB burden
  countries <- rr_top10$country
  
  rr_top10 <- rr_top10 %>%
    dplyr::mutate(
      country = factor(country, levels = countries)
    )
  
  # ==========================================
  # 2. SL TREATMENT OUTCOME IMPROVEMENTS
  # ==========================================
  
  tx_results <- purrr::map_dfr(countries, function(ctry) {
    
    mat <- tx_outcomes_all[[ctry]]
    
    if (is.null(mat)) {
      stop(paste("Missing treatment data:", ctry))
    }
    
    mdr_counts <- as.numeric(
      mat["SL", c("Success", "Failed", "Died", "LTFU")]
    )
    
    if (anyNA(mdr_counts) ||
        sum(mdr_counts) <= 0 ||
        sum(mdr_counts[1:2]) <= 0) {
      stop(paste("Invalid SL counts:", ctry))
    }
    
    tmp1 <- mdr_counts / sum(mdr_counts)
    
    tmp2_completed <- sum(tmp1[1:2])
    tmp2_died <- tmp1[3]
    tmp2_ltfu <- tmp1[4]
    
    p_cure <- tmp1[1] / tmp2_completed
    
    # Current rates
    rTx_sl <- r_Tx2_fixed
    
    muTx_sl_current <-
      (tmp2_died / tmp2_completed) * rTx_sl
    
    ltfu_sl_current <-
      (tmp2_ltfu / tmp2_completed) * rTx_sl
    
    vec_current <- c(
      rTx_sl,
      ltfu_sl_current,
      muTx_sl_current
    )
    
    propvec <- vec_current / sum(vec_current)
    
    # Same improvement equations as original
    newvec <- c(
      propvec[1] + alpha * (1 - propvec[1]),
      (1 - alpha) * propvec[2],
      (1 - alpha) * propvec[3]
    )
    
    newvec <- newvec / sum(newvec)
    
    cure_sl_new <- (p_cure + 1) * 0.5
    
    tibble::tibble(
      country = ctry,
      Died_current = propvec[3],
      Died_new = newvec[3],
      Failure_current =
        propvec[1] * (1 - p_cure),
      Failure_new =
        newvec[1] * (1 - cure_sl_new),
      LTFU_current = propvec[2],
      LTFU_new = newvec[2]
    )
  })
  
  tx_plot_df <- tx_results %>%
    tidyr::pivot_longer(
      cols = -country,
      names_to = c("outcome", "scenario"),
      names_pattern =
        "^(Died|Failure|LTFU)_(current|new)$",
      values_to = "prop"
    ) %>%
    tidyr::pivot_wider(
      names_from = scenario,
      values_from = prop
    ) %>%
    dplyr::mutate(
      country = factor(
        country,
        levels = rev(countries)
      ),
      outcome = factor(
        outcome,
        levels = c("Died", "Failure", "LTFU")
      )
    )
  
  # ==========================================
  # 3. COMMON STYLING
  # ==========================================
  
  current_col <- "#DDA448"
  new_col_sl <- "#14591D"
  
  common_theme <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
  
  sl_colors <- c(
    "Current" = current_col,
    "Improved" = new_col_sl
  )
  
  # ==========================================
  # 4. PANEL A: ABSOLUTE RR-TB BURDEN
  # ==========================================
  
  rr_top10 <- rr_top10 %>%
    dplyr::arrange(dplyr::desc(e_inc_rr_num)) %>%
    dplyr::mutate(
      country = factor(
        as.character(country),
        levels = as.character(country)
      ),
      cumulative_prop = cumsum(e_inc_rr_num) /
        dplyr::first(total_rr_inc)
    )
  
  max_inc <- max(rr_top10$e_inc_rr_num, na.rm = TRUE)
  scale_factor <- max_inc
  
  p_rr <- ggplot2::ggplot(
    rr_top10,
    ggplot2::aes(x = country)
  ) +
    
    ggplot2::geom_col(
      ggplot2::aes(y = e_inc_rr_num),
      fill = "#7B9EA8",
      width = 0.7
    ) +
    
    ggplot2::geom_line(
      ggplot2::aes(
        y = cumulative_prop * scale_factor,
        group = 1
      ),
      colour = "#5D3A00",
      linewidth = 1.1
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(y = cumulative_prop * scale_factor),
      colour = "#5D3A00",
      size = 2.5
    ) +
    
    ggplot2::geom_text(
      ggplot2::aes(
        y = cumulative_prop * scale_factor,
        label = scales::percent(
          cumulative_prop,
          accuracy = 1
        )
      ),
      vjust = -0.8,
      size = 3
    ) +
    
    ggplot2::scale_y_continuous(
      name = "RR-TB cases",
      labels = scales::label_comma(),
      sec.axis = ggplot2::sec_axis(
        ~ . / scale_factor,
        name = "Cumulative proportion of global burden",
        labels = scales::label_percent(accuracy = 1)
      ),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    
    ggplot2::labs(x = NULL) +
    
    common_theme +
    
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 9
      ),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank()
    )
  
  # ==========================================
  # 5. PANEL B: RR-TB PROPORTION OF ALL TB
  # ==========================================
  
  p_prop <- ggplot2::ggplot(
    rr_top10,
    ggplot2::aes(x = country, y = rr_tb_prop)
  ) +
    
    ggplot2::geom_col(
      fill = "#7B9EA8",
      width = 0.7
    ) +
    
    ggplot2::scale_y_continuous(
      name = "RR-TB as a proportion of all TB (%)",
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    
    ggplot2::labs(x = NULL) +
    
    common_theme +
    
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 9
      ),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank()
    )
  
  # ==========================================
  # 6. SEPARATE SL TREATMENT OUTCOMES PLOT
  # ==========================================
  
  p_tx <- ggplot2::ggplot(
    tx_plot_df,
    ggplot2::aes(y = country)
  ) +
    
    ggplot2::geom_segment(
      ggplot2::aes(
        x = current,
        xend = new,
        yend = country
      ),
      colour = "grey65",
      linewidth = 0.8
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        x = current,
        colour = "Current"
      ),
      size = 2.7
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        x = new,
        colour = "Improved"
      ),
      size = 2.7
    ) +
    
    ggplot2::facet_wrap(
      ~ outcome,
      nrow = 1,
      scales = "free_x"
    ) +
    
    ggplot2::scale_colour_manual(
      values = sl_colors,
      breaks = c("Current", "Improved")
    ) +
    
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0.05, 0.1))
    ) +
    
    ggplot2::labs(
      x = "SL treatment outcome (%)",
      y = NULL
    ) +
    
    common_theme
  
  # ==========================================
  # 7. COMBINE RR-TB BURDEN PANELS A AND B
  # ==========================================
  
  combined_plot <- cowplot::plot_grid(
    p_rr,
    p_prop,
    labels = c("A", "B"),
    ncol = 1,
    rel_heights = c(1, 1),
    align = "v",
    axis = "l"
  )
  
  # ==========================================
  # 8. RETURN TWO SEPARATE FIGURES
  # ==========================================
  
  return(list(
    burden_plot = combined_plot,
    sl_plot = p_tx,
    rr_data = rr_top10,
    tx_data = tx_plot_df,
    rr_plot = p_rr,
    proportion_plot = p_prop
  ))
}


plot_inc <- function(filename = "data/TB_burden_countries_2025-11-27.csv", 
                     filename_rr = "data/MDR_RR_TB_burden_estimates_2025-11-27.csv", 
                     countries) {
  
  df <- read.csv(filename)
  df_rr <- read.csv(filename_rr)
  
  df_plot <- df %>%
    filter(country %in% countries) %>%
    arrange(country, year) %>%
    dplyr::select(country, year, e_inc_100k, e_inc_100k_lo, e_inc_100k_hi)
  
  df_rr_plot <- df_rr %>%
    filter(country %in% countries) %>%
    merge(df %>% dplyr::select(country, year, e_pop_num)) %>%
    mutate(e_inc_rr_100k = e_inc_rr_num/e_pop_num*100000,
           e_inc_rr_100k_lo = e_inc_rr_num_lo/e_pop_num*100000,
           e_inc_rr_100k_hi = e_inc_rr_num_hi/e_pop_num*100000) %>%
    dplyr::select(country, year, e_inc_rr_100k, e_inc_rr_100k_lo, e_inc_rr_100k_hi)
  
  df_plot_combined <- left_join(df_plot, df_rr_plot)
  
  p <- ggplot(df_plot_combined, aes(x = year)) +
    geom_ribbon(aes(ymin = e_inc_100k_lo, ymax = e_inc_100k_hi), 
                alpha = 0.2) +
    geom_line(aes(y = e_inc_100k)) +
    geom_ribbon(aes(ymin = e_inc_rr_100k_lo, ymax = e_inc_rr_100k_hi),
                alpha = 0.2, fill = "blue") +
    geom_line(aes(y = e_inc_rr_100k), color = "blue") +
    facet_wrap(~ country, scales = "free_y") +
    labs(x = "Year",
         y = "Incidence per 100k") +
    theme_minimal()
  
  
  p_rr <- ggplot(df_plot_combined %>% na.omit(), aes(x = year)) +
    geom_ribbon(aes(ymin = e_inc_rr_100k_lo, ymax = e_inc_rr_100k_hi),
                alpha = 0.2, fill = "blue") +
    geom_line(aes(y = e_inc_rr_100k), color = "blue") +
    facet_wrap(~ country, scales = "free_y") +
    labs(x = "Year",
         y = "Incidence per 100k") +
    theme_minimal()
  p_rr
  return(list(p = p, p_rr = p_rr))
}


plot_hiv_incidence <- function(hiv_sa, hiv_ua) {
  
  hiv_df <- dplyr::bind_rows(
    tibble::tibble(
      country = "South Africa",
      year = hiv_sa$HIV_years,
      hiv_incidence = hiv_sa$rHIV
    ),
    tibble::tibble(
      country = "Ukraine",
      year = hiv_ua$HIV_years,
      hiv_incidence = hiv_ua$rHIV
    )
  )
  
  p <- ggplot2::ggplot(
    hiv_df,
    ggplot2::aes(
      x = year,
      y = hiv_incidence,
      colour = country
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(
      ~ country,
      nrow = 1,
      scales = "free_y"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "South Africa" = "#a63a50",
        "Ukraine" = "#235789"
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(1980, 2035, 10)
    ) +
    ggplot2::labs(
      x = "Year",
      y = "Annual HIV incidence rate"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(
        accuracy = 0.0001)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
  
  
  return(p)
}


