get_hiv_data <- function(filename, country_name) {
  df <- read.csv(filename)
  
  # dat <- df %>%
  #   filter(Area == country_name) %>%
  #   filter(Indicator %in% c("HIV Incidence per 1000 population", "HIV Prevalence", "Coverage of people living with HIV receiving ART")) %>%
  #   filter(Subgroup %in% c("All ages estimate", "All ages lower estimate", "All ages upper estimate", "Adults (15+) estimate")) %>%
  #   filter(`Time.Period` == 2024)
  
  dat <- df %>%
    filter(Area == country_name)
  
  pick <- function(ind, year = NULL) {
    x <- dat %>% filter(Indicator_GId == ind)
    
    if (!is.null(year)) {
      x <- x %>% filter(`Time.Period` == year)
    }
    
    y <- x %>% filter(Subgroup_Val_GId == "ESTIMATE_ALL_AGES")
    if (nrow(y) == 0) y <- x %>% filter(Subgroup_Val_GId == "ESTIMATE_15_AND_UP")
    y
  }
  
  pick_ci <- function(ind, year = NULL) {
    x <- dat %>% filter(Indicator_GId == ind)
    
    if (!is.null(year)) {
      x <- x %>% filter(`Time.Period` == year)
    }
    
    lab <- if ("Subgroup" %in% names(x)) x$Subgroup else x$Subgroup_Val_GId
    
    lo <- x$Data.value[grepl("lower", lab, ignore.case = TRUE)][1]
    mid <- x$Data.value[grepl("estimate", lab, ignore.case = TRUE) &
                          !grepl("lower|upper", lab, ignore.case = TRUE)][1]
    hi <- x$Data.value[grepl("upper", lab, ignore.case = TRUE)][1]
    
    c(lo, mid, hi)
  }
  
  inc <- pick("HIV_INCIDENCE") %>%
    transmute(year = as.integer(`Time.Period`),
              value = Data.value / 1000) %>%
    arrange(year)
  
  hiv_years <- 1980:2035
  
  # observed years and values
  x <- inc$year
  y <- inc$value
  
  # start with NA
  rHIV <- rep(NA_real_, length(hiv_years))
  
  # 1) ramp from 1980 to first observed year
  first_year <- min(x)
  first_val  <- y[which.min(x)]
  
  pre <- hiv_years < first_year
  rHIV[pre] <- first_val * (hiv_years[pre] - 1980) / (first_year - 1980)
  
  # 2) interpolate within observed years
  mid <- hiv_years >= first_year & hiv_years <= max(x)
  rHIV[mid] <- approx(x = x, y = y, xout = hiv_years[mid], rule = 1)$y
  
  # 3) linear extrapolation after last observed year
  last2 <- tail(inc, 2)
  m <- diff(last2$value) / diff(last2$year)
  
  future <- hiv_years > max(x)
  rHIV[future] <- last2$value[2] + m * (hiv_years[future] - last2$year[2])
  
  # 4) never allow negative incidence
  rHIV <- pmax(rHIV, 0)
  
  # HIV prevalence in 2024
  HIVprev_2024 <- pick_ci("HIV_PREVALENCE", 2024) / 100
  
  # ART coverage in 2024 and start year
  art <- pick("PERCENT_PLWH") %>%
    transmute(year = as.integer(`Time.Period`),
              value = Data.value / 100) %>%
    arrange(year)
  
  ARTcovg_2024 <- pick_ci("PERCENT_PLWH", 2024) / 100
  
  art_nz <- art %>%
    filter(value > 0)
  
  xp <- c(first(art_nz$year), last(art_nz$year))
  yp <- c(first(art_nz$value), last(art_nz$value))
  
  m <- (yp[2] - yp[1]) / (xp[2] - xp[1])
  b <- yp[1] - m * xp[1]
  
  ART_start <- floor(-b / m)
  
  list(
    rHIV = rHIV,
    HIVprev_2024 = HIVprev_2024,
    ARTcovg_2024 = ARTcovg_2024,
    ART_start = ART_start,
    HIV_years = hiv_years
  )
  
}


# Using data from Thembisa model
library(readxl)

get_num <- function(x) {
  x <- as.character(x)
  x <- sub("\\s*\\[.*$", "", x)   # keep estimate before CI
  x <- gsub("[^0-9]", "", x)      # remove spaces, <, commas, etc.
  ifelse(x == "", NA_real_, as.numeric(x))
}

get_row <- function(file, country_name) {
  df <- read_xlsx(file)
  row <- df[df[[1]] == country_name, ]
  if (nrow(row) == 0) stop("Country not found: ", country_name)
  row
}


get_num_ci <- function(x) {
  x <- as.character(x)
  x <- sub("\\s*\\(.*$", "", x)   # remove trailing notes
  
  est <- sub("\\s*\\[.*$", "", x)
  rng <- sub(".*\\[(.*)\\].*", "\\1", x)
  
  lo <- sub("^(.*) - .*", "\\1", rng)
  hi <- sub("^.* - (.*)$", "\\1", rng)
  
  num <- function(z) {
    z <- gsub("[^0-9]", "", z)
    ifelse(z == "", NA_real_, as.numeric(z))
  }
  
  c(estimate = num(est), lower = num(lo), upper = num(hi))
}

extend_linear <- function(x, y, xout) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  
  out <- approx(x, y, xout = xout, rule = 1)$y
  
  if (length(x) >= 2) {
    if (min(xout) < min(x)) {
      first2 <- head(data.frame(x = x, y = y), 2)
      m <- diff(first2$y) / diff(first2$x)
      pre <- xout < min(x)
      out[pre] <- first2$y[1] + m * (xout[pre] - first2$x[1])
    }
    
    if (max(xout) > max(x)) {
      last2 <- tail(data.frame(x = x, y = y), 2)
      m <- diff(last2$y) / diff(last2$x)
      post <- xout > max(x)
      out[post] <- last2$y[2] + m * (xout[post] - last2$x[2])
    }
  }
  
  out
}

build_series <- function(x, y, years, start_year = 1980, start_zero = FALSE) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  
  out <- rep(NA_real_, length(years))
  if (length(x) == 0) return(out)
  
  first_year <- min(x)
  last_year  <- max(x)
  
  mid <- years >= first_year & years <= last_year
  out[mid] <- approx(x, y, xout = years[mid])$y
  
  pre <- years < first_year
  if (any(pre)) {
    first_val <- y[which.min(x)]
    if (start_zero) {
      out[pre] <- first_val * (years[pre] - start_year) / (first_year - start_year)
    } else {
      out[pre] <- first_val * (years[pre] - start_year) / (first_year - start_year)
    }
  }
  
  post <- years > last_year
  if (any(post) && length(x) >= 2) {
    last2 <- tail(data.frame(x = x, y = y), 2)
    m <- diff(last2$y) / diff(last2$x)
    out[post] <- last2$y[2] + m * (years[post] - last2$x[2])
  }
  
  pmax(out, 0)
}

get_hiv_data_thembisa <- function(country_name = "Ukraine",
                                  infections_file = "data/hiv_thembisa/New_infections.xlsx",
                                  prevalence_file = "data/hiv_thembisa/HIV prevalence.xlsx",
                                  art_file = "data/hiv_thembisa/ART coverage.xlsx",
                                  pop_file = "data/TB_burden_countries_2025-11-27.csv") {
  
  all_years <- 1980:2035
  
  # population
  pop <- read.csv(pop_file)
  pop <- pop[pop$country == country_name, c("year", "e_pop_num")]
  pop_full <- extend_linear(pop$year, pop$e_pop_num, all_years)
  
  # infections
  inf <- get_row(infections_file, country_name)
  inf_years <- as.integer(names(inf)[-1])
  inf_vals  <- sapply(inf[-1], get_num)
  inf_rate  <- (inf_vals / pop_full[match(inf_years, all_years)])
  rHIV <- build_series(inf_years, inf_rate, all_years, start_zero = TRUE)
  
  # prevalence with uncertainty intervals
  prev <- get_row(prevalence_file, country_name)
  prev_years <- as.integer(names(prev)[-1])
  prev_raw <- as.character(unlist(prev[-1], use.names = FALSE))
  
  prev_q <- t(sapply(prev[-1], get_num_ci))
  bad_prev <- grepl("^<", prev_raw)
  prev_q[bad_prev, ] <- NA_real_
  
  prev_rate_est <- prev_q[, "estimate"] / pop_full[match(prev_years, all_years)]
  prev_rate_lo  <- prev_q[, "lower"]    / pop_full[match(prev_years, all_years)]
  prev_rate_hi  <- prev_q[, "upper"]    / pop_full[match(prev_years, all_years)]
  
  last_idx <- tail(which(!is.na(prev_rate_est)), 1)
  
  HIVprev_year <- prev_years[last_idx]
  HIVprev <- c(prev_rate_lo[last_idx], prev_rate_est[last_idx], prev_rate_hi[last_idx]) %>% as.numeric()
  
  
  # ART coverage with uncertainty intervals
  art <- get_row(art_file, country_name)
  art_years <- as.integer(names(art)[-1])
  
  art_q <- t(sapply(art[-1], get_num_ci)) / 100
  
  ART_est <- build_series(art_years, art_q[, "estimate"], all_years)
  ART_lo  <- build_series(art_years, art_q[, "lower"],    all_years)
  ART_hi  <- build_series(art_years, art_q[, "upper"],    all_years)
  
  ART_est <- pmin(pmax(ART_est, 0), 1)
  ART_lo  <- pmin(pmax(ART_lo, 0), 1)
  ART_hi  <- pmin(pmax(ART_hi, 0), 1)
  
  ARTcovg_2024 <- c(ART_lo[all_years == 2024], ART_est[all_years == 2024], ART_hi[all_years == 2024])
  
  list(rHIV = rHIV,
       HIVprev = HIVprev,
       HIVprev_year = HIVprev_year,
       ARTcovg_2024 = ARTcovg_2024,
       HIV_years = all_years)
}
