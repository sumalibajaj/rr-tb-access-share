get_addresses <- function(groups, i = list(), s = list(), d = list(), lim = 0) {
  # Initialize fields in 's' if they are not already present
  if (length(s) > 0) {
    existing <- names(s)
  } else {
    existing <- character()
  }
  
  # Ensure each group name is initialized
  for (gp in groups) {
    for (key in gp) {
      if (!key %in% existing) {
        s[[key]] <- integer()
      }
    }
  }
  
  # Depending on how many groupings, we construct combinations
  if (length(groups) == 1) {
    gp1 <- groups[[1]]
    for (g1 in gp1) {
      lim <- lim + 1
      i[[g1]] <- lim
      s[[g1]] <- c(s[[g1]], lim)
      d[[lim]] <- g1
    }
  }
  
  if (length(groups) == 2) {
    gp1 <- groups[[1]]
    gp2 <- groups[[2]]
    for (g1 in gp1) {
      for (g2 in gp2) {
        lim <- lim + 1
        if (is.null(i[[g1]])) i[[g1]] <- list()
        i[[g1]][[g2]] <- lim
        s[[g1]] <- c(s[[g1]], lim)
        s[[g2]] <- c(s[[g2]], lim)
        d[[lim]] <- paste(g1, g2)
      }
    }
  }
  
  if (length(groups) == 3) {
    gp1 <- groups[[1]]
    gp2 <- groups[[2]]
    gp3 <- groups[[3]]
    for (g1 in gp1) {
      for (g2 in gp2) {
        for (g3 in gp3) {
          lim <- lim + 1
          if (is.null(i[[g1]])) i[[g1]] <- list()
          if (is.null(i[[g1]][[g2]])) i[[g1]][[g2]] <- list()
          i[[g1]][[g2]][[g3]] <- lim
          s[[g1]] <- c(s[[g1]], lim)
          s[[g2]] <- c(s[[g2]], lim)
          s[[g3]] <- c(s[[g3]], lim)
          d[[lim]] <- paste(g1, g2, g3)
        }
      }
    }
  }
  
  if (length(groups) == 4) {
    gp1 <- groups[[1]]
    gp2 <- groups[[2]]
    gp3 <- groups[[3]]
    gp4 <- groups[[4]]
    for (g1 in gp1) {
      for (g2 in gp2) {
        for (g3 in gp3) {
          for (g4 in gp4) {
            lim <- lim + 1
            if (is.null(i[[g1]])) i[[g1]] <- list()
            if (is.null(i[[g1]][[g2]])) i[[g1]][[g2]] <- list()
            if (is.null(i[[g1]][[g2]][[g3]])) i[[g1]][[g2]][[g3]] <- list()
            i[[g1]][[g2]][[g3]][[g4]] <- lim
            s[[g1]] <- c(s[[g1]], lim)
            s[[g2]] <- c(s[[g2]], lim)
            s[[g3]] <- c(s[[g3]], lim)
            s[[g4]] <- c(s[[g4]], lim)
            d[[lim]] <- paste(g1, g2, g3, g4)
          }
        }
      }
    }
  }
  
  i$nstates <- lim
  return(list(i = i, s = s, d = d, lim = lim))
}
