allocate_parameters_m1_1 <- function(theta, thetai, p, r) {
  r[["beta"]]           <- theta[thetai[["r_beta"]]]
  r[["Tx_init_new"]]    <- theta[thetai[["r_Tx_init_new"]]]
  p[["pu"]]             <- theta[thetai[["p_pu"]]]
  p[["pr_en"]]          <- theta[thetai[["p_pr_en"]]]
  
  r[["Tx_init_old"]]    <- r[["Tx_init_new"]]*theta[thetai[["rf_Tx_init_old"]]]   
  r[["mort_TB"]]        <- r[["mort_TB0"]]*theta[thetai[["rf_mort_TB"]]] 
  r[["self_cure"]]      <- r[["self_cure0"]]*theta[thetai[["rf_self_cure"]]]  
  r[["progression"]]    <- r[["progression0"]]*theta[thetai[["rf_progression"]]]
  r[["stabilisation"]]  <- r[["stabilisation0"]]*theta[thetai[["rf_stabilisation"]]]
  r[["reactivation"]]   <- r[["reactivation0"]]*theta[thetai[["rf_reactivation"]]]
  r[["relapse"]]        <- r[["relapse0"]]*theta[thetai[["rf_relapse"]]] 
  
  p[["imm"]]            <- theta[thetai[["p_imm"]]]
  
  p[["MDR_rec"]]        <- theta[thetai[["p_MDR_rec"]]]
  
  return(list(p = p, r = r))
}

allocate_parameters_m2_1 <- function(theta, thetai, p, r) {
  r[["beta"]]           <- theta[thetai[["r_beta"]]]
  r[["Tx_init_new"]]    <- theta[thetai[["r_Tx_init_new"]]]
  p[["pu"]]             <- theta[thetai[["p_pu"]]]
  
  r[["Tx_init_old"]]    <- r[["Tx_init_new"]]*theta[thetai[["rf_Tx_init_old"]]]   
  r[["mort_TB"]]        <- r[["mort_TB0"]]*theta[thetai[["rf_mort_TB"]]] 
  r[["self_cure"]]      <- r[["self_cure0"]]*theta[thetai[["rf_self_cure"]]]  
  r[["progression"]]    <- r[["progression0"]]*theta[thetai[["rf_progression"]]]
  r[["stabilisation"]]  <- r[["stabilisation0"]]*theta[thetai[["rf_stabilisation"]]]
  r[["reactivation"]]   <- r[["reactivation0"]]*theta[thetai[["rf_reactivation"]]]
  r[["relapse"]]        <- r[["relapse0"]]*theta[thetai[["rf_relapse"]]] 
  
  p[["imm"]]            <- theta[thetai[["p_imm"]]]
  
  p[["MDR_rec"]]        <- theta[thetai[["p_MDR_rec"]]]
  p[["MDR_rec_pr_nen"]]        <- theta[thetai[["p_MDR_rec_pr_nen"]]]
  p[["SL_trans_pr_nen"]]        <- theta[thetai[["p_SL_trans_pr_nen"]]]
  
  return(list(p = p, r = r))
}

allocate_parameters_m3_1 <- function(theta, thetai, p, r) {
  r[["beta"]]           <- theta[thetai[["r_beta"]]]
  r[["Tx_init_new"]]    <- theta[thetai[["r_Tx_init_new"]]]
  
  r[["Tx_init_old"]]    <- r[["Tx_init_new"]]*theta[thetai[["rf_Tx_init_old"]]]
  r[["mort_TB"]]        <- r[["mort_TB0"]]*theta[thetai[["rf_mort_TB"]]] 
  r[["self_cure"]]      <- r[["self_cure0"]]*theta[thetai[["rf_self_cure"]]]  
  r[["progression"]]    <- r[["progression0"]]*theta[thetai[["rf_progression"]]]
  r[["stabilisation"]]  <- r[["stabilisation0"]]*theta[thetai[["rf_stabilisation"]]]
  r[["reactivation"]]   <- r[["reactivation0"]]*theta[thetai[["rf_reactivation"]]]
  r[["relapse"]]        <- r[["relapse0"]]*theta[thetai[["rf_relapse"]]] 
  
  p[["imm"]]            <- theta[thetai[["p_imm"]]]
  
  p[["MDR_rec"]]        <- theta[thetai[["p_MDR_rec"]]]
  
  p[["mort_Tx_TB"]]     <- theta[thetai[["p_mort_Tx_TB"]]]
  
  return(list(p = p, r = r))
}

allocate_parameters_m4_1 <- function(theta, thetai, p, r) {
  theta <- as.numeric(theta)
  r[["beta"]]           <- theta[thetai[["r_beta"]]]
  r[["Tx_init_new"]]    <- theta[thetai[["r_Tx_init_new"]]]
  
  r[["Tx_init_old"]]    <- r[["Tx_init_new"]]*theta[thetai[["rf_Tx_init_old"]]]
  r[["mort_TB"]]        <- r[["mort_TB0"]]*theta[thetai[["rf_mort_TB"]]] 
  r[["self_cure"]]      <- r[["self_cure0"]]*theta[thetai[["rf_self_cure"]]]  
  r[["progression"]]    <- r[["progression0"]]*theta[thetai[["rf_progression"]]]
  r[["stabilisation"]]  <- r[["stabilisation0"]]*theta[thetai[["rf_stabilisation"]]]
  r[["reactivation"]]   <- r[["reactivation0"]]*theta[thetai[["rf_reactivation"]]]
  r[["relapse"]]        <- r[["relapse0"]]*theta[thetai[["rf_relapse"]]] 
  
  p[["MDR_rec"]]        <- theta[thetai[["p_MDR_rec"]]]
  
  # HIV related parameters
  r[["ART_init"]]       <- theta[thetai[["r_ART_init"]]]
  r[["HIV_mort"]]       <- r[["mort"]]*theta[thetai[["rf_HIV_mort"]]] 
  rel                   <- theta[thetai[["p_HIV_relrate"]]]
  
  # HIV-adjusted natural history (order- h0, h1, hart)
  r[["self_cure"]]      <- r[["self_cure"]]*c(1, 0, 1)
  r[["progression"]]    <- r[["progression"]]*c(1, rel, rel*0.4)
  r[["stabilisation"]]  <- r[["stabilisation"]]*c(1, 0, 1)
  r[["reactivation"]]   <- r[["reactivation"]]*c(1, rel, rel*0.4)
  
  p_imm                 <- theta[thetai[["p_imm"]]]
  p[["imm"]]            <- c(p_imm, 0, p_imm)
  
  return(list(p = p, r = r))
}

