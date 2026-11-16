goveqs_m1_1 <- function(t, state, parameters) {
  r <- parameters[["r"]]
  p <- parameters[["p"]]
  
  # rates and proportions
  r_beta_ds              <- r[["beta"]][1]
  r_beta_rr              <- r[["beta"]][2]
  r_Tx_init              <- r[["Tx_init"]]
  r_progression          <- r[["progression"]]
  r_stabilisation        <- r[["stabilisation"]]
  r_reactivation         <- r[["reactivation"]]
  r_relapse_Rlo          <- r[["relapse"]][1]
  r_relapse_Rhi          <- r[["relapse"]][2]
  r_relapse_R            <- r[["relapse"]][3]
  r_self_cure            <- r[["self_cure"]]
  r_MDR_acqu             <- r[["MDR_acqu"]]
  r_Tx                   <- r[["Tx"]]
  r_Tx2                  <- r[["Tx2"]]
  r_ltfu1                <- r[["ltfu"]][1]
  r_ltfu2_pu             <- r[["ltfu"]][2]
  r_ltfu2_pr_en          <- r[["ltfu"]][3]
  r_ltfu2_pr_nen         <- r[["ltfu"]][4]
  
  p_pu                   <- p[["pu"]]
  p_pr_en                <- p[["pr_en"]]
  p_imm                  <- p[["imm"]]
  p_MDR_rec_pu           <- p[["MDR_rec"]][1]
  p_MDR_rec_pr_en        <- p[["MDR_rec"]][2]
  p_MDR_rec_pr_nen       <- p[["MDR_rec_pr_nen"]]
  p_SL_trans_pu          <- p[["SL_trans"]][1]
  p_SL_trans_pr_en       <- p[["SL_trans"]][1]
  p_SL_trans_pr_nen      <- p[["SL_trans"]][2]
  p_cure1                <- p[["cure"]][1]
  p_cure2_pu             <- p[["cure"]][2]
  p_cure2_pr_en          <- p[["cure"]][3]
  p_cure2_pr_nen         <- p[["cure"]][4]
  
  r_mort                 <- r[["mort"]]
  r_mort_TB              <- r[["mort_TB"]]
  r_muTx1                <- r[["muTx"]][1]
  r_muTx2_pu             <- r[["muTx"]][2]
  r_muTx2_pr_en          <- r[["muTx"]][3]
  r_muTx2_pr_nen         <- r[["muTx"]][4]
  
  # states
  U                      <- state[1]
  Lf_ds                  <- state[2]
  Lf_rr                  <- state[3]
  Ls_ds                  <- state[4]
  Ls_rr                  <- state[5]
  I_ds                   <- state[6]
  I_rr                   <- state[7]
  Rlo_ds                 <- state[8]
  Rlo_rr                 <- state[9]
  Rhi_ds                 <- state[10]
  Rhi_rr                 <- state[11]
  R_ds                   <- state[12]
  R_rr                   <- state[13]
  Tx_ds_pu               <- state[14]
  Tx_ds_pr_en            <- state[15]
  Tx_ds_pr_nen           <- state[16]
  Tx_rr_pu               <- state[17]
  Tx_rr_pr_en            <- state[18]
  Tx_rr_pr_nen           <- state[19]
  Tx2_ds_pu              <- state[20]
  Tx2_ds_pr_en           <- state[21]
  Tx2_ds_pr_nen          <- state[22]
  Tx2_rr_pu              <- state[23]
  Tx2_rr_pr_en           <- state[24]
  Tx2_rr_pr_nen          <- state[25]
  
  # auxiliaries
  cu_inc_ds              <- state[26]
  cu_inc_rr              <- state[27]
  cu_noti_pu             <- state[28]
  cu_mort_tb             <- state[29]
  cu_mort_tb_ds          <- state[30]
  cu_mort_tb_rr          <- state[31]
  cu_noti_pr_en          <- state[32]
  cu_noti_rr_pu          <- state[33]
  cu_noti_rr_pr_en       <- state[34]
  
  # compute some useful values
  N <- sum(state[1:25])
  lam_ds <- r_beta_ds*(I_ds)/N  # ds FOI
  lam_rr <- r_beta_rr*(I_rr + Tx_rr_pu + Tx_rr_pr_en + Tx_rr_pr_nen)/N  # rr FOI
  births <- (r_mort)*(U + Lf_ds + Lf_rr + Ls_ds + Ls_rr + 
                        Rlo_ds + Rlo_rr + Rhi_ds + Rhi_rr + 
                        R_ds + R_rr) + 
    (r_mort + r_mort_TB)*(I_ds + I_rr + 
                            Tx_rr_pu + Tx_rr_pr_en + Tx_rr_pr_nen) +
    (r_mort + r_muTx1)*(Tx_ds_pu + Tx_ds_pr_en + Tx_ds_pr_nen) +
    (r_mort + r_muTx2_pu)*Tx2_rr_pu + (r_mort + r_muTx2_pr_en)*Tx2_rr_pr_en + (r_mort + r_muTx2_pr_nen)*Tx2_rr_pr_nen
  
  # rates of change
  dU <- births - (lam_ds + lam_rr)*U - (r_mort)*U
  
  dLf_ds <- lam_ds*U + lam_ds*p_imm*(Ls_ds + Ls_rr + Rlo_ds + Rlo_rr + Rhi_ds + Rhi_rr + R_ds + R_rr) - 
    (r_stabilisation + r_progression)*Lf_ds - 
    (r_mort)*Lf_ds
  
  dLf_rr <- lam_rr*U + lam_rr*p_imm*(Ls_ds + Ls_rr + Rlo_ds + Rlo_rr + Rhi_ds + Rhi_rr + R_ds + R_rr) - 
    (r_stabilisation + r_progression)*Lf_rr - 
    (r_mort)*Lf_rr
  
  dLs_ds <- r_stabilisation*Lf_ds - 
    (r_reactivation + (lam_ds + lam_rr)*p_imm)*Ls_ds - 
    (r_mort)*Ls_ds
  
  dLs_rr <- r_stabilisation*Lf_rr - 
    (r_reactivation + (lam_ds + lam_rr)*p_imm)*Ls_rr - 
    (r_mort)*Ls_rr
  
  dI_ds <- (r_progression*Lf_ds + r_reactivation*Ls_ds + r_relapse_Rlo*Rlo_ds + r_relapse_Rhi*Rhi_ds + r_relapse_R*R_ds) + 
    r_Tx*(1-p_cure1)*(Tx_ds_pu + Tx_ds_pr_en + Tx_ds_pr_nen) - 
    (r_Tx_init + r_self_cure)*I_ds - 
    (r_mort + r_mort_TB)*I_ds
  
  dI_rr <- (r_progression*Lf_rr + r_reactivation*Ls_rr + r_relapse_Rlo*Rlo_rr + r_relapse_Rhi*Rhi_rr + r_relapse_R*R_rr) + 
    (r_Tx*(1-p_SL_trans_pu) + r_ltfu1)*Tx_rr_pu + (r_Tx*(1-p_SL_trans_pr_en) + r_ltfu1)*Tx_rr_pr_en + (r_Tx*(1-p_SL_trans_pr_nen) + r_ltfu1)*Tx_rr_pr_nen + 
    (r_Tx2*(1-p_cure2_pu))*Tx2_rr_pu + (r_Tx2*(1-p_cure2_pr_en))*Tx2_rr_pr_en + (r_Tx2*(1-p_cure2_pr_nen))*Tx2_rr_pr_nen - 
    (r_Tx_init + r_self_cure)*I_rr - 
    (r_mort + r_mort_TB)*I_rr
  
  dRlo_ds <- (r_Tx*p_cure1)*(Tx_ds_pu + Tx_ds_pr_en + Tx_ds_pr_nen) - 
    (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm)*Rlo_ds - 
    (r_mort)*Rlo_ds
  
  dRlo_rr <- (r_Tx2*p_cure2_pu)*Tx2_rr_pu + (r_Tx2*p_cure2_pr_en)*Tx2_rr_pr_en + (r_Tx2*p_cure2_pr_nen)*Tx2_rr_pr_nen - 
    (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm)*Rlo_rr - 
    (r_mort)*Rlo_rr
  
  dRhi_ds <- r_self_cure*(I_ds) + r_ltfu1*(Tx_ds_pu + Tx_ds_pr_en +Tx_ds_pr_nen) - 
    (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm)*Rhi_ds -
    (r_mort)*Rhi_ds
  
  dRhi_rr <- r_self_cure*(I_rr + Tx_rr_pu + Tx_rr_pr_en + Tx_rr_pr_nen) + r_ltfu2_pu*Tx2_rr_pu + r_ltfu2_pr_en*Tx2_rr_pr_en + r_ltfu2_pr_nen*Tx2_rr_pr_nen - 
    (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm)*Rhi_rr - 
    (r_mort)*Rhi_rr
  
  dR_ds <- 0.5*(Rlo_ds + Rhi_ds) - 
    (r_relapse_R + (lam_ds + lam_rr)*p_imm)*R_ds - 
    (r_mort)*R_ds
  
  dR_rr <- 0.5*(Rlo_rr + Rhi_rr) - 
    (r_relapse_R + (lam_ds + lam_rr)*p_imm)*R_rr - 
    (r_mort)*R_rr
  
  dTx_ds_pu <- r_Tx_init*p_pu*I_ds - 
    (r_Tx + r_ltfu1 + r_MDR_acqu)*Tx_ds_pu - 
    (r_mort + r_muTx1)*Tx_ds_pu
  
  dTx_ds_pr_en <- r_Tx_init*(1-p_pu)*p_pr_en*I_ds - 
    (r_Tx + r_ltfu1 + r_MDR_acqu)*Tx_ds_pr_en - 
    (r_mort + r_muTx1)*Tx_ds_pr_en
  
  dTx_ds_pr_nen <- r_Tx_init*(1-p_pu)*(1-p_pr_en)*I_ds - 
    (r_Tx + r_ltfu1 + r_MDR_acqu)*Tx_ds_pr_nen - 
    (r_mort + r_muTx1)*Tx_ds_pr_nen
  
  dTx_rr_pu <- r_Tx_init*p_pu*(1-p_MDR_rec_pu)*I_rr + r_MDR_acqu*Tx_ds_pu - 
    (r_Tx + r_ltfu1 + r_self_cure)*Tx_rr_pu - 
    (r_mort + r_mort_TB)*Tx_rr_pu
  
  dTx_rr_pr_en <- r_Tx_init*(1-p_pu)*p_pr_en*(1-p_MDR_rec_pr_en)*I_rr + r_MDR_acqu*Tx_ds_pr_en - 
    (r_Tx + r_ltfu1 + r_self_cure)*Tx_rr_pr_en - 
    (r_mort + r_mort_TB)*Tx_rr_pr_en
  
  dTx_rr_pr_nen <- r_Tx_init*(1-p_pu)*(1-p_pr_en)*(1-p_MDR_rec_pr_nen)*I_rr + r_MDR_acqu*Tx_ds_pr_nen -
    (r_Tx + r_ltfu1 + r_self_cure)*Tx_rr_pr_nen - 
    (r_mort + r_mort_TB)*Tx_rr_pr_nen
  
  dTx2_ds_pu <- 0
  
  dTx2_ds_pr_en <- 0
  
  dTx2_ds_pr_nen <- 0
  
  dTx2_rr_pu <- r_Tx_init*p_pu*p_MDR_rec_pu*I_rr + r_Tx*p_SL_trans_pu*Tx_rr_pu - 
    (r_Tx2 + r_ltfu2_pu)*Tx2_rr_pu - 
    (r_mort + r_muTx2_pu)*Tx2_rr_pu
  
  dTx2_rr_pr_en <- r_Tx_init*(1-p_pu)*p_pr_en*p_MDR_rec_pr_en*I_rr + r_Tx*p_SL_trans_pr_en*Tx_rr_pr_en - 
    (r_Tx2 + r_ltfu2_pr_en)*Tx2_rr_pr_en - 
    (r_mort + r_muTx2_pr_en)*Tx2_rr_pr_en
  
  dTx2_rr_pr_nen <- r_Tx_init*(1-p_pu)*(1-p_pr_en)*p_MDR_rec_pr_nen*I_rr + r_Tx*p_SL_trans_pr_nen*Tx_rr_pr_nen -
    (r_Tx2 + r_ltfu2_pr_nen)*Tx2_rr_pr_nen - 
    (r_mort + r_muTx2_pr_nen)*Tx2_rr_pr_nen
  
  # auxiliaries
  dcu_inc_ds <- (r_progression*Lf_ds + r_reactivation*Ls_ds + 
                   r_relapse_Rlo*Rlo_ds + r_relapse_Rhi*Rhi_ds + 
                   r_relapse_R*R_ds)
  
  dcu_inc_rr <- (r_progression*Lf_rr + r_reactivation*Ls_rr + 
                   r_relapse_Rlo*Rlo_rr + r_relapse_Rhi*Rhi_rr + 
                   r_relapse_R*R_rr)
  
  # (not counting people going from Tx_rr_pu to Tx2_rr_pu and in pr for now. So we are undercounting the rr notifications)
  noti_ds_pu <- r_Tx_init*p_pu*I_ds + r_Tx_init*p_pu*(1-p_MDR_rec_pu)*I_rr
  noti_rr_pu <- r_Tx_init*p_pu*p_MDR_rec_pu*I_rr
  noti_ds_pr_en <- r_Tx_init*(1-p_pu)*p_pr_en*I_ds + r_Tx_init*(1-p_pu)*p_pr_en*(1-p_MDR_rec_pr_en)*I_rr
  noti_rr_pr_en <- r_Tx_init*(1-p_pu)*p_pr_en*(p_MDR_rec_pr_en)*I_rr
  
  dcu_mort_tb <- r_mort_TB*(I_ds + I_rr + Tx_rr_pu + Tx_rr_pr_en + Tx_rr_pr_nen) + 
    r_muTx1*(Tx_ds_pu + Tx_ds_pr_en + Tx_ds_pr_nen) + 
    r_muTx2_pu*Tx2_rr_pu + r_muTx2_pr_en*Tx2_rr_pr_en + r_muTx2_pr_nen*Tx2_rr_pr_nen
  
  dcu_mort_tb_ds <- r_mort_TB*(I_ds) + 
    r_muTx1*(Tx_ds_pu + Tx_ds_pr_en + Tx_ds_pr_nen)
  
  dcu_mort_tb_rr <- r_mort_TB*(I_rr + Tx_rr_pu + Tx_rr_pr_en + Tx_rr_pr_nen) + 
    r_muTx2_pu*Tx2_rr_pu + r_muTx2_pr_en*Tx2_rr_pr_en + r_muTx2_pr_nen*Tx2_rr_pr_nen
  
  dcu_noti_rr_pu <- noti_rr_pu
  dcu_noti_rr_pr_en <- noti_rr_pr_en
  dcu_noti_pu <- noti_ds_pu + noti_rr_pu
  dcu_noti_pr_en <- noti_ds_pr_en + noti_rr_pr_en
  
  return(list(c(dU, dLf_ds, dLf_rr, dLs_ds, dLs_rr, dI_ds, dI_rr,
                dRlo_ds, dRlo_rr, dRhi_ds, dRhi_rr, dR_ds, dR_rr,
                dTx_ds_pu, dTx_ds_pr_en, dTx_ds_pr_nen, dTx_rr_pu, dTx_rr_pr_en, dTx_rr_pr_nen, 
                dTx2_ds_pu, dTx2_ds_pr_en, dTx2_ds_pr_nen, dTx2_rr_pu, dTx2_rr_pr_en, dTx2_rr_pr_nen,
                dcu_inc_ds, dcu_inc_rr, dcu_noti_pu, dcu_mort_tb, dcu_mort_tb_ds, dcu_mort_tb_rr, 
                dcu_noti_pr_en, dcu_noti_rr_pu, dcu_noti_rr_pr_en)))
}

goveqs_m2_1 <- function(t, state, parameters) {
  r <- parameters[["r"]]
  p <- parameters[["p"]]
  
  # rates and proportions
  r_beta_ds              <- r[["beta"]][1]
  r_beta_rr              <- r[["beta"]][2]
  r_Tx_init              <- r[["Tx_init"]]
  r_progression          <- r[["progression"]]
  r_stabilisation        <- r[["stabilisation"]]
  r_reactivation         <- r[["reactivation"]]
  r_relapse_Rlo          <- r[["relapse"]][1]
  r_relapse_Rhi          <- r[["relapse"]][2]
  r_relapse_R            <- r[["relapse"]][3]
  r_self_cure            <- r[["self_cure"]]
  r_MDR_acqu             <- r[["MDR_acqu"]]
  r_Tx                   <- r[["Tx"]]
  r_Tx2                  <- r[["Tx2"]]
  r_ltfu1                <- r[["ltfu"]][1]
  r_ltfu2_pu             <- r[["ltfu"]][2]
  r_ltfu2_pr_nen         <- r[["ltfu"]][3]
  
  p_pu                   <- p[["pu"]]
  p_imm                  <- p[["imm"]]
  p_MDR_rec_pu           <- p[["MDR_rec"]][1]
  p_MDR_rec_pr_nen       <- p[["MDR_rec_pr_nen"]]
  p_SL_trans_pu          <- p[["SL_trans"]]
  p_SL_trans_pr_nen      <- p[["SL_trans_pr_nen"]]
  p_cure1                <- p[["cure"]][1]
  p_cure2_pu             <- p[["cure"]][2]
  p_cure2_pr_nen         <- p[["cure"]][3]
  
  r_mort                 <- r[["mort"]]
  r_mort_TB              <- r[["mort_TB"]]
  r_muTx1                <- r[["muTx"]][1]
  r_muTx2_pu             <- r[["muTx"]][2]
  r_muTx2_pr_nen         <- r[["muTx"]][3]
  
  # states
  U                      <- state[1]
  Lf_ds                  <- state[2]
  Lf_rr                  <- state[3]
  Ls_ds                  <- state[4]
  Ls_rr                  <- state[5]
  I_ds                   <- state[6]
  I_rr                   <- state[7]
  Rlo_ds                 <- state[8]
  Rlo_rr                 <- state[9]
  Rhi_ds                 <- state[10]
  Rhi_rr                 <- state[11]
  R_ds                   <- state[12]
  R_rr                   <- state[13]
  Tx_ds_pu               <- state[14]
  Tx_ds_pr_nen           <- state[15]
  Tx_rr_pu               <- state[16]
  Tx_rr_pr_nen           <- state[17]
  Tx2_ds_pu              <- state[18]
  Tx2_ds_pr_nen          <- state[19]
  Tx2_rr_pu              <- state[20]
  Tx2_rr_pr_nen          <- state[21]
  
  # auxiliaries
  cu_inc_ds              <- state[22]
  cu_inc_rr              <- state[23]
  cu_noti_pu             <- state[24]
  cu_mort_tb             <- state[25]
  cu_mort_tb_ds          <- state[26]
  cu_mort_tb_rr          <- state[27]
  cu_noti_rr_pu          <- state[28]
  
  # compute some useful values
  N <- sum(state[1:21])
  lam_ds <- r_beta_ds*(I_ds)/N  # ds FOI
  lam_rr <- r_beta_rr*(I_rr + Tx_rr_pu + Tx_rr_pr_nen)/N  # rr FOI
  births <- (r_mort)*(U + Lf_ds + Lf_rr + Ls_ds + Ls_rr + 
                        Rlo_ds + Rlo_rr + Rhi_ds + Rhi_rr + 
                        R_ds + R_rr) + 
    (r_mort + r_mort_TB)*(I_ds + I_rr + 
                            Tx_rr_pu + Tx_rr_pr_nen) +
    (r_mort + r_muTx1)*(Tx_ds_pu + Tx_ds_pr_nen) +
    (r_mort + r_muTx2_pu)*Tx2_rr_pu + (r_mort + r_muTx2_pr_nen)*Tx2_rr_pr_nen
  
  # rates of change
  dU <- births - (lam_ds + lam_rr)*U - (r_mort)*U
  
  dLf_ds <- lam_ds*U + lam_ds*p_imm*(Ls_ds + Ls_rr + Rlo_ds + Rlo_rr + Rhi_ds + Rhi_rr + R_ds + R_rr) - 
    (r_stabilisation + r_progression)*Lf_ds - 
    (r_mort)*Lf_ds
  
  dLf_rr <- lam_rr*U + lam_rr*p_imm*(Ls_ds + Ls_rr + Rlo_ds + Rlo_rr + Rhi_ds + Rhi_rr + R_ds + R_rr) - 
    (r_stabilisation + r_progression)*Lf_rr - 
    (r_mort)*Lf_rr
  
  dLs_ds <- r_stabilisation*Lf_ds - 
    (r_reactivation + (lam_ds + lam_rr)*p_imm)*Ls_ds - 
    (r_mort)*Ls_ds
  
  dLs_rr <- r_stabilisation*Lf_rr - 
    (r_reactivation + (lam_ds + lam_rr)*p_imm)*Ls_rr - 
    (r_mort)*Ls_rr
  
  dI_ds <- (r_progression*Lf_ds + r_reactivation*Ls_ds + r_relapse_Rlo*Rlo_ds + r_relapse_Rhi*Rhi_ds + r_relapse_R*R_ds) + 
    r_Tx*(1-p_cure1)*(Tx_ds_pu + Tx_ds_pr_nen) - 
    (r_Tx_init + r_self_cure)*I_ds - 
    (r_mort + r_mort_TB)*I_ds
  
  dI_rr <- (r_progression*Lf_rr + r_reactivation*Ls_rr + r_relapse_Rlo*Rlo_rr + r_relapse_Rhi*Rhi_rr + r_relapse_R*R_rr) + 
    (r_Tx*(1-p_SL_trans_pu) + r_ltfu1)*Tx_rr_pu + (r_Tx*(1-p_SL_trans_pr_nen) + r_ltfu1)*Tx_rr_pr_nen + 
    (r_Tx2*(1-p_cure2_pu))*Tx2_rr_pu + (r_Tx2*(1-p_cure2_pr_nen))*Tx2_rr_pr_nen - 
    (r_Tx_init + r_self_cure)*I_rr - 
    (r_mort + r_mort_TB)*I_rr
  
  dRlo_ds <- (r_Tx*p_cure1)*(Tx_ds_pu + Tx_ds_pr_nen) - 
    (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm)*Rlo_ds - 
    (r_mort)*Rlo_ds
  
  dRlo_rr <- (r_Tx2*p_cure2_pu)*Tx2_rr_pu + (r_Tx2*p_cure2_pr_nen)*Tx2_rr_pr_nen - 
    (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm)*Rlo_rr - 
    (r_mort)*Rlo_rr
  
  dRhi_ds <- r_self_cure*(I_ds) + r_ltfu1*(Tx_ds_pu +Tx_ds_pr_nen) - 
    (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm)*Rhi_ds -
    (r_mort)*Rhi_ds
  
  dRhi_rr <- r_self_cure*(I_rr + Tx_rr_pu + Tx_rr_pr_nen) + r_ltfu2_pu*Tx2_rr_pu + r_ltfu2_pr_nen*Tx2_rr_pr_nen - 
    (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm)*Rhi_rr - 
    (r_mort)*Rhi_rr
  
  dR_ds <- 0.5*(Rlo_ds + Rhi_ds) - 
    (r_relapse_R + (lam_ds + lam_rr)*p_imm)*R_ds - 
    (r_mort)*R_ds
  
  dR_rr <- 0.5*(Rlo_rr + Rhi_rr) - 
    (r_relapse_R + (lam_ds + lam_rr)*p_imm)*R_rr - 
    (r_mort)*R_rr
  
  dTx_ds_pu <- r_Tx_init*p_pu*I_ds - 
    (r_Tx + r_ltfu1 + r_MDR_acqu)*Tx_ds_pu - 
    (r_mort + r_muTx1)*Tx_ds_pu
  
  dTx_ds_pr_nen <- r_Tx_init*(1-p_pu)*I_ds - 
    (r_Tx + r_ltfu1 + r_MDR_acqu)*Tx_ds_pr_nen - 
    (r_mort + r_muTx1)*Tx_ds_pr_nen
  
  dTx_rr_pu <- r_Tx_init*p_pu*(1-p_MDR_rec_pu)*I_rr + r_MDR_acqu*Tx_ds_pu - 
    (r_Tx + r_ltfu1 + r_self_cure)*Tx_rr_pu - 
    (r_mort + r_mort_TB)*Tx_rr_pu
  
  dTx_rr_pr_nen <- r_Tx_init*(1-p_pu)*(1-p_MDR_rec_pr_nen)*I_rr + r_MDR_acqu*Tx_ds_pr_nen -
    (r_Tx + r_ltfu1 + r_self_cure)*Tx_rr_pr_nen - 
    (r_mort + r_mort_TB)*Tx_rr_pr_nen
  
  dTx2_ds_pu <- 0
  
  dTx2_ds_pr_nen <- 0
  
  dTx2_rr_pu <- r_Tx_init*p_pu*p_MDR_rec_pu*I_rr + r_Tx*p_SL_trans_pu*Tx_rr_pu - 
    (r_Tx2 + r_ltfu2_pu)*Tx2_rr_pu - 
    (r_mort + r_muTx2_pu)*Tx2_rr_pu
  
  dTx2_rr_pr_nen <- r_Tx_init*(1-p_pu)*p_MDR_rec_pr_nen*I_rr + r_Tx*p_SL_trans_pr_nen*Tx_rr_pr_nen -
    (r_Tx2 + r_ltfu2_pr_nen)*Tx2_rr_pr_nen - 
    (r_mort + r_muTx2_pr_nen)*Tx2_rr_pr_nen
  
  # auxiliaries
  dcu_inc_ds <- (r_progression*Lf_ds + r_reactivation*Ls_ds + 
                   r_relapse_Rlo*Rlo_ds + r_relapse_Rhi*Rhi_ds + 
                   r_relapse_R*R_ds)
  
  dcu_inc_rr <- (r_progression*Lf_rr + r_reactivation*Ls_rr + 
                   r_relapse_Rlo*Rlo_rr + r_relapse_Rhi*Rhi_rr + 
                   r_relapse_R*R_rr)
  
  # (not counting people going from Tx_rr_pu to Tx2_rr_pu and in pr for now. So we are undercounting the rr notifications)
  noti_ds_pu <- r_Tx_init*p_pu*I_ds + r_Tx_init*p_pu*(1-p_MDR_rec_pu)*I_rr
  noti_rr_pu <- r_Tx_init*p_pu*p_MDR_rec_pu*I_rr
  
  dcu_mort_tb <- r_mort_TB*(I_ds + I_rr + Tx_rr_pu + Tx_rr_pr_nen) + 
    r_muTx1*(Tx_ds_pu + Tx_ds_pr_nen) + 
    r_muTx2_pu*Tx2_rr_pu + r_muTx2_pr_nen*Tx2_rr_pr_nen
  
  dcu_mort_tb_ds <- r_mort_TB*(I_ds) + 
    r_muTx1*(Tx_ds_pu + Tx_ds_pr_nen)
  
  dcu_mort_tb_rr <- r_mort_TB*(I_rr + Tx_rr_pu + Tx_rr_pr_nen) + 
    r_muTx2_pu*Tx2_rr_pu + r_muTx2_pr_nen*Tx2_rr_pr_nen
  
  dcu_noti_rr_pu <- noti_rr_pu
  dcu_noti_pu <- noti_ds_pu + noti_rr_pu
  
  return(list(c(dU, dLf_ds, dLf_rr, dLs_ds, dLs_rr, dI_ds, dI_rr,
                dRlo_ds, dRlo_rr, dRhi_ds, dRhi_rr, dR_ds, dR_rr,
                dTx_ds_pu, dTx_ds_pr_nen, dTx_rr_pu, dTx_rr_pr_nen, 
                dTx2_ds_pu, dTx2_ds_pr_nen, dTx2_rr_pu, dTx2_rr_pr_nen,
                dcu_inc_ds, dcu_inc_rr, dcu_noti_pu, dcu_mort_tb, dcu_mort_tb_ds, dcu_mort_tb_rr, 
                dcu_noti_rr_pu)))
}

goveqs_m3_1 <- function(t, state, parameters) {
  r               <- parameters[["r"]]
  p               <- parameters[["p"]]
  
  # Rates and proportions
  r_beta_ds       <- r[["beta"]][1]
  r_beta_rr       <- r[["beta"]][2]
  r_Tx_init       <- r[["Tx_init"]]
  r_progression   <- r[["progression"]]
  r_stabilisation <- r[["stabilisation"]]
  r_reactivation  <- r[["reactivation"]]
  r_relapse_Rlo   <- r[["relapse"]][1]
  r_relapse_Rhi   <- r[["relapse"]][2]
  r_relapse_R     <- r[["relapse"]][3]
  r_self_cure     <- r[["self_cure"]]
  r_MDR_acqu      <- r[["MDR_acqu"]]
  r_Tx            <- r[["Tx"]]
  r_Tx2           <- r[["Tx2"]]
  r_ltfu1         <- r[["ltfu"]][1]
  r_ltfu2         <- r[["ltfu"]][2]
  
  p_imm           <- p[["imm"]]
  p_MDR_rec       <- p[["MDR_rec"]]
  p_SL_trans      <- p[["SL_trans"]]
  p_cure1         <- p[["cure"]][1]
  p_cure2         <- p[["cure"]][2]
  p_mort_Tx_TB    <- p[["mort_Tx_TB"]]
  
  r_mort          <- r[["mort"]]
  r_mort_TB       <- r[["mort_TB"]]
  r_muTx1         <- r[["muTx"]][1]
  r_muTx2         <- r[["muTx"]][2]
  
  # states
  U               <- state[1]
  Lf_ds           <- state[2]
  Lf_rr           <- state[3]
  Ls_ds           <- state[4]
  Ls_rr           <- state[5]
  I_ds            <- state[6]
  I_rr            <- state[7]
  Rlo_ds          <- state[8]
  Rlo_rr          <- state[9]
  Rhi_ds          <- state[10]
  Rhi_rr          <- state[11]
  R_ds            <- state[12]
  R_rr            <- state[13]
  Tx_ds           <- state[14]
  Tx_rr           <- state[15]
  Tx2_ds          <- state[16]
  Tx2_rr          <- state[17]
  
  # auxiliaries
  cu_inc_ds       <- state[18]
  cu_inc_rr       <- state[19]
  cu_noti         <- state[20]
  cu_mort_tb      <- state[21]
  cu_mort_tb_ds   <- state[22]
  cu_mort_tb_rr   <- state[23]
  cu_noti_rr      <- state[24]
  
  # compute some useful values
  N <- sum(state[1:17])
  lam_ds <- r_beta_ds*(I_ds)/N  # ds FOI
  lam_rr <- r_beta_rr*(I_rr + Tx_rr)/N  # rr FOI
  births <- (r_mort)*(U + Lf_ds + Lf_rr + Ls_ds + Ls_rr + 
                        Rlo_ds + Rlo_rr + Rhi_ds + Rhi_rr + R_ds + R_rr) + 
    (r_mort + r_mort_TB)*(I_ds + I_rr + Tx_rr) +
    (r_mort + r_muTx1)*(Tx_ds) +
    (r_mort + r_muTx2)*(Tx2_rr)
  
  # Rates of change
  dU <- births - 
    (lam_ds + lam_rr)*U - (r_mort)*U
  
  dLf_ds <- lam_ds*U + lam_ds*p_imm*(Ls_ds + Ls_rr + Rlo_ds + Rlo_rr + Rhi_ds + Rhi_rr + R_ds + R_rr) - 
    (r_stabilisation + r_progression)*Lf_ds - 
    (r_mort)*Lf_ds
  
  dLf_rr <- lam_rr*U + lam_rr*p_imm*(Ls_ds + Ls_rr + Rlo_ds + Rlo_rr + Rhi_ds + Rhi_rr + R_ds + R_rr) - 
    (r_stabilisation + r_progression)*Lf_rr - 
    (r_mort)*Lf_rr
  
  dLs_ds <- r_stabilisation*Lf_ds - 
    (r_reactivation + (lam_ds + lam_rr)*p_imm)*Ls_ds - 
    (r_mort)*Ls_ds
  
  dLs_rr <- r_stabilisation*Lf_rr - 
    (r_reactivation + (lam_ds + lam_rr)*p_imm)*Ls_rr - 
    (r_mort)*Ls_rr
  
  dI_ds <- (r_progression*Lf_ds + r_reactivation*Ls_ds + r_relapse_Rlo*Rlo_ds + r_relapse_Rhi*Rhi_ds + r_relapse_R*R_ds) +
    r_Tx*(1-p_cure1)*(Tx_ds) -
    (r_Tx_init + r_self_cure)*I_ds - 
    (r_mort + r_mort_TB)*I_ds
  
  dI_rr <- (r_progression*Lf_rr + r_reactivation*Ls_rr + r_relapse_Rlo*Rlo_rr + r_relapse_Rhi*Rhi_rr + r_relapse_R*R_rr) +
    (r_Tx*(1-p_SL_trans) + r_ltfu1)*Tx_rr +
    (r_Tx2*(1-p_cure2))*(Tx2_rr) -
    (r_Tx_init + r_self_cure)*I_rr - 
    (r_mort + r_mort_TB)*I_rr
  
  dRlo_ds <- (r_Tx*p_cure1)*(Tx_ds) - 
    (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm)*Rlo_ds - 
    (r_mort)*Rlo_ds
  
  dRlo_rr <- (r_Tx2*p_cure2)*(Tx2_rr) - 
    (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm)*Rlo_rr - 
    (r_mort)*Rlo_rr
  
  dRhi_ds <- r_self_cure*(I_ds) + 
    r_ltfu1*(Tx_ds) - 
    (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm)*Rhi_ds - 
    (r_mort)*Rhi_ds
  
  dRhi_rr <- r_self_cure*(I_rr + Tx_rr) +
    r_ltfu2*(Tx2_rr)- 
    (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm)*Rhi_rr - 
    (r_mort)*Rhi_rr
  
  dR_ds <- 0.5*(Rlo_ds + Rhi_ds) - 
    (r_relapse_R + (lam_ds + lam_rr)*p_imm)*R_ds - 
    (r_mort)*R_ds
  
  dR_rr <- 0.5*(Rlo_rr + Rhi_rr) - 
    (r_relapse_R + (lam_ds + lam_rr)*p_imm)*R_rr - 
    (r_mort)*R_rr
  
  dTx_ds <- r_Tx_init*I_ds - 
    (r_Tx + r_ltfu1 + r_MDR_acqu)*Tx_ds - 
    (r_mort + r_muTx1)*Tx_ds
  
  dTx_rr <- r_Tx_init*(1-p_MDR_rec)*I_rr + r_MDR_acqu*Tx_ds  - 
    (r_Tx + r_ltfu1 + r_self_cure)*Tx_rr - 
    (r_mort + r_mort_TB)*Tx_rr
  
  dTx2_ds <- 0
  
  dTx2_rr <- r_Tx_init*p_MDR_rec*I_rr + r_Tx*p_SL_trans*Tx_rr - 
    (r_Tx2 + r_ltfu2)*Tx2_rr - 
    (r_mort + r_muTx2)*Tx2_rr
  
  # Auxiliaries
  dcu_inc_ds <- (r_progression*Lf_ds + r_reactivation*Ls_ds + 
                   r_relapse_Rlo*Rlo_ds + r_relapse_Rhi*Rhi_ds + 
                   r_relapse_R*R_ds)
  
  dcu_inc_rr <- (r_progression*Lf_rr + r_reactivation*Ls_rr + 
                   r_relapse_Rlo*Rlo_rr + r_relapse_Rhi*Rhi_rr + 
                   r_relapse_R*R_rr)
  
  # not counting people going from Tx_rr_pu to Tx2_rr_pu and in pr for now. So we are undercounting the rr notifications
  dcu_noti_ds <- r_Tx_init*I_ds + r_Tx_init*(1-p_MDR_rec)*I_rr
  dcu_noti_rr <- r_Tx_init*p_MDR_rec*I_rr
  dcu_noti <- dcu_noti_ds + dcu_noti_rr
  
  dcu_mort_tb <- r_mort_TB*(I_ds + I_rr + Tx_rr) +
    r_muTx1*p_mort_Tx_TB*(Tx_ds) + 
    r_muTx2*p_mort_Tx_TB*(Tx2_rr)
  
  dcu_mort_tb_ds <- r_mort_TB*(I_ds) +
    r_muTx1*p_mort_Tx_TB*(Tx_ds)
  
  dcu_mort_tb_rr <- r_mort_TB*(I_rr + Tx_rr) +
    r_muTx2*p_mort_Tx_TB*(Tx2_rr)
  
  return(list(c(dU, dLf_ds, dLf_rr, dLs_ds, dLs_rr, dI_ds, dI_rr,
                dRlo_ds, dRlo_rr, dRhi_ds, dRhi_rr, dR_ds, dR_rr,
                dTx_ds, dTx_rr, 
                dTx2_ds, dTx2_rr,
                dcu_inc_ds, dcu_inc_rr, dcu_noti, dcu_mort_tb, dcu_mort_tb_ds, dcu_mort_tb_rr, 
                dcu_noti_rr)))
}


goveqs_m4_1 <- function(t, state, parameters) {
  r               <- parameters[["r"]]
  p               <- parameters[["p"]]
  i               <- parameters[["i"]]
  
  # Extract h0, h1, hart blocks using i
  mk_blk <- function(k) {
    list(
      U      = state[i$U[[k]]],
      Lf_ds  = state[i$Lf$ds[[k]]],  Lf_rr  = state[i$Lf$rr[[k]]],
      Ls_ds  = state[i$Ls$ds[[k]]],  Ls_rr  = state[i$Ls$rr[[k]]],
      I_ds   = state[i$I$ds[[k]]],   I_rr   = state[i$I$rr[[k]]],
      Rlo_ds = state[i$Rlo$ds[[k]]], Rlo_rr = state[i$Rlo$rr[[k]]],
      Rhi_ds = state[i$Rhi$ds[[k]]], Rhi_rr = state[i$Rhi$rr[[k]]],
      R_ds   = state[i$R$ds[[k]]],   R_rr   = state[i$R$rr[[k]]],
      Tx_ds  = state[i$Tx$ds[[k]]],  Tx_rr  = state[i$Tx$rr[[k]]],
      Tx2_ds = state[i$Tx2$ds[[k]]], Tx2_rr = state[i$Tx2$rr[[k]]]
    )
  }
  
  h0   <- mk_blk("h0")
  h1   <- mk_blk("h1")
  hart <- mk_blk("hart")
  
  cu_inc_ds     <- state[i$aux$inc[1]]
  cu_inc_rr     <- state[i$aux$inc[2]]
  cu_noti       <- state[i$aux$noti]
  cu_mort_tb    <- state[i$aux$mort_tb]
  cu_mort_tb_ds <- state[i$aux$mort_tb_ds]
  cu_mort_tb_rr <- state[i$aux$mort_tb_rr]
  cu_noti_rr    <- state[i$aux$noti_rr]
  cu_inc_h1hart <- state[i$aux$inc_h1hart]
  cu_mort_h0    <- state[i$aux$mort_tb_h0]
  cu_mort_h1hart<- state[i$aux$mort_tb_h1hart]
  
  
  # Rates and proportions
  r_beta_ds       <- r[["beta"]][1]
  r_beta_rr       <- r[["beta"]][2]
  r_Tx_init       <- r[["Tx_init"]]
  r_MDR_acqu      <- r[["MDR_acqu"]]
  r_relapse_Rlo   <- r[["relapse"]][1]
  r_relapse_Rhi   <- r[["relapse"]][2]
  r_relapse_R     <- r[["relapse"]][3]
  r_Tx            <- r[["Tx"]]
  r_Tx2           <- r[["Tx2"]]
  r_ltfu1         <- r[["ltfu"]][1]
  r_ltfu2         <- r[["ltfu"]][2]
  r_ltfu1_h1hart  <- r[["ltfu"]][3]
  r_ltfu2_h1hart  <- r[["ltfu"]][2]
  
  
  p_MDR_rec       <- p[["MDR_rec"]]
  p_SL_trans      <- p[["SL_trans"]]
  p_cure1         <- p[["cure"]][1]
  p_cure2         <- p[["cure"]][2]
  p_cure1_h1hart  <- p[["cure"]][3]
  p_cure2_h1hart  <- p[["cure"]][2]
  
  r_muTx1         <- r[["muTx"]][1]
  r_muTx2         <- r[["muTx"]][2]
  r_muTx1_h1hart  <- r[["muTx"]][3]
  r_muTx2_h1hart  <- r[["muTx"]][2]
  
  r_progression   <- r[["progression"]]     # h0, h1, hart
  r_stabilisation <- r[["stabilisation"]]   # h0, h1, hart
  r_reactivation  <- r[["reactivation"]]    # h0, h1, hart
  r_self_cure     <- r[["self_cure"]]       # h0, h1, hart
  p_imm           <- p[["imm"]]             # h0, h1, hart, with h1 = 0
  r_ART_init      <- r[["ART_init"]]
  r_mort_h0       <- r[["mort"]]
  r_mort_h1       <- r[["HIV_mort"]]
  r_mort_hart     <- r[["mort"]]
  r_mort_TB_h0    <- r[["mort_TB"]][1]
  r_mort_TB_h1    <- r[["mort_TB"]][2]
  r_mort_TB_hart  <- r[["mort_TB"]][1]
  
  # externally supplied annual HIV incidence rate used to move people from HIV negative to HIV positive
  r_HIV <- if (is.function(p[["rHIV"]])) p[["rHIV"]](t) else p[["rHIV"]]
  
  # compute some useful values
  N <- sum(state[1:i$nstates])
  lam_ds <- r_beta_ds*sum(c(h0$I_ds, h1$I_ds, hart$I_ds))/N # ds FOI
  lam_rr <- r_beta_rr*sum(c(h0$I_rr, h1$I_rr, hart$I_rr,
                            h0$Tx_rr, h1$Tx_rr, hart$Tx_rr))/N # rr FOI
  
  # natural mortality
  mort_part_nat <- (r_mort_h0)*(h0$U + h0$Lf_ds + h0$Lf_rr + h0$Ls_ds + h0$Ls_rr +
                                  h0$Rlo_ds + h0$Rlo_rr + h0$Rhi_ds + h0$Rhi_rr +
                                  h0$R_ds + h0$R_rr) +
    (r_mort_h1)*(h1$U + h1$Lf_ds + h1$Lf_rr + h1$Ls_ds + h1$Ls_rr +
                   h1$Rlo_ds + h1$Rlo_rr + h1$Rhi_ds + h1$Rhi_rr +
                   h1$R_ds + h1$R_rr) +
    (r_mort_hart)*(hart$U + hart$Lf_ds + hart$Lf_rr + hart$Ls_ds + hart$Ls_rr +
                     hart$Rlo_ds + hart$Rlo_rr + hart$Rhi_ds + hart$Rhi_rr +
                     hart$R_ds + hart$R_rr)
  
  # active TB mortality
  mort_part_TB <- (r_mort_h0 + r_mort_TB_h0)*(h0$I_ds + h0$I_rr + h0$Tx_rr) +
    (r_mort_h1 + r_mort_TB_h1)*(h1$I_ds + h1$I_rr + h1$Tx_rr) +
    (r_mort_hart + r_mort_TB_hart)*(hart$I_ds + hart$I_rr + hart$Tx_rr)
  
  # mortality during treatment
  mort_part_treat <- (r_mort_h0 + r_muTx1)*h0$Tx_ds +
    (r_mort_h1 + r_muTx1_h1hart)*h1$Tx_ds +
    (r_mort_hart + r_muTx1_h1hart)*hart$Tx_ds +
    (r_mort_h0 + r_muTx2)*h0$Tx2_rr +
    (r_mort_h1 + r_muTx2_h1hart)*h1$Tx2_rr +
    (r_mort_hart + r_muTx2_h1hart)*hart$Tx2_rr 
  
  births <-  mort_part_nat + mort_part_TB + mort_part_treat 
  
  # Compartments for reinfection (-> Lf) among partially immune states
  prev_h0   <- h0$Ls_ds + h0$Ls_rr + h0$Rlo_ds + h0$Rlo_rr + h0$Rhi_ds + h0$Rhi_rr + h0$R_ds + h0$R_rr
  prev_h1   <- h1$Ls_ds + h1$Ls_rr + h1$Rlo_ds + h1$Rlo_rr + h1$Rhi_ds + h1$Rhi_rr + h1$R_ds + h1$R_rr
  prev_hart <- hart$Ls_ds + hart$Ls_rr + hart$Rlo_ds + hart$Rlo_rr + hart$Rhi_ds + hart$Rhi_rr + hart$R_ds + hart$R_rr
  
  # Rates of change
  # TB-only derivatives by HIV stratum
  tb_U <- c(
    births - (lam_ds + lam_rr + r_mort_h0)*h0$U,
    - (lam_ds + lam_rr + r_mort_h1)*h1$U,
    - (lam_ds + lam_rr + r_mort_hart)*hart$U
  )
  
  tb_Lf_ds <- c(
    lam_ds*h0$U + lam_ds*p_imm[1]*prev_h0 - (r_progression[1] + r_stabilisation[1] + r_mort_h0)*h0$Lf_ds,
    lam_ds*h1$U + lam_ds*p_imm[2]*prev_h1 - (r_progression[2] + r_stabilisation[2] + r_mort_h1)*h1$Lf_ds,
    lam_ds*hart$U + lam_ds*p_imm[3]*prev_hart - (r_progression[3] + r_stabilisation[3] + r_mort_hart)*hart$Lf_ds
  )
  
  tb_Lf_rr <- c(
    lam_rr*h0$U + lam_rr*p_imm[1]*prev_h0 - (r_progression[1] + r_stabilisation[1] + r_mort_h0)*h0$Lf_rr,
    lam_rr*h1$U + lam_rr*p_imm[2]*prev_h1 - (r_progression[2] + r_stabilisation[2] + r_mort_h1)*h1$Lf_rr,
    lam_rr*hart$U + lam_rr*p_imm[3]*prev_hart - (r_progression[3] + r_stabilisation[3] + r_mort_hart)*hart$Lf_rr
  )
  
  tb_Ls_ds <- c(
    r_stabilisation[1]*h0$Lf_ds - (r_reactivation[1] + (lam_ds + lam_rr)*p_imm[1] + r_mort_h0)*h0$Ls_ds,
    r_stabilisation[2]*h1$Lf_ds - (r_reactivation[2] + (lam_ds + lam_rr)*p_imm[2] + r_mort_h1)*h1$Ls_ds,
    r_stabilisation[3]*hart$Lf_ds - (r_reactivation[3] + (lam_ds + lam_rr)*p_imm[3] + r_mort_hart)*hart$Ls_ds
  )
  
  tb_Ls_rr <- c(
    r_stabilisation[1]*h0$Lf_rr - (r_reactivation[1] + (lam_ds + lam_rr)*p_imm[1] + r_mort_h0)*h0$Ls_rr,
    r_stabilisation[2]*h1$Lf_rr - (r_reactivation[2] + (lam_ds + lam_rr)*p_imm[2] + r_mort_h1)*h1$Ls_rr,
    r_stabilisation[3]*hart$Lf_rr - (r_reactivation[3] + (lam_ds + lam_rr)*p_imm[3] + r_mort_hart)*hart$Ls_rr
  )
  
  tb_I_ds <- c(
    r_progression[1]*h0$Lf_ds + r_reactivation[1]*h0$Ls_ds +
      r_relapse_Rlo*h0$Rlo_ds + r_relapse_Rhi*h0$Rhi_ds + r_relapse_R*h0$R_ds +
      r_Tx*(1 - p_cure1)*h0$Tx_ds -
      (r_Tx_init + r_self_cure[1] + r_mort_h0 + r_mort_TB_h0)*h0$I_ds,
    
    r_progression[2]*h1$Lf_ds + r_reactivation[2]*h1$Ls_ds +
      r_relapse_Rlo*h1$Rlo_ds + r_relapse_Rhi*h1$Rhi_ds + r_relapse_R*h1$R_ds +
      r_Tx*(1 - p_cure1_h1hart)*h1$Tx_ds -
      (r_Tx_init + r_self_cure[2] + r_mort_h1 + r_mort_TB_h1)*h1$I_ds,
    
    r_progression[3]*hart$Lf_ds + r_reactivation[3]*hart$Ls_ds +
      r_relapse_Rlo*hart$Rlo_ds + r_relapse_Rhi*hart$Rhi_ds + r_relapse_R*hart$R_ds +
      r_Tx*(1 - p_cure1_h1hart)*hart$Tx_ds -
      (r_Tx_init + r_self_cure[3] + r_mort_hart + r_mort_TB_hart)*hart$I_ds
  )
  
  tb_I_rr <- c(
    r_progression[1]*h0$Lf_rr + r_reactivation[1]*h0$Ls_rr +
      r_relapse_Rlo*h0$Rlo_rr + r_relapse_Rhi*h0$Rhi_rr + r_relapse_R*h0$R_rr +
      (r_Tx*(1 - p_SL_trans) + r_ltfu1)*h0$Tx_rr +
      r_Tx2 * (1 - p_cure2)*h0$Tx2_rr -
      (r_Tx_init + r_self_cure[1] + r_mort_h0 + r_mort_TB_h0)*h0$I_rr,
    
    r_progression[2]*h1$Lf_rr + r_reactivation[2]*h1$Ls_rr +
      r_relapse_Rlo*h1$Rlo_rr + r_relapse_Rhi*h1$Rhi_rr + r_relapse_R*h1$R_rr +
      (r_Tx*(1 - p_SL_trans) + r_ltfu1_h1hart)*h1$Tx_rr +
      r_Tx2*(1 - p_cure2_h1hart)*h1$Tx2_rr -
      (r_Tx_init + r_self_cure[2] + r_mort_h1 + r_mort_TB_h1)*h1$I_rr,
    
    r_progression[3]*hart$Lf_rr + r_reactivation[3]*hart$Ls_rr +
      r_relapse_Rlo*hart$Rlo_rr + r_relapse_Rhi*hart$Rhi_rr + r_relapse_R*hart$R_rr +
      (r_Tx*(1 - p_SL_trans) + r_ltfu1_h1hart)*hart$Tx_rr +
      r_Tx2*(1 - p_cure2_h1hart)*hart$Tx2_rr -
      (r_Tx_init + r_self_cure[3] + r_mort_hart + r_mort_TB_hart)*hart$I_rr
  )
  
  tb_Rlo_ds <- c(
    r_Tx*p_cure1*h0$Tx_ds - (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm[1] + r_mort_h0)*h0$Rlo_ds,
    r_Tx*p_cure1_h1hart*h1$Tx_ds - (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm[2] + r_mort_h1)*h1$Rlo_ds,
    r_Tx*p_cure1_h1hart*hart$Tx_ds - (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm[3] + r_mort_hart)*hart$Rlo_ds
  )
  
  tb_Rlo_rr <- c(
    r_Tx2*p_cure2*h0$Tx2_rr - (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm[1] + r_mort_h0)*h0$Rlo_rr,
    r_Tx2*p_cure2_h1hart*h1$Tx2_rr - (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm[2] + r_mort_h1)*h1$Rlo_rr,
    r_Tx2*p_cure2_h1hart*hart$Tx2_rr - (0.5 + r_relapse_Rlo + (lam_ds + lam_rr)*p_imm[3] + r_mort_hart)*hart$Rlo_rr
  )
  
  tb_Rhi_ds <- c(
    r_self_cure[1]*h0$I_ds + r_ltfu1*h0$Tx_ds - (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm[1] + r_mort_h0)*h0$Rhi_ds,
    r_self_cure[2]*h1$I_ds + r_ltfu1_h1hart*h1$Tx_ds - (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm[2] + r_mort_h1)*h1$Rhi_ds,
    r_self_cure[3]*hart$I_ds + r_ltfu1_h1hart*hart$Tx_ds - (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm[3] + r_mort_hart)*hart$Rhi_ds
  )
  
  tb_Rhi_rr <- c(
    r_self_cure[1]*(h0$I_rr + h0$Tx_rr) + r_ltfu2*h0$Tx2_rr - (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm[1] + r_mort_h0)*h0$Rhi_rr,
    r_self_cure[2]*(h1$I_rr + h1$Tx_rr) + r_ltfu2_h1hart*h1$Tx2_rr - (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm[2] + r_mort_h1)*h1$Rhi_rr,
    r_self_cure[3]*(hart$I_rr + hart$Tx_rr) + r_ltfu2_h1hart*hart$Tx2_rr - (0.5 + r_relapse_Rhi + (lam_ds + lam_rr)*p_imm[3] + r_mort_hart)*hart$Rhi_rr
  )
  
  tb_R_ds <- c(
    0.5*(h0$Rlo_ds + h0$Rhi_ds) - (r_relapse_R + (lam_ds + lam_rr)*p_imm[1] + r_mort_h0)*h0$R_ds,
    0.5*(h1$Rlo_ds + h1$Rhi_ds) - (r_relapse_R + (lam_ds + lam_rr)*p_imm[2] + r_mort_h1)*h1$R_ds,
    0.5*(hart$Rlo_ds + hart$Rhi_ds) - (r_relapse_R + (lam_ds + lam_rr)*p_imm[3] + r_mort_hart)*hart$R_ds
  )
  
  tb_R_rr <- c(
    0.5*(h0$Rlo_rr + h0$Rhi_rr) - (r_relapse_R + (lam_ds + lam_rr)*p_imm[1] + r_mort_h0)*h0$R_rr,
    0.5*(h1$Rlo_rr + h1$Rhi_rr) - (r_relapse_R + (lam_ds + lam_rr)*p_imm[2] + r_mort_h1)*h1$R_rr,
    0.5*(hart$Rlo_rr + hart$Rhi_rr) - (r_relapse_R + (lam_ds + lam_rr)*p_imm[3] + r_mort_hart)*hart$R_rr
  )
  
  tb_Tx_ds <- c(
    r_Tx_init*h0$I_ds - (r_Tx + r_ltfu1 + r_MDR_acqu + r_mort_h0 + r_muTx1)*h0$Tx_ds,
    r_Tx_init*h1$I_ds - (r_Tx + r_ltfu1_h1hart + r_MDR_acqu + r_mort_h1 + r_muTx1_h1hart)*h1$Tx_ds,
    r_Tx_init*hart$I_ds - (r_Tx + r_ltfu1_h1hart + r_MDR_acqu + r_mort_hart + r_muTx1_h1hart)*hart$Tx_ds
  )
  
  tb_Tx_rr <- c(
    r_Tx_init*(1 - p_MDR_rec)*h0$I_rr + r_MDR_acqu*h0$Tx_ds -
      (r_Tx + r_ltfu1 + r_self_cure[1] + r_mort_h0 + r_mort_TB_h0)*h0$Tx_rr,
    r_Tx_init*(1 - p_MDR_rec)*h1$I_rr + r_MDR_acqu*h1$Tx_ds -
      (r_Tx + r_ltfu1_h1hart + r_self_cure[2] + r_mort_h1 + r_mort_TB_h1)*h1$Tx_rr,
    r_Tx_init*(1 - p_MDR_rec)*hart$I_rr + r_MDR_acqu*hart$Tx_ds -
      (r_Tx + r_ltfu1_h1hart + r_self_cure[3] + r_mort_hart + r_mort_TB_hart)*hart$Tx_rr
  )
  
  tb_Tx2_ds <- c(0, 0, 0)
  
  tb_Tx2_rr <- c(
    r_Tx_init*p_MDR_rec*h0$I_rr + r_Tx*p_SL_trans*h0$Tx_rr -
      (r_Tx2 + r_ltfu2 + r_mort_h0 + r_muTx2)*h0$Tx2_rr,
    r_Tx_init*p_MDR_rec*h1$I_rr + r_Tx*p_SL_trans*h1$Tx_rr -
      (r_Tx2 + r_ltfu2_h1hart + r_mort_h1 + r_muTx2_h1hart)*h1$Tx2_rr,
    r_Tx_init*p_MDR_rec*hart$I_rr + r_Tx*p_SL_trans*hart$Tx_rr -
      (r_Tx2 + r_ltfu2_h1hart + r_mort_hart + r_muTx2_h1hart)*hart$Tx2_rr
  )
  
  # HIV transition helper: h0 -> h1 -> hart
  apply_hiv <- function(tb, x) {
    c(tb[1] - r_HIV*x[1],
      tb[2] + r_HIV*x[1] - r_ART_init*x[2],
      tb[3] + r_ART_init*x[2])}
  
  dU       <- apply_hiv(tb_U,       c(h0$U,      h1$U,      hart$U))
  dLf_ds   <- apply_hiv(tb_Lf_ds,   c(h0$Lf_ds,  h1$Lf_ds,  hart$Lf_ds))
  dLf_rr   <- apply_hiv(tb_Lf_rr,   c(h0$Lf_rr,  h1$Lf_rr,  hart$Lf_rr))
  dLs_ds   <- apply_hiv(tb_Ls_ds,   c(h0$Ls_ds,  h1$Ls_ds,  hart$Ls_ds))
  dLs_rr   <- apply_hiv(tb_Ls_rr,   c(h0$Ls_rr,  h1$Ls_rr,  hart$Ls_rr))
  dI_ds    <- apply_hiv(tb_I_ds,    c(h0$I_ds,   h1$I_ds,   hart$I_ds))
  dI_rr    <- apply_hiv(tb_I_rr,    c(h0$I_rr,   h1$I_rr,   hart$I_rr))
  dRlo_ds  <- apply_hiv(tb_Rlo_ds,  c(h0$Rlo_ds, h1$Rlo_ds, hart$Rlo_ds))
  dRlo_rr  <- apply_hiv(tb_Rlo_rr,  c(h0$Rlo_rr, h1$Rlo_rr, hart$Rlo_rr))
  dRhi_ds  <- apply_hiv(tb_Rhi_ds,  c(h0$Rhi_ds, h1$Rhi_ds, hart$Rhi_ds))
  dRhi_rr  <- apply_hiv(tb_Rhi_rr,  c(h0$Rhi_rr, h1$Rhi_rr, hart$Rhi_rr))
  dR_ds    <- apply_hiv(tb_R_ds,    c(h0$R_ds,   h1$R_ds,   hart$R_ds))
  dR_rr    <- apply_hiv(tb_R_rr,    c(h0$R_rr,   h1$R_rr,   hart$R_rr))
  dTx_ds   <- apply_hiv(tb_Tx_ds,   c(h0$Tx_ds,  h1$Tx_ds,  hart$Tx_ds))
  dTx_rr   <- apply_hiv(tb_Tx_rr,   c(h0$Tx_rr,  h1$Tx_rr,  hart$Tx_rr))
  dTx2_ds  <- apply_hiv(tb_Tx2_ds,  c(h0$Tx2_ds, h1$Tx2_ds, hart$Tx2_ds))
  dTx2_rr  <- apply_hiv(tb_Tx2_rr,  c(h0$Tx2_rr, h1$Tx2_rr, hart$Tx2_rr))
  
  # Auxiliaries
  dcu_inc_ds <- r_progression[1]*h0$Lf_ds + r_progression[2]*h1$Lf_ds + r_progression[3]*hart$Lf_ds +
    r_reactivation[1]*h0$Ls_ds + r_reactivation[2]*h1$Ls_ds + r_reactivation[3]*hart$Ls_ds +
    r_relapse_Rlo*(h0$Rlo_ds + h1$Rlo_ds + hart$Rlo_ds) +
    r_relapse_Rhi*(h0$Rhi_ds + h1$Rhi_ds + hart$Rhi_ds) +
    r_relapse_R*(h0$R_ds + h1$R_ds + hart$R_ds)
  
  dcu_inc_rr <- r_progression[1]*h0$Lf_rr + r_progression[2]*h1$Lf_rr + r_progression[3]*hart$Lf_rr +
    r_reactivation[1]*h0$Ls_rr + r_reactivation[2]*h1$Ls_rr + r_reactivation[3]*hart$Ls_rr +
    r_relapse_Rlo*(h0$Rlo_rr + h1$Rlo_rr + hart$Rlo_rr) +
    r_relapse_Rhi*(h0$Rhi_rr + h1$Rhi_rr + hart$Rhi_rr) +
    r_relapse_R*(h0$R_rr + h1$R_rr + hart$R_rr)
  
  
  # not counting people going from Tx_rr_pu to Tx2_rr_pu and in pr for now. So we are undercounting the rr notifications
  # SUMALI NOTE
  dcu_noti_ds <- r_Tx_init*(h0$I_ds + h1$I_ds + hart$I_ds) + r_Tx_init*(1 - p_MDR_rec)*(h0$I_rr + h1$I_rr + hart$I_rr)
  dcu_noti_rr <- r_Tx_init*p_MDR_rec*(h0$I_rr + h1$I_rr + hart$I_rr)
  dcu_noti <- dcu_noti_ds + dcu_noti_rr
  
  dcu_mort_tb <- r_mort_TB_h0*(h0$I_ds + h0$I_rr + h0$Tx_rr) +
    r_mort_TB_h1*(h1$I_ds + h1$I_rr + h1$Tx_rr) +
    r_mort_TB_hart*(hart$I_ds + hart$I_rr + hart$Tx_rr) +
    r_muTx1*(h0$Tx_ds) +
    r_muTx1_h1hart*(h1$Tx_ds + hart$Tx_ds) +
    r_muTx2*(h0$Tx2_rr ) +
    r_muTx2_h1hart*(h1$Tx2_rr + hart$Tx2_rr)
  
  dcu_mort_tb_ds <- r_mort_TB_h0*h0$I_ds + r_mort_TB_h1*h1$I_ds + r_mort_TB_hart*hart$I_ds +
    r_muTx1*(h0$Tx_ds) + r_muTx1_h1hart*(h1$Tx_ds + hart$Tx_ds)
  
  dcu_mort_tb_rr <- r_mort_TB_h0*(h0$I_rr + h0$Tx_rr) +
    r_mort_TB_h1*(h1$I_rr + h1$Tx_rr) + 
    r_mort_TB_hart*(hart$I_rr + hart$Tx_rr) +
    r_muTx2*(h0$Tx2_rr) +
    r_muTx2_h1hart*(h1$Tx2_rr + hart$Tx2_rr)

  dcu_inc_h1hart <- r_progression[2]*(h1$Lf_ds + h1$Lf_rr) +
    r_reactivation[2]*(h1$Ls_ds + h1$Ls_rr) +
    r_relapse_Rlo*(h1$Rlo_ds + h1$Rlo_rr) +
    r_relapse_Rhi*(h1$Rhi_ds + h1$Rhi_rr) +
    r_relapse_R *(h1$R_ds + h1$R_rr) +
    r_progression[3]*(hart$Lf_ds + hart$Lf_rr) +
    r_reactivation[3]*(hart$Ls_ds + hart$Ls_rr) +
    r_relapse_Rlo*(hart$Rlo_ds + hart$Rlo_rr) +
    r_relapse_Rhi*(hart$Rhi_ds + hart$Rhi_rr) +
    r_relapse_R *(hart$R_ds   + hart$R_rr)
  
  dcu_mort_tb_h0 <- r_mort_TB_h0*(h0$I_ds + h0$I_rr + h0$Tx_rr) +
    r_muTx1*h0$Tx_ds +
    r_muTx2*h0$Tx2_rr
  
  dcu_mort_tb_h1hart <- r_mort_TB_h1*(h1$I_ds + h1$I_rr + h1$Tx_rr) +
    r_mort_TB_hart*(hart$I_ds + hart$I_rr + hart$Tx_rr) +
    r_muTx1_h1hart*(h1$Tx_ds + hart$Tx_ds) +
    r_muTx2_h1hart*(h1$Tx2_rr + hart$Tx2_rr)
  
  # Pack output using i
  deriv <- rep(0, i$nx)
  
  deriv[i$U$h0]       <- dU[1];      deriv[i$U$h1]       <- dU[2];      deriv[i$U$hart]       <- dU[3]
  deriv[i$Lf$ds$h0]   <- dLf_ds[1];  deriv[i$Lf$ds$h1]   <- dLf_ds[2];  deriv[i$Lf$ds$hart]   <- dLf_ds[3]
  deriv[i$Lf$rr$h0]   <- dLf_rr[1];  deriv[i$Lf$rr$h1]   <- dLf_rr[2];  deriv[i$Lf$rr$hart]   <- dLf_rr[3]
  deriv[i$Ls$ds$h0]   <- dLs_ds[1];  deriv[i$Ls$ds$h1]   <- dLs_ds[2];  deriv[i$Ls$ds$hart]   <- dLs_ds[3]
  deriv[i$Ls$rr$h0]   <- dLs_rr[1];  deriv[i$Ls$rr$h1]   <- dLs_rr[2];  deriv[i$Ls$rr$hart]   <- dLs_rr[3]
  deriv[i$I$ds$h0]    <- dI_ds[1];   deriv[i$I$ds$h1]    <- dI_ds[2];   deriv[i$I$ds$hart]    <- dI_ds[3]
  deriv[i$I$rr$h0]    <- dI_rr[1];   deriv[i$I$rr$h1]    <- dI_rr[2];   deriv[i$I$rr$hart]    <- dI_rr[3]
  deriv[i$Rlo$ds$h0]  <- dRlo_ds[1]; deriv[i$Rlo$ds$h1]  <- dRlo_ds[2]; deriv[i$Rlo$ds$hart]  <- dRlo_ds[3]
  deriv[i$Rlo$rr$h0]  <- dRlo_rr[1]; deriv[i$Rlo$rr$h1]  <- dRlo_rr[2]; deriv[i$Rlo$rr$hart]  <- dRlo_rr[3]
  deriv[i$Rhi$ds$h0]  <- dRhi_ds[1]; deriv[i$Rhi$ds$h1]  <- dRhi_ds[2]; deriv[i$Rhi$ds$hart]  <- dRhi_ds[3]
  deriv[i$Rhi$rr$h0]  <- dRhi_rr[1]; deriv[i$Rhi$rr$h1]  <- dRhi_rr[2]; deriv[i$Rhi$rr$hart]  <- dRhi_rr[3]
  deriv[i$R$ds$h0]    <- dR_ds[1];   deriv[i$R$ds$h1]    <- dR_ds[2];   deriv[i$R$ds$hart]    <- dR_ds[3]
  deriv[i$R$rr$h0]    <- dR_rr[1];   deriv[i$R$rr$h1]    <- dR_rr[2];   deriv[i$R$rr$hart]    <- dR_rr[3]
  deriv[i$Tx$ds$h0]   <- dTx_ds[1];  deriv[i$Tx$ds$h1]   <- dTx_ds[2];  deriv[i$Tx$ds$hart]   <- dTx_ds[3]
  deriv[i$Tx$rr$h0]   <- dTx_rr[1];  deriv[i$Tx$rr$h1]   <- dTx_rr[2];  deriv[i$Tx$rr$hart]   <- dTx_rr[3]
  deriv[i$Tx2$ds$h0]  <- dTx2_ds[1]; deriv[i$Tx2$ds$h1]  <- dTx2_ds[2]; deriv[i$Tx2$ds$hart]  <- dTx2_ds[3]
  deriv[i$Tx2$rr$h0]  <- dTx2_rr[1]; deriv[i$Tx2$rr$h1]  <- dTx2_rr[2]; deriv[i$Tx2$rr$hart]  <- dTx2_rr[3]
  
  # Auxiliaries
  deriv[i$aux$inc[1]]       <- dcu_inc_ds
  deriv[i$aux$inc[2]]       <- dcu_inc_rr
  deriv[i$aux$noti]         <- dcu_noti
  deriv[i$aux$mort_tb]      <- dcu_mort_tb
  deriv[i$aux$mort_tb_ds]   <- dcu_mort_tb_ds
  deriv[i$aux$mort_tb_rr]   <- dcu_mort_tb_rr
  deriv[i$aux$noti_rr]      <- dcu_noti_rr
  deriv[i$aux$inc_h1hart]   <- dcu_inc_h1hart
  deriv[i$aux$mort_tb_h0]      <- dcu_mort_tb_h0
  deriv[i$aux$mort_tb_h1hart]  <- dcu_mort_tb_h1hart
  
  return(list(deriv))
}



goveqs_scaleup <- function(t, state, parameters) {
  # Extract start/end parameter sets and scale time
  goveqs <- parameters$goveqs
  p_start <- parameters$p_start
  p_end <- parameters$p_end
  r_start <- parameters$r_start
  r_end <- parameters$r_end
  times <- parameters$times
  
  # Calculate scale factor (0 to 1)
  scale <- min(max((t - times[1]) / (times[2] - times[1]), 0), 1)
  
  # Interpolate p
  p_t <- p_start
  for (name in names(p_start)) {
    if (is.numeric(p_start[[name]]) && length(p_start[[name]]) == length(p_end[[name]])) {
      p_t[[name]] <- p_start[[name]] + scale * (p_end[[name]] - p_start[[name]])
    }
  }
  
  # Interpolate r
  r_t <- r_start
  for (name in names(r_start)) {
    if (is.numeric(r_start[[name]]) && length(r_start[[name]]) == length(r_end[[name]])) {
      r_t[[name]] <- r_start[[name]] + scale * (r_end[[name]] - r_start[[name]])
    }
  }
  
  # Construct updated parameters list
  updated_parameters <- list(r = r_t, p = p_t)
  
  # Call main ODE function
  result <- goveqs(t, state, updated_parameters)
  return(result)
}


goveqs_scaleup_hiv <- function(t, state, parameters) {
  # Extract start/end parameter sets and scale time
  goveqs <- parameters$goveqs
  p_start <- parameters$p_start
  p_end <- parameters$p_end
  r_start <- parameters$r_start
  r_end <- parameters$r_end
  times <- parameters$times
  i <- parameters$i
  
  # Calculate scale factor (0 to 1)
  scale <- min(max((t - times[1]) / (times[2] - times[1]), 0), 1)
  
  # Interpolate p
  p_t <- p_start
  for (name in names(p_start)) {
    if (is.numeric(p_start[[name]]) && length(p_start[[name]]) == length(p_end[[name]])) {
      p_t[[name]] <- p_start[[name]] + scale * (p_end[[name]] - p_start[[name]])
    }
  }
  
  # Interpolate r
  r_t <- r_start
  for (name in names(r_start)) {
    if (is.numeric(r_start[[name]]) && length(r_start[[name]]) == length(r_end[[name]])) {
      r_t[[name]] <- r_start[[name]] + scale * (r_end[[name]] - r_start[[name]])
    }
  }
  
  # Construct updated parameters list
  updated_parameters <- list(r = r_t, p = p_t, i = i)
  
  # Call main ODE function
  result <- goveqs(t, state, updated_parameters)
  return(result)
}


