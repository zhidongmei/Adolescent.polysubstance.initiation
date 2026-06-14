
rm(list=ls());

proj_dir = "C:/Users/DZ1011/Dropbox/2024.PSU.lancet/"

library(tidyverse)

# ========================================================
# combine PGS, environment, PSU, covariates, together
# ========================================================
{
  # ==================
  # load PGS dataset
  # ==================
  Basic_information = read.csv(paste0(proj_dir, "data/psu/All_ABCD_samples_11868_basic_information.data__updated.25.01.19.csv")) 
  Basic_information$coreID <- gsub("NDAR_INV", "", Basic_information$IID)
  
  # ==================
  # load substance use
  # ==================
  load(paste0(proj_dir, "data/psu/abcd.6.0.substance_use.lifetime.NA.03.05.2026.Rdb"))  # substance_use_lifetime
  substance_use_lifetime$coreID <- gsub("sub-", "", substance_use_lifetime$IID)

  # ==================
  # load environment exposure
  # ==================
  environment_cov_definition = read.csv(paste0(proj_dir, "data/psu/abcd.6.0.environment.definition.time.final.08.26.2025.csv")) 
  environment_cov_definition$coreID <- gsub("sub-", "", environment_cov_definition$participant_id)
  
  # --------------- combine all participants based on coreID, long format
  pgs_environ_psu = environment_cov_definition %>% 
    left_join(Basic_information[, c("height_prs", "genetic_ancestry", "sud_prs2", "aud_prs2", "cod_prs2", "cud_prs2", "nid_prs2", "oud_prs2", "aud_uniq_prs2", 
                                    "cod_uniq_prs2", "cud_uniq_prs2", "nid_uniq_prs2", "oud_uniq_prs2", "geno_unrelated", "coreID")], by = "coreID") %>%
    left_join(substance_use_lifetime, by = c("coreID"))
  
  
  # --------------- included participants, have psu2 definition, genetic_unrelated, three ancestry
  pgs_environ_psu$analysis_sample_y0_5_all = 0
  pgs_environ_psu$analysis_sample_y0_5_all[which(!is.na(pgs_environ_psu$psu2_0_5) & !is.na(pgs_environ_psu$genetic_ancestry) &
                                                   (pgs_environ_psu$genetic_ancestry!="EAS") & (pgs_environ_psu$genetic_ancestry!="SAS"))] = 1
  
  pgs_environ_psu$analysis_sample_y0_6_5_all = 0
  pgs_environ_psu$analysis_sample_y0_6_5_all[which(!is.na(pgs_environ_psu$psu2_0_6_5) & !is.na(pgs_environ_psu$genetic_ancestry) &
                                                     (pgs_environ_psu$genetic_ancestry!="EAS") & (pgs_environ_psu$genetic_ancestry!="SAS"))] = 1
  
  # ==================
  # load transition, survived time
  # ==================
  transition_time = read.csv(paste0(proj_dir, "data/psu/abcd.psu.transition.survived.time.y0_6.NA.03.05.2026.csv"))
  transition_time$coreID <- gsub("sub-", "", transition_time$participant_id)
  pgs_environ_psu = pgs_environ_psu %>%
    left_join(transition_time %>% dplyr::select(c("first_single", "first_multi", "non_to_single", 
                                           "non_to_multi", "coreID")), by = "coreID")
  
  # ---------------  from age to calculate the survival time
  datadir = "C:/Users/DZ1011/Partners HealthCare Dropbox/Dongmei Zhi/abcd6.0.tsv/"
  age_M_file = read_tsv(paste0(datadir, "su_y_mypi.tsv"));
  age_M_wide <- age_M_file %>%
    dplyr::select(participant_id, session_id, su_y_mypi_age) %>%
    mutate(varname = paste0("Age_", session_id)) %>%
    pivot_wider(id_cols = participant_id, names_from = varname, values_from = su_y_mypi_age)
  
  age_M_wide$coreID <- gsub("sub-", "", age_M_wide$participant_id)
  pgs_environ_psu = pgs_environ_psu %>%
    left_join(age_M_wide %>% dplyr::select(-"participant_id"), by = "coreID")
  
  # ---------from NSU to PSU
  pgs_environ_psu = pgs_environ_psu %>%
    rowwise() %>%
    mutate(NSU2PSU = case_when(
      non_to_multi == 0 ~ NA, non_to_multi == 0.5 ~ `Age_ses-00M` - Age_ses.00A,
      non_to_multi == 1 ~ Age_ses.01A - Age_ses.00A, non_to_multi == 1.5 ~ `Age_ses-01M` - Age_ses.00A,
      non_to_multi == 2 ~ Age_ses.02A - Age_ses.00A, non_to_multi == 2.5 ~ `Age_ses-02M` - Age_ses.00A,
      non_to_multi == 3 ~ Age_ses.03A - Age_ses.00A, non_to_multi == 3.5 ~ `Age_ses-03M` - Age_ses.00A,
      non_to_multi == 4 ~ Age_ses.04A - Age_ses.00A, non_to_multi == 4.5 ~ `Age_ses-04M` - Age_ses.00A,
      non_to_multi == 5 ~ Age_ses.05A - Age_ses.00A, non_to_multi == 5.5 ~ `Age_ses-05M` - Age_ses.00A,
      non_to_multi == 6 ~ Age_ses.06A - Age_ses.00A, 
      is.na(non_to_multi) & psu2_0_6==0 ~ Age_ses.06A - Age_ses.00A,
      is.na(non_to_multi) & is.na(psu2_0_6) & psu2_0_5.5==0 ~ `Age_ses-05M` - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_multi, psu2_0_6, psu2_0_5.5))) == 3) & psu2_0_5==0 ~ Age_ses.05A - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_multi, psu2_0_6, psu2_0_5.5, psu2_0_5))) == 4) & psu2_0_4.5==0 ~ `Age_ses-04M` - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_multi, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5))) == 5) & psu2_0_4==0 ~ Age_ses.04A - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_multi, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4))) == 6) & psu2_0_3.5==0 ~ `Age_ses-03M` - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_multi, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5))) == 7) & psu2_0_3==0 ~ Age_ses.03A - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_multi, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5, psu2_0_3))) == 8) & psu2_0_2.5==0 ~ `Age_ses-02M` - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_multi, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5, psu2_0_3, psu2_0_2.5))) == 9) & psu2_0_2==0 ~ Age_ses.02A - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_multi, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5, psu2_0_3, psu2_0_2.5, psu2_0_2))) == 10) & psu2_0_1.5==0 ~ `Age_ses-01M` - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_multi, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5, psu2_0_3, psu2_0_2.5, psu2_0_2, psu2_0_1.5))) == 11) & psu2_0_1==0 ~ Age_ses.01A - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_multi, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5, psu2_0_3, psu2_0_2.5, psu2_0_2, psu2_0_1.5, psu2_0_1))) == 12) & psu2_0_0.5==0 ~ `Age_ses-00M` - Age_ses.00A,
      
      TRUE ~ NA_real_ 
    )) %>%  ungroup()
  
  # ---------from NSU to SSU
  pgs_environ_psu = pgs_environ_psu %>%
    rowwise() %>%
    mutate(NSU2SSU = case_when(
      non_to_single == 0 ~ NA, non_to_single == 0.5 ~ `Age_ses-00M` - Age_ses.00A,
      non_to_single == 1 ~ Age_ses.01A - Age_ses.00A, non_to_single == 1.5 ~ `Age_ses-01M` - Age_ses.00A,
      non_to_single == 2 ~ Age_ses.02A - Age_ses.00A, non_to_single == 2.5 ~ `Age_ses-02M` - Age_ses.00A,
      non_to_single == 3 ~ Age_ses.03A - Age_ses.00A, non_to_single == 3.5 ~ `Age_ses-03M` - Age_ses.00A,
      non_to_single == 4 ~ Age_ses.04A - Age_ses.00A, non_to_single == 4.5 ~ `Age_ses-04M` - Age_ses.00A,
      non_to_single == 5 ~ Age_ses.05A - Age_ses.00A, non_to_single == 5.5 ~ `Age_ses-05M` - Age_ses.00A,
      non_to_single == 6 ~ Age_ses.06A - Age_ses.00A,
      
      is.na(non_to_single) & psu2_0_6==0 ~ Age_ses.06A - Age_ses.00A,
      is.na(non_to_single) & is.na(psu2_0_6) & psu2_0_5.5==0 ~ `Age_ses-05M` - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_single, psu2_0_6, psu2_0_5.5))) == 3) & psu2_0_5==0 ~ Age_ses.05A - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_single, psu2_0_6, psu2_0_5.5, psu2_0_5))) == 4) & psu2_0_4.5==0 ~ `Age_ses-04M` - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_single, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5))) == 5) & psu2_0_4==0 ~ Age_ses.04A - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_single, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4))) == 6) & psu2_0_3.5==0 ~ `Age_ses-03M` - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_single, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5))) == 7) & psu2_0_3==0 ~ Age_ses.03A - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_single, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5, psu2_0_3))) == 8) & psu2_0_2.5==0 ~ `Age_ses-02M` - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_single, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5, psu2_0_3, psu2_0_2.5))) == 9) & psu2_0_2==0 ~ Age_ses.02A - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_single, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5, psu2_0_3, psu2_0_2.5, psu2_0_2))) == 10) & psu2_0_1.5==0 ~ `Age_ses-01M` - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_single, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5, psu2_0_3, psu2_0_2.5, psu2_0_2, psu2_0_1.5))) == 11) & psu2_0_1==0 ~ Age_ses.01A - Age_ses.00A,
      (rowSums(is.na(data.frame(non_to_single, psu2_0_6, psu2_0_5.5, psu2_0_5, psu2_0_4.5, psu2_0_4, psu2_0_3.5, psu2_0_3, psu2_0_2.5, psu2_0_2, psu2_0_1.5, psu2_0_1))) == 12) & psu2_0_0.5==0 ~ `Age_ses-00M` - Age_ses.00A,
      
      TRUE ~ NA_real_
    )) %>%  ungroup()
  
  # ---------from SSU to PSU
  pgs_environ_psu <- pgs_environ_psu %>%
    mutate(
      # 1. 先定义 non -> single 的起始年龄
      single_age = case_when(
        non_to_single == 0   ~ Age_ses.00A,
        non_to_single == 0.5 ~ `Age_ses-00M`,
        non_to_single == 1   ~ Age_ses.01A,
        non_to_single == 1.5 ~ `Age_ses-01M`,
        non_to_single == 2   ~ Age_ses.02A,
        non_to_single == 2.5 ~ `Age_ses-02M`,
        non_to_single == 3   ~ Age_ses.03A,
        non_to_single == 3.5 ~ `Age_ses-03M`,
        non_to_single == 4   ~ Age_ses.04A,
        non_to_single == 4.5 ~ `Age_ses-04M`,
        non_to_single == 5   ~ Age_ses.05A,
        non_to_single == 5.5 ~ `Age_ses-05M`,
        non_to_single == 6   ~ Age_ses.06A,
        TRUE ~ NA_real_
      ),
      
      # 2. 默认值：如果NSU2PSU > NSU2SSU，则直接差值
      SSU2PSU_default = case_when(
        is.finite(NSU2SSU) & is.finite(NSU2PSU) & NSU2PSU > NSU2SSU ~ NSU2PSU - NSU2SSU,
        TRUE ~ NA_real_
      ),
      
      # 3. 按优先级选择 PSU 发生时点，并计算跨度
      SSU2PSU = case_when(
        psu2_0_6 == 1 ~ Age_ses.06A - single_age,
        is.na(psu2_0_6) & psu2_0_5.5 == 1 ~ `Age_ses-05M` - single_age,
        is.na(psu2_0_6) & is.na(psu2_0_5.5) & psu2_0_5 == 1 ~ Age_ses.05A - single_age,
        TRUE ~ SSU2PSU_default
      )
    ) %>%
    dplyr::select(-single_age, -SSU2PSU_default)
  
  # --------------- included participants in survival time analysis, censoring
  pgs_environ_psu$analysis_sample_y0_6_5_survive = 0
  pgs_environ_psu$analysis_sample_y0_6_5_survive[which((!is.na(pgs_environ_psu$NSU2SSU) | !is.na(pgs_environ_psu$SSU2PSU) | !is.na(pgs_environ_psu$NSU2PSU)) & 
                                                       !is.na(pgs_environ_psu$genetic_ancestry) &
                                                     (pgs_environ_psu$genetic_ancestry!="EAS") & (pgs_environ_psu$genetic_ancestry!="SAS"))] = 1
  
  # --------------- flag for final flag in survival analysis
  pgs_environ_psu$psu2_0_6_5_survive = NA
  pgs_environ_psu$psu2_0_6_5_survive[which(pgs_environ_psu$analysis_sample_y0_6_5_survive==1)] = pgs_environ_psu$psu_sum_0_6[which(pgs_environ_psu$analysis_sample_y0_6_5_survive==1)]
  pgs_environ_psu$psu2_0_6_5_survive[which(pgs_environ_psu$psu2_0_6_5_survive>2)] = 2
  
  #  --------------- need based on the updated psu2_0_6_5_survive
  pgs_environ_psu$SSU2PSU_flag = ifelse(is.finite(pgs_environ_psu$SSU2PSU) & pgs_environ_psu$psu2_0_6_5>1, 1, NA)
  pgs_environ_psu$SSU2PSU_flag[which(pgs_environ_psu$psu2_0_6_5==1 & pgs_environ_psu$psu_sum_0_5.5!=0)] = 0
  
  # ==================
  # add covariate of INR and family history of substance use
  # ==================
  INR = read.csv(paste0(proj_dir, "data/psu/Basic_info_prenatal.complete.mean.thk.04.29.2025.csv")) 
  INR = INR[which(INR$eventname=="BL"), ]
  INR$coreID = gsub("NDAR_INV", "", INR$src_subject_id)
  pgs_environ_psu = pgs_environ_psu %>%
    left_join(INR[, c("coreID", "INR.BL")], by = c("coreID"))

  family_history = read_tsv(paste0("C:/Users/DZ1011/Partners HealthCare Dropbox/Dongmei Zhi/abcd6.0.tsv/mh_p_famhx.tsv"))
  family_history = family_history[which(family_history$session_id=="ses-00A"), c("participant_id", "mh_p_famhx__alc_001")]

  pgs_environ_psu = pgs_environ_psu %>%
    left_join(family_history, by = c("participant_id"))
  
  # ==================
  # subtract the substance use from CBCL total scores
  # ==================
  datadir = "C:/Users/DZ1011/Partners HealthCare Dropbox/Dongmei Zhi/abcd6.0.tsv/"
  substance_item_CBCL = read_tsv(paste0(datadir, "mh_p_cbcl.tsv"))
  substance_item_CBCL = substance_item_CBCL %>%
    dplyr::select(participant_id, session_id, mh_p_cbcl__rule_001, mh_p_cbcl__rule_005, mh_p_cbcl__rule_006)
  
  substance_item_CBCL$sum_substance_CBCL <- rowSums(substance_item_CBCL[, c("mh_p_cbcl__rule_001",  "mh_p_cbcl__rule_005",  "mh_p_cbcl__rule_006")], na.rm = TRUE)
  
  substance_wide <- substance_item_CBCL %>%
    dplyr::select(participant_id, session_id, sum_substance_CBCL) %>%
    pivot_wider(names_from = session_id, values_from = sum_substance_CBCL, names_prefix = "substance_item_CBCL_")
  
  # subtract the substance use from CBCL total scores
  pgs_environ_psu <- pgs_environ_psu %>% left_join(substance_wide, by = c("participant_id"))
  
  sessions <- c("00A","01A","02A","03A","04A","05A","06A")
  for(s_time in sessions){
    pgs_environ_psu[[paste0("Total.problems_minus_substance_", s_time)]] <-
      pgs_environ_psu[[paste0("Total.problems_r_ses.", s_time)]] -
      pgs_environ_psu[[paste0("substance_item_CBCL_ses-", s_time)]]
  }
  
  # --------------- save into csv file
  outdir = "C:/Users/DZ1011/Dropbox/2024.PSU.lancet/data/psu/" 
  write.csv(pgs_environ_psu, file= paste0(outdir, "abcd.6.0.pgs.environment.psu.cross.sectional.NA.03.05.2025.csv"), row.names = FALSE)
}
