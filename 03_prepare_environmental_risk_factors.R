
rm(list=ls());

library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(mice)

# ========================================================
# Definition for environment variables
# ========================================================
environment_variable <- list(
  list(var = "Age", file = "ab_p_demo", column = "ab_p_demo_age", index = c(), cat = ""),
  list(var = "Sex", file = "ab_g_stc", column = "ab_g_stc__cohort_sex", index = c(), cat = ""),
  list(var = "Race", file = "ab_g_stc", column = "ab_g_stc__cohort_ethnrace__leg", index = c(), cat = ""),
  list(var = "Handedness", file = "nc_y_ehis", column = "nc_y_ehis_score", index = c(), cat = ""),
  list(var = "Height", file = "ph_y_anthr", column = "ph_y_anthr__height_mean", index = c(), cat = ""),
  list(var = "Weight", file = "ph_y_anthr", column = "ph_y_anthr__weight_mean", index = c(), cat = ""),
  list(var = "Waist", file = "ph_y_anthr", column = "ph_y_anthr__waist_001", index = c(), cat = ""),
  list(var = "Pubertal", file = c("ph_p_pds", "ph_y_pds"), column = c("ph_p_pds__f_mean", "ph_p_pds__m_mean", "ph_y_pds__f_mean", "ph_y_pds__m_mean"), index = c(), cat = "mean"),
  list(var = "Site", file = "ab_g_dyn", column = "ab_g_dyn__design_site", index = c(), cat = ""),
  list(var = "Scanner", file = "ab_g_dyn", column = "ab_g_dyn__design_mr__manufact", index = c(), cat = ""),
  list(var = "FamilyID", file = "ab_g_stc", column = "ab_g_stc__design_id__fam", index = c(), cat = ""),
  
  # ---------------------- image quality control
  list(var = "mr_y_qc__incl__smri__t1_indicator", file = "mr_y_qc__incl", column = "mr_y_qc__incl__smri__t1_indicator", index = c(), cat = ""),
  list(var = "mr_y_qc__incl__smri__t2_indicator", file = "mr_y_qc__incl", column = "mr_y_qc__incl__smri__t2_indicator", index = c(), cat = ""),
  list(var = "mr_y_qc__incl__dmri_indicator", file = "mr_y_qc__incl", column = "mr_y_qc__incl__dmri_indicator", index = c(), cat = ""),
  list(var = "mr_y_qc__incl__rsfmri_indicator", file = "mr_y_qc__incl", column = "mr_y_qc__incl__rsfmri_indicator", index = c(), cat = ""),
  list(var = "mr_y_qc__incl__tfmri__mid_indicator", file = "mr_y_qc__incl", column = "mr_y_qc__incl__tfmri__mid_indicator", index = c(), cat = ""),
  list(var = "mr_y_qc__incl__tfmri__nback_indicator", file = "mr_y_qc__incl", column = "mr_y_qc__incl__tfmri__nback_indicaor", index = c(), cat = ""),
  list(var = "mr_y_qc__incl__tfmri__sst_indicator", file = "mr_y_qc__incl", column = "mr_y_qc__incl__tfmri__sst_indicator", index = c(), cat = ""),
  
  list(var = "mr_y_qc__mot__dmri__mot_mean", file = "mr_y_qc__mot", column = "mr_y_qc__mot__dmri__mot_mean", index = c(), cat = ""),
  list(var = "mr_y_qc__mot__rsfmri__mot_mean", file = "mr_y_qc__mot", column = "mr_y_qc__mot__rsfmri__mot_mean", index = c(), cat = ""),
  list(var = "mr_y_qc__mot__tfmri__mot_mean", file = "mr_y_qc__mot", column = "mr_y_qc__mot__tfmri__mot_mean", index = c(), cat = ""),
  list(var = "mr_y_qc__mot__tfmri__nback__mot_mean", file = "mr_y_qc__mot", column = "mr_y_qc__mot__tfmri__nback__mot_mean", index = c(), cat = ""),
  list(var = "mr_y_qc__mot__tfmri__sst__mot_mean", file = "mr_y_qc__mot", column = "mr_y_qc__mot__tfmri__sst__mot_mean", index = c(), cat = ""),
  list(var = "mr_y_smri__vol__aseg__icv_sum", file = "mr_y_smri__vol__aseg", column = "mr_y_smri__vol__aseg__icv_sum", index = c(), cat = ""),
  
  # ---------------------- environment exposures
  list(var = "Maternal_age", file = "ph_p_dhx", column = "ph_p_dhx_003__01", index = c(), cat = ""),
  list(var = "Paternal_age", file = "ph_p_dhx", column = "ph_p_dhx_004__01", index = c(), cat = ""),
  list(var = "Twin_birth", file = "ph_p_dhx", column = "ph_p_dhx_005", index = c(), cat = ""),
  list(var = "Planned_pregnancy", file = "ph_p_dhx", column = "ph_p_dhx_006", index = c(), cat = ""),
  list(var = "Premature_birth", file = "ph_p_dhx", column = "ph_p_dhx__birth_001", index = c(), cat = ""), #  -----1
  list(var = "Weeks_premature", file = "ph_p_dhx", column = "ph_p_dhx__birth_001__01", index = c(), cat = ""),
  list(var = "Cesarean_delivery", file = "ph_p_dhx", column = "ph_p_dhx__birth_002", index = c(), cat = ""),
  list(var = "Birth_complication", file = "ph_p_dhx", column = c("ph_p_dhx__birth_003", "ph_p_dhx__birth_004", "ph_p_dhx__birth_005", 
             "ph_p_dhx__birth_006", "ph_p_dhx__birth_007", "ph_p_dhx__birth_008", "ph_p_dhx__birth_009", "ph_p_dhx__birth_010"), index = c(), cat = "logic_sum"),
  list(var = "Months_breastfed", file = "ph_p_dhx", column = "ph_p_dhx_012", index = c(), cat = ""),
  list(var = "Birth_weight", file = "ph_p_dhx", column = "ph_p_dhx_002__01", index = c(), cat = ""),
  list(var = "Delayed_motor_development", file = "ph_p_dhx", column = "ph_p_dhx_013", index = c(), cat = ""),
  list(var = "Delayed_verbal_development", file = "ph_p_dhx", column = "ph_p_dhx_014", index = c(), cat = ""),
  list(var = "Prenatal_other_substance_use", file = "ph_p_dhx", column = c("ph_p_dhx__coc_001a", "ph_p_dhx__opi_001a", "ph_p_dhx__rxpain_001a", 
             "ph_p_dhx__nic_001b", "ph_p_dhx__alc_001b", "ph_p_dhx__mj_001b", "ph_p_dhx__coc_001b", 
             "ph_p_dhx__opi_001b", "ph_p_dhx__rxpain_001b"), index = c(), cat = "logic_sum"),  
  list(var = "Prenatal_tobacco_use", file = "ph_p_dhx", column = "ph_p_dhx__nic_001a", index = c(), cat = ""),
  list(var = "Prenatal_alcohol_use", file = "ph_p_dhx", column = "ph_p_dhx__alc_001a", index = c(), cat = ""),
  list(var = "Prenatal_marijuana_use", file = "ph_p_dhx", column = "ph_p_dhx__mj_001a", index = c(), cat = ""),

  list(var = "Prenatal_medication_conditions", file = "ph_p_dhx", column = c("ph_p_dhx__med_001", "ph_p_dhx__med_002", "ph_p_dhx__med_003", 
             "ph_p_dhx__med_004", "ph_p_dhx__med_005", "ph_p_dhx__med_006", "ph_p_dhx__med_007", "ph_p_dhx__med_008", "ph_p_dhx__med_009", 
             "ph_p_dhx__med_010", "ph_p_dhx__med_011", "ph_p_dhx__med_012", "ph_p_dhx__med_013"), index = c(), cat = "logic_sum"), 
  list(var = "Traumatic_events_parent", file = "mh_p_ksads__ptsd", column = "mh_p_ksads__ptsd__trma__evnt__pres_sx", index = c(), cat = ""),
  
  list(var = "Suicidal_ideation_youth", file = "mh_y_ksads__suic", column = "mh_y_ksads__suic__idea__pres_sx", index = c(), cat = ""),
  
  list(var = "Cyberbullying", file = "mh_y_cb", column = "mh_y_cb_001a", index = c(), cat = ""),
  list(var = "Serious_medical_history", file = "ph_p_mhx", column = c("ph_p_mhx__er_001__01", "ph_p_mhx__er__lt_001__01"), index = c(), cat = "sum"),
  list(var = "Religiosity", file = "ab_p_demo", column = "ab_p_demo__relig_002", index = c(), cat = ""),
  list(var = "Physical_activity", file = "ph_y_pa", column = "ph_y_pa_001", index = c(), cat = ""),
  list(var = "Physical_health", file = "nt_y_fitbq", column = "nt_y_fitbq__pre_001", index = c(), cat = ""),
  list(var = "Sedentary_health", file = "nt_y_fitbq", column = "nt_y_fitbq__pre_003", index = c(), cat = ""),
  list(var = "Sleep_health", file = "nt_y_fitbq", column = "nt_y_fitbq__pre_005", index = c(), cat = ""),
  list(var = "Sleep_duration_youth", file = "ph_y_mctq", column = "ph_y_mctq__sleep_dur", index = c(), cat = ""),
  
  list(var = "TV_screen_use", file = "nt_y_stq", column = c("nt_y_stq__screen__wkdy_001", "nt_y_stq__screen__wknd_001"), index = c(), cat = "sum"),
  list(var = "Video_screen_use", file = "nt_y_stq", column = c("nt_y_stq__screen__wkdy_002", "nt_y_stq__screen__wknd_002"), index = c(), cat = "sum"), 
  list(var = "Games_screen_use", file = "nt_y_stq", column = c("nt_y_stq__screen__wkdy_003", "nt_y_stq__screen__wknd_003"), index = c(), cat = "sum"),
  list(var = "Text_screen_use", file = "nt_y_stq", column = c("nt_y_stq__screen__wkdy_004", "nt_y_stq__screen__wknd_004"), index = c(), cat = "sum"),
  list(var = "Social_media_screen_use", file = "nt_y_stq", column = c("nt_y_stq__screen__wkdy_005", "nt_y_stq__screen__wknd_005"), index = c(), cat = "sum"),
  list(var = "Chat_screen_use", file = "nt_y_stq", column = c("nt_y_stq__screen__wkdy_006", "nt_y_stq__screen__wknd_006"), index = c(), cat = "sum"),
  list(var = "Mature_games_screen_use", file = "nt_y_stq", column = "nt_y_stq__screen__mature_001", index = c(), cat = ""),
  list(var = "R_rated_screen_use", file = "nt_y_stq", column = "nt_y_stq__screen__mature_002", index = c(), cat = ""),
  
  list(var = "Life_events_parent", file = "mh_p_ple", column = "mh_p_ple_count", index = c(), cat = ""),
  
  list(var = "Life_event_severity_parent", file = "mh_p_ple", column = "mh_p_ple__severity_mean", index = c(), cat = ""),
  

  list(var = "Life_events_youth", file = "mh_y_ple", column = "mh_y_ple_count", index = c(), cat = ""),
  
  list(var = "Life_event_severity_youth", file = "mh_y_ple", column = "mh_y_ple__severity_mean", index = c(), cat = ""),
  
  list(var = "Family_conflict_parent", file = "fc_p_fes", column = "fc_p_fes__confl_mean", index = c(), cat = ""),
  list(var = "Family_conflict_youth", file = "fc_y_fes", column = "fc_y_fes__confl_mean", index = c(), cat = ""),
  list(var = "Family_intellectual_culture", file = "fc_p_fes", column = "fc_p_fes__intelcult_mean", index = c(), cat = ""),
  list(var = "Family_activity_recreational", file = "fc_p_fes", column = "fc_p_fes__rec_mean", index = c(), cat = ""),
  list(var = "Family_organization", file = "fc_p_fes", column = "fc_p_fes__org_mean", index = c(), cat = ""),
  list(var = "Family_cohesion_youth", file = "fc_y_fes", column = "fc_y_fes__cohes_mean", index = c(), cat = ""),
  list(var = "Family_cohesion_parent", file = "fc_p_fes", column = "fc_p_fes__cohes_mean", index = c(), cat = ""),
  list(var = "Family_expression", file = "fc_p_fes", column = "fc_p_fes__expr_mean", index = c(), cat = ""),
  list(var = "Family_support", file = "fc_y_vs", column = "fc_y_vs__supp_mean", index = c(), cat = ""), 
  list(var = "Parental_monitoring", file = "fc_p_pk", column = "fc_p_pk__knowl_mean", index = c(), cat = ""),
  list(var = "Primary_caregiver_warmth", file = "fc_y_crpbi", column = "fc_y_crpbi__cg1_mean", index = c(), cat = ""),
  list(var = "Secondary_caregiver_warmth", file = "fc_y_crpbi", column = "fc_y_crpbi__cg2_mean", index = c(), cat = ""),
  list(var = "Family_income", file = "ab_p_demo", column = "ab_p_demo__income__hhold_001", index = c(), cat = ""),
  list(var = "Caregiver_education", file = "ab_p_demo", column = c("ab_p_demo__edu__slf_001", "ab_p_demo__edu__prtnr_001"), index = c(), cat = "mean"),
  list(var = "Caregiver_employment", file = "ab_p_demo", column = c("ab_p_demo__empl__slf_001", "ab_p_demo__empl__prtnr_001"), index = c(), cat = "mean"),
  list(var = "Number_of_people_living", file = "ab_p_demo", column = "ab_p_demo__roster_001", index = c(), cat = ""),
  list(var = "Caregiver_marital_status", file = "ab_p_demo", column = "ab_p_demo__marital__slf_001", index = c(), cat = ""),
  list(var = "Severe_financial_difficulty", file = "ab_p_demo", column = c("ab_p_demo__exp__fam_001", "ab_p_demo__exp__fam_002", "ab_p_demo__exp__fam_003", 
             "ab_p_demo__exp__fam_004", "ab_p_demo__exp__fam_005", "ab_p_demo__exp__fam_006", "ab_p_demo__exp__fam_007"), index = c(), cat = "logic_sum"),
  list(var = "Parental_alcohol_issues", file = "mh_p_famhx", column = "mh_p_famhx__alc_001", index = c(), cat = ""),
  list(var = "Parental_drug_issues", file = "mh_p_famhx", column = "mh_p_famhx__drg_001", index = c(), cat = ""),
  list(var = "Parental_depression_issues", file = "mh_p_famhx", column = "mh_p_famhx__dep_001", index = c(), cat = ""),
  list(var = "Parental_mania_issues", file = "mh_p_famhx", column = "mh_p_famhx__mania_001", index = c(), cat = ""),
  list(var = "Parental_visions_issues", file = "mh_p_famhx", column = "mh_p_famhx__halluc_001", index = c(), cat = ""),
  list(var = "Parental_trouble_issues", file = "mh_p_famhx", column = "mh_p_famhx__troub_001", index = c(), cat = ""),
  list(var = "Parental_nerves_issues", file = "mh_p_famhx", column = "mh_p_famhx__nerve_001", index = c(), cat = ""),
  list(var = "Parental_emotion_issues", file = "mh_p_famhx", column = "mh_p_famhx__doc_001", index = c(), cat = ""),
  list(var = "Parental_hospital_issues", file = "mh_p_famhx", column = "mh_p_famhx__hosp_001", index = c(), cat = ""),
  list(var = "Parental_suicide_issues", file = "mh_p_famhx", column = "mh_p_famhx__suic_001", index = c(), cat = ""),
  list(var = "Parental_supervision", file = "fc_y_mnbs", column = "fc_y_mnbs__superv_mean", index = c(), cat = ""),
  list(var = "Substance_availiability", file = "su_p_crpf", column = c("su_p_crpf__alc_001", "su_p_crpf__nic__cig_001", "su_p_crpf__nic__vape_001", 
             "su_p_crpf__mj_001", "su_p_crpf__illdrg_001"), index = c(), cat = "mean"),
  list(var = "Neighborhood_security_parent", file = "fc_p_nsc", column = "fc_p_nsc__ns_mean", index = c(), cat = ""),
  list(var = "Neighborhood_security_youth", file = "fc_y_nsc", column = "fc_y_nsc__ns_003", index = c(), cat = ""),
  list(var = "Neighborhood_crime", file = "le_l_crime", column = "le_l_crime__addr1_count", index = c(), cat = ""),
  list(var = "Neighborhood_collective_capacity", file = "fc_p_nce", column = "fc_p_nce_mean", index = c(), cat = ""),
  list(var = "Area_deprivation_index", file = "le_l_adi", column = "le_l_adi__addr1__national_prcnt", index = c(), cat = ""),
  list(var = "Pollution", file = "le_l_pm25", column = "le_l_pm25__addr1__pm25_mean__2016", index = c(), cat = ""),
  list(var = "Lead_risk", file = "le_l_leadrisk", column = "le_l_leadrisk__addr1_idx", index = c(), cat = ""),
  list(var = "Noise", file = "le_l_noise", column = "le_l_noise__addr1__24hmean__total_soundlvl", index = c(), cat = ""),
  list(var = "School_environment", file = "fc_y_srpf", column = "fc_y_srpf__env_mean", index = c(), cat = ""),
  list(var = "Positive_school_involvement", file = "fc_y_srpf", column = "fc_y_srpf__involv_mean", index = c(), cat = ""),
  list(var = "School_disengagement", file = "fc_y_srpf", column = "fc_y_srpf__dis_mean", index = c(), cat = ""),
  list(var = "Close_friends", file = "mh_y_resil", column = "mh_y_resil__close_sum", index = c(), cat = ""),
  
  list(var = "Peer_influence", file = "fc_y_rpi", column = "fc_y_rpi_mean", index = c(), cat = ""),
  list(var = "Peer_relation_victimization", file = "mh_y_peq", column = "mh_y_peq__rel__vict_sum", index = c(), cat = ""),
  list(var = "Peer_reputation_aggression", file = "mh_y_peq", column = "mh_y_peq__rep__agg_sum", index = c(), cat = ""),
  list(var = "Peer_reputation_victimization", file = "mh_y_peq", column = "mh_y_peq__rep__vict_sum", index = c(), cat = ""),
  list(var = "Peer_overt_aggression", file = "mh_y_peq", column = "mh_y_peq__overt__agg_sum", index = c(), cat = ""),
  list(var = "Peer_overt_victimization", file = "mh_y_peq", column = "mh_y_peq__overt__vict_sum", index = c(), cat = ""),
  list(var = "Peer_relational_aggression", file = "mh_y_peq", column = "mh_y_peq__rel__agg_sum", index = c(), cat = ""),
  list(var = "Peer_tobacco_use", file = "su_y_pgd", column = "su_y_pgd_001", index = c(), cat = ""),
  list(var = "Peer_alcohol_use", file = "su_y_pgd", column = "su_y_pgd_002", index = c(), cat = ""),
  
  # ---------------------- CBCL behaviors
  list(var = "Anxious/depressed_r", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__anxdep_sum", index = c(), cat = ""),
  list(var = "Anxious/depressed_t", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__anxdep_tscore", index = c(), cat = ""),
  list(var = "Withdrawn/depressed_r", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__wthdep_sum", index = c(), cat = ""),
  list(var = "Withdrawn/depressed_t", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__wthdep_tscore", index = c(), cat = ""),
  list(var = "Somatic complaints_r", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__som_sum", index = c(), cat = ""),
  list(var = "Somatic complaints_t", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__som_tscore", index = c(), cat = ""),
  list(var = "Social problems_r", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__soc_sum", index = c(), cat = ""),
  list(var = "Social problems_t", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__soc_tscore", index = c(), cat = ""),
  list(var = "Thought problems_r", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__tho_sum", index = c(), cat = ""),
  list(var = "Thought problems_t", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__tho_tscore", index = c(), cat = ""),
  list(var = "Attention problems_r", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__attn_sum", index = c(), cat = ""),
  list(var = "Attention problems_t", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__attn_tscore", index = c(), cat = ""),
  list(var = "Rule-breaking behaviors_r", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__rule_sum", index = c(), cat = ""),
  list(var = "Rule-breaking behaviors_t", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__rule_tscore", index = c(), cat = ""),
  list(var = "Aggressive behaviors_r", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__aggr_sum", index = c(), cat = ""),
  list(var = "Aggressive behaviors_t", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__aggr_tscore", index = c(), cat = ""),
  list(var = "Internalizing problems_r", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__int_sum", index = c(), cat = ""),
  list(var = "Internalizing problems_t", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__int_tscore", index = c(), cat = ""),
  list(var = "Externalizing problems_r", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__ext_sum", index = c(), cat = ""),
  list(var = "Externalizing problems_t", file = "mh_p_cbcl", column = "mh_p_cbcl__synd__ext_tscore", index = c(), cat = ""),
  list(var = "Total problems_r", file = "mh_p_cbcl", column = "mh_p_cbcl_sum", index = c(), cat = ""),
  list(var = "Total problems_t", file = "mh_p_cbcl", column = "mh_p_cbcl_tscore", index = c(), cat = ""),
  list(var = "Depressive symptoms_r", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__dep_sum", index = c(), cat = ""),
  list(var = "Depressive symptoms_t", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__dep_tscore", index = c(), cat = ""),
  list(var = "Anxiety symptoms_r", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__anx_sum", index = c(), cat = ""),
  list(var = "Anxiety symptoms_t", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__anx_tscore", index = c(), cat = ""),
  list(var = "Somatic symptoms_r", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__somat_sum", index = c(), cat = ""),
  list(var = "Somatic symptoms_t", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__somat_tscore", index = c(), cat = ""),
  list(var = "ADHD symptoms_r", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__adhd_sum", index = c(), cat = ""),
  list(var = "ADHD symptoms_t", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__adhd_tscore", index = c(), cat = ""),
  list(var = "Oppositional defiant symptoms_r", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__opp_sum", index = c(), cat = ""),
  list(var = "Oppositional defiant symptoms_t", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__opp_tscore", index = c(), cat = ""),
  list(var = "Conduct symptoms_r", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__cond_sum", index = c(), cat = ""),
  list(var = "Conduct symptoms_t", file = "mh_p_cbcl", column = "mh_p_cbcl__dsm__cond_tscore", index = c(), cat = "")
)


# ========================================================
# extract each variable from tsv file
# ========================================================
datadir = "C:/Users/DZ1011/Partners HealthCare Dropbox/Dongmei Zhi/abcd6.0.tsv/"
age_file = read_tsv(paste0(datadir,  environment_variable[[1]]$file, ".tsv"));
environment_all = age_file[, c(1,2)]

for (num_evn in (1:length(environment_variable))){
  
  # ---------------- load corresponding tsv file
  if (length(environment_variable[[num_evn]]$file)==1){
    data_file = read_tsv(paste0(datadir,  environment_variable[[num_evn]]$file, ".tsv"));
    
  }else{
    data_file = read_tsv(paste0(datadir,  environment_variable[[num_evn]]$file[1], ".tsv"));
    
    for (name_file in c(environment_variable[[num_evn]]$file[2:length(environment_variable[[num_evn]]$file)])){
      data_add = read_tsv(paste0(datadir,  name_file, ".tsv"));
      data_file = full_join(data_file, data_add, by = c("participant_id", "session_id"))
      
    }
  }
  
  
  if ((length(environment_variable[[num_evn]]$column)==1) & (environment_variable[[num_evn]]$file[1] != "ab_g_stc")){ # only one column, not static
    data_column = data_file[, c("participant_id", "session_id", environment_variable[[num_evn]]$column)]
    environment_all = full_join(environment_all, data_column, by = c("participant_id", "session_id"))
    colnames(environment_all)[ncol(environment_all)] = environment_variable[[num_evn]]$var
    
  }else if (((length(environment_variable[[num_evn]]$column)==1) & (environment_variable[[num_evn]]$file[1] == "ab_g_stc"))){ # only one column, static
    data_column = data_file[, c("participant_id", environment_variable[[num_evn]]$column)]
    environment_all = full_join(environment_all, data_column, by = c("participant_id"))
    colnames(environment_all)[ncol(environment_all)] = environment_variable[[num_evn]]$var
    
  }else if (length(environment_variable[[num_evn]]$column)>1) {
    data_column = data_file[, c("participant_id", "session_id", environment_variable[[num_evn]]$column)]
    
    # remove outliers
    if (typeof(data_column[, 3])=="list"){
      data_column[data_column == "n/a"] <- NA
      data_column[, c(environment_variable[[num_evn]]$column)] <- lapply(data_column[, environment_variable[[num_evn]]$column], as.numeric)
    }
    
    data_column[] <- lapply(data_column, function(col) {
      if (is.numeric(col)) {
        col[col %in% c(222, 444, 555, 666, 777, 888, 999)] <- NA
      }
      return(col)})
    
    if (environment_variable[[num_evn]]$cat == "mean"){
      data_column[, paste0(environment_variable[[num_evn]]$var, "_mean")] = rowMeans(data_column[, environment_variable[[num_evn]]$column], na.rm = TRUE)
      
      # -------------set NA if all is NA
      all_na_rows <- apply(data_column[, environment_variable[[num_evn]]$column], 1, function(x) all(is.na(x)))
      data_column[all_na_rows, paste0(environment_variable[[num_evn]]$var, "_mean")] <- NA
      
      environment_all = full_join(environment_all, data_column[, c("participant_id", "session_id", paste0(environment_variable[[num_evn]]$var, "_mean"))], by = c("participant_id", "session_id"))
      
    }else if (environment_variable[[num_evn]]$cat == "logic_sum"){
      data_column[, paste0(environment_variable[[num_evn]]$var, "_logic")] = rowSums(data_column[, environment_variable[[num_evn]]$column] == 1, na.rm = TRUE)
      data_column[which(data_column[, paste0(environment_variable[[num_evn]]$var, "_logic")]>0), paste0(environment_variable[[num_evn]]$var, "_logic")] = 1
      
      # -------------set NA if all is NA
      all_na_rows <- apply(data_column[, environment_variable[[num_evn]]$column], 1, function(x) all(is.na(x)))
      data_column[all_na_rows, paste0(environment_variable[[num_evn]]$var, "_logic")] <- NA
      
      environment_all = full_join(environment_all, data_column[, c("participant_id", "session_id", paste0(environment_variable[[num_evn]]$var, "_logic"))], by = c("participant_id", "session_id"))
    }else if (environment_variable[[num_evn]]$cat == "sum"){
      data_column[, paste0(environment_variable[[num_evn]]$var, "_sum")] = rowSums(data_column[, environment_variable[[num_evn]]$column], na.rm = TRUE)
      
      # -------------set NA if all is NA
      all_na_rows <- apply(data_column[, environment_variable[[num_evn]]$column], 1, function(x) all(is.na(x)))
      data_column[all_na_rows, paste0(environment_variable[[num_evn]]$var, "_sum")] <- NA
      
      environment_all = full_join(environment_all, data_column[, c("participant_id", "session_id", paste0(environment_variable[[num_evn]]$var, "_sum"))], by = c("participant_id", "session_id"))
      
    }
    rm(all_na_rows)
    
  }
  
  rm(data_file, data_column)
  
}


# ========================================================
# set outlier as NA
# ========================================================
for (num_evn in c(3:dim(environment_all)[2])){
  if (typeof(environment_all[, num_evn])=="list"){
    environment_all[which(environment_all[, num_evn] == "n/a"), num_evn] <- NA
  }
}

environment_all[, c(3:11, 13:dim(environment_all)[2])] <- lapply(environment_all[, c(3:11, 13:dim(environment_all)[2])], as.numeric)  # exclude scanner
environment_all[] <- lapply(environment_all, function(col) {
  if (is.numeric(col)) {
    col[col %in% c(222, 444, 555, 666, 777, 888, 999)] <- NA
  }
  return(col)
  })

outdir = "C:/Users/DZ1011/Dropbox/2024.PSU.lancet/data/psu/" 
save(environment_all, file= paste0(outdir, "abcd.6.0.environment.definition.session.all.08.26.2025.Rdb"))

# ========================================================
# extract environment exposures with the best time point
# ========================================================
{
  environmental_1point = environment_all[, c(1:2, 27:116)]
  
  vars_to_check <- setdiff(names(environmental_1point), c("participant_id", "session_id"))
  
  # count the number of NA and non-NA for each variable
  summary_df <- environmental_1point %>%
    group_by(session_id) %>%
    summarise(across(all_of(vars_to_check),
                     ~sum(!is.na(.)),
                     .names = "{.col}"),
              .groups = "drop")
  
  max_timepoints <- apply(summary_df[,-1], 2, function(x) summary_df$session_id[which.max(x)])
  
  result_list <- lapply(names(max_timepoints), function(varname) {
    max_time <- max_timepoints[[varname]]
    
    # Filter the data of the variable at the maximum time point
    temp <- environmental_1point %>%
      filter(session_id == max_time) %>%
      select(participant_id, value = all_of(varname))
    
    # Renaming variables
    colnames(temp)[2] <- paste0(varname, "_", max_time)
    return(temp)
  })
}

# ========================================================
# extract covariate with all time point
# ========================================================
covariate_wide <- environment_all %>%
  pivot_wider(
    id_cols = participant_id,
    names_from = session_id,
    values_from = c(colnames(environment_all)[c(3:26, 117:150)])
    names_glue = "{.value}_{session_id}"  
  )

# Merge all variable data (by participant_id)
environment_time <- Reduce(function(x, y) full_join(x, y, by = "participant_id"), 
         c(result_list, list(covariate_wide)))
environment_time$`Weeks_premature_ses-00A`[which(is.na(environment_time$`Weeks_premature_ses-00A`) & (environment_time$`Premature_birth_ses-00A`==0))] = 0
environment_time <- environment_time %>%
  select(-`Premature_birth_ses-00A`)

outdir = "C:/Users/DZ1011/Dropbox/2024.PSU.lancet/data/psu/" 
write.csv(environment_time, file= paste0(outdir, "abcd.6.0.environment.definition.time.final.08.26.2025.csv"), row.names = FALSE)

rm(result_list)


# ========================================================
# using VIF to reduce multicollinearity among variables
# ========================================================
if (1==0){
  
  proj_dir = "/local_mount/space/arya/3/users/dmzhi/project/PRS.PSU.environment/"
  
  environment_time = read.csv(paste0(proj_dir, "data/abcd.6.0.environment.definition.time.final.08.15.2025.csv"))
  
  environment_time = environment_time[, c(2:51, 53:90)]
  
  # identify the participants with a high proportion of missing data
  na_ratio <- colSums(is.na(environment_time)) / nrow(environment_time)
  print(na_ratio)
  length(which(na_ratio > 0.3))
  environment_time_vif = environment_time[, which(na_ratio < 0.5)]
  
  # type of data
  data_type = sapply(environment_time_vif, class)
  
  method_vec <- sapply(environment_time_vif, function(col) {
    unique_vals <- unique(col[!is.na(col)]) 
    n_unique <- length(unique_vals)
    
    if (any(is.na(col))) { 
      if (n_unique == 2) {
        return("logreg")   
      } else {
        return("pmm")     
      }
    } else {
      return("")  
    }
  })
  
  # missing data with mice
  num_impute = 10
  data_fea_imputed <- mice(environment_time_vif, m = num_impute, method = method_vec, seed = 123)
  
  save(data_fea_imputed, file = paste0(proj_dir, "experiment/environment.all.time.imputation10.08.15.2025.RData"))
  
  data_fea_imputed_exp <- complete(data_fea_imputed, 1)
  aa = cor(data_fea_imputed_exp[, sapply(data_fea_imputed_exp, is.numeric)], use = "pairwise.complete.obs")
  diag(aa) <- NA
  which_large <- which(aa > 0.8, arr.ind = TRUE) 
  
  #------------------------ using each imputation to perform VIF
  vif_list <- list()
  for (num_imp in c(1:num_impute)){
    data_fea_imputed_exp <- complete(data_fea_imputed, num_imp)
    data_fea_imputed_exp = data_fea_imputed_exp[, 
                                                !(names(data_fea_imputed_exp) %in% c("participant_id", "Screen_use_weekdays_ses-00A", "Screen_use_weekend_ses-00A"))]
    
    calculate_vif_per_variable <- function(df) {
      vif_values <- numeric(ncol(df))
      names(vif_values) <- colnames(df)
      
      for (i in 1:ncol(df)) {
        y <- df[[i]]
        x <- df[,-i]
        fit <- lm(y ~ ., data = x)
        
        R2 <- summary(fit)$r.squared
        vif_values[i] <- 1 / (1 - R2)
      }
      return(vif_values)
    }
    vif_list[[num_imp]] <- calculate_vif_per_variable(data_fea_imputed_exp)
  }
  
  vif_across_all <- do.call(cbind, vif_list)
  colnames(vif_across_all) <- paste0("imp", 1:num_impute)
  vif_mean <- rowMeans(vif_across_all)
  
  high_vif_vars <- rownames(vif_across_all)[vif_mean > 5]
  print(high_vif_vars)
  
  write.csv(vif_across_all, file= paste0(proj_dir, "experiment/VIF.value.environment.exposures.03.17.2026.csv"), row.names = FALSE)
}






