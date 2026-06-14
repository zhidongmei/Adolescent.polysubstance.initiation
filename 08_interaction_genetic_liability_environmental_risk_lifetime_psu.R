
rm(list=ls());

library(car)
library(glmmTMB)
library(dplyr)
library(lmtest)

proj_dir = "/local_mount/space/arya/3/users/dmzhi/project/PRS.PSU.environment/"

# ========================================================
# environment*prs interactions on PSU
# ========================================================
Basic_information = read.csv(paste0(proj_dir, "data/abcd.6.0.pgs.environment.psu.cross.sectional.NA.03.05.2025.csv")) 
Basic_information = Basic_information[which(Basic_information$analysis_sample_y0_6_5_all==1),]
Basic_information[, "psu_ssu_0_6_5"] = Basic_information$psu_0_6_5
Basic_information$psu_ssu_0_6_5 = ifelse(is.na(Basic_information$psu_ssu_0_6_5),  0,  # 原来是 NA 的 → 变成 0
                                ifelse(Basic_information$psu_ssu_0_6_5 == 0, NA, Basic_information$psu_ssu_0_6_5))  # 原来是 0 的 → 变成 NA；其他保持不变

#------------------------ logistic regression with psu as outcome
out = c()

for (prs_thresh in c(3/4)){ #, 3/4
  for (type_race in c("EUR", "AMR", "AFR")){ #, "AMR","AFR"
    
    data_fea = Basic_information[which(Basic_information$genetic_ancestry == type_race), ]
    
    data_fea$prs_group <- cut(data_fea$sud_prs2, breaks = quantile(data_fea$sud_prs2, probs = c(0, prs_thresh, 1), na.rm = TRUE),
                              include.lowest = TRUE, labels = c("Low", "High"))
    
    for (type_psu in c("psu_0_6_5", "ssu_0_6_5", "psu_ssu_0_6_5")){#"psu_0_6_5", "ssu_0_6_5", 
      
      
      num_environ = c(2:51, 53:90)
      
      for (num_fea in num_environ){    
        
        print(num_fea)
        
        data_combine = data.frame(IID = data_fea$participant_id, environ = data_fea[, num_fea], psu_out = data_fea[, type_psu], prs_group = data_fea$prs_group,
                                  age = data_fea$Age_ses.00A, sex = data_fea$Sex_ses.00A, race = data_fea$genetic_ancestry,
                                  site=data_fea$Site_ses.00A, familyID = data_fea$FamilyID_ses.00A, 
                                  Tanner_stage = data_fea$Pubertal_mean_ses.00A, INR = data_fea$INR.BL)
        
        data_combine = na.omit(data_combine[, c('IID', 'environ', 'psu_out', 'prs_group', 'age', 'sex', 'Tanner_stage', 'INR', 'site', 'familyID')])
        
        data_combine$psu_out <- as.numeric(data_combine$psu_out)
        data_combine$sex <- as.factor(data_combine$sex)
        data_combine$site <- as.factor(data_combine$site)
        data_combine$familyID <- as.factor(data_combine$familyID)
        
        data_combine$prs_group <- as.factor(data_combine$prs_group)
        
        data_combine$environ <- as.numeric(scale(data_combine$environ))
        data_combine$age <- as.numeric(scale(data_combine$age))
        data_combine$Tanner_stage <- as.numeric(scale(data_combine$Tanner_stage))
        data_combine$INR <- scale(as.numeric(data_combine$INR))
        
        # adjust the number for each group
        num_count = data_combine %>% count(prs_group, environ) %>% arrange(prs_group, environ)
        
        if (dim(data_combine)[1]>500 & (length(num_count$n) >= 8) || (length(num_count$n) < 8 && min(num_count$n) > 10)){
          
          model <-  glmmTMB(psu_out ~ prs_group*environ + prs_group*(age + sex + Tanner_stage + INR) + environ*(age + sex + Tanner_stage + INR) + 
                              (1 | familyID) + (1 | site), data = data_combine, family = binomial(link = "logit"))
          
          coeff = summary(model)$coefficients$cond
          
          F_stat = car::Anova(model, type = 3)
          
          model_reduced = glmmTMB(psu_out ~ prs_group + environ + prs_group*(age + sex + Tanner_stage + INR) + environ*(age + sex + Tanner_stage + INR) + 
                                    (1 | familyID) + (1 | site), data = data_combine, family = binomial(link = "logit"))
          
          lr_test_stat = lrtest(model_reduced, model)
          
          # extract each term, and interaction term
          out1 = c(round(coeff["environ", c("Estimate", "Std. Error", "z value")], 3), sprintf("%.2e", coeff["environ", "Pr(>|z|)"]),
                   round(coeff["prs_groupHigh", c("Estimate", "Std. Error", "z value")], 3), sprintf("%.2e", coeff["prs_groupHigh", "Pr(>|z|)"]),
                   
                   round(coeff["prs_groupHigh:environ", c("Estimate", "Std. Error", "z value")], 3), sprintf("%.2e", coeff["prs_groupHigh:environ", "Pr(>|z|)"]),
                   
                   round(F_stat["prs_group", "Chisq"], 3), sprintf("%.2e", F_stat["prs_group", "Pr(>Chisq)"]),
                   round(F_stat["prs_group:environ", "Chisq"], 3), sprintf("%.2e", F_stat["prs_group:environ", "Pr(>Chisq)"]),
                   
                   round(lr_test_stat$Chisq[[2]], 3), sprintf("%.2e", lr_test_stat$`Pr(>Chisq)`[[2]]))
          
          # count number of cases
          total_num = length(unique(data_combine$IID))
          group0_num = length(unique(data_combine$IID[which(data_combine$psu_out==0)]))
          group1_num = length(unique(data_combine$IID[which(data_combine$psu_out==1)]))
          
          out1 = c(prs_thresh, type_race, type_psu, colnames(data_fea)[num_fea], out1, total_num, group0_num, group1_num, 
                   table(data_combine$psu_out, data_combine$prs_group))
          
        }else{
          out1 = rep(NA, dim(out)[2])
        }
        out = rbind(out, out1)
        rm(model, coeff, out1)
      }
    }
  }
}

colnames(out) = c("prs_thresh", "type_race", "type_psu", "envrionment","Beta", "std", "Z-value", "P_value", 
                  "Beta_prs_High", "std_prs_High", "Z-value_prs_High", "P_value_prs_High", 
                  "Beta_prs_High_inter", "std_prs_High_inter", "Z-value_prs_High_inter", "P_value_prs_High_inter", 
                  "anova_Chisq_prs", "anova_Pr(>Chisq)_prs",
                  "anova_Chisq_inter", "anova_Pr(>Chisq)_inter",
                  "LRT_Chisq_inter", "LRT__Pr(>Chisq)_inter", 
                  "Total", "Group0", "Group1", "psu0_low", "psu1_low", "psu0_High", "psu1_High")
out = data.frame(out)

write.csv(out, file = paste0(proj_dir, 'experiment/Table.Step.6.4.prs.environment.interaction.psu.lifetime.two.levels.overall.ses.03.17.2026.csv'))



