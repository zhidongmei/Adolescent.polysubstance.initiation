
rm(list=ls());

library(glmmTMB)
library(dplyr)

proj_dir = "/local_mount/space/arya/3/users/dmzhi/project/PRS.PSU.environment/"

# ========================================================
# psu and environment associations
# ========================================================
data_fea = read.csv(paste0(proj_dir, "data/abcd.6.0.pgs.environment.psu.cross.sectional.NA.03.05.2025.csv")) 
data_fea = data_fea[which(data_fea$analysis_sample_y0_6_5_all==1), ]

#------------------------ logistic regression with psu as outcome
for (type_psu in c("psu_0_6_5", "ssu_0_6_5", "psu_ssu_0_6_5")){
  
  if (type_psu=="psu_ssu_0_6_5"){
    data_fea[, "psu_ssu_0_6_5"] = data_fea$psu_0_6_5
    data_fea$psu_ssu_0_6_5 = ifelse(is.na(data_fea$psu_ssu_0_6_5),  0,  
                                    ifelse(data_fea$psu_ssu_0_6_5 == 0, NA, data_fea$psu_ssu_0_6_5)) 
  }
  
  out = c()
  num_environ = c(2:51, 53:90)
  
  for (num_fea in num_environ){    
    
    print(num_fea)
    
    data_combine = data.frame(IID = data_fea$participant_id, environ = data_fea[, num_fea], psu_out = data_fea[, type_psu], 
                              age = data_fea$Age_ses.00A, sex = data_fea$Sex_ses.00A, race = data_fea$genetic_ancestry,
                              site = data_fea$Site_ses.00A, familyID = data_fea$FamilyID_ses.00A, 
                              Tanner_stage = data_fea$Pubertal_mean_ses.00A, INR = data_fea$INR.BL)
    
    data_combine = na.omit(data_combine[, c('IID', 'environ', 'psu_out', 'age', 'sex', 'Tanner_stage', 'INR', 'race', 'site', 'familyID')])
    
    data_combine$psu_out <- as.numeric(data_combine$psu_out)
    data_combine$sex <- as.factor(data_combine$sex)
    data_combine$race <- as.factor(data_combine$race)
    data_combine$site <- as.factor(data_combine$site)
    data_combine$familyID <- as.factor(data_combine$familyID)
    
    data_combine$environ <- as.numeric(scale(data_combine$environ))
    data_combine$age <- as.numeric(scale(data_combine$age))
    data_combine$Tanner_stage <- as.numeric(scale(data_combine$Tanner_stage))
    data_combine$INR <- scale(as.numeric(data_combine$INR))
    
    if (dim(data_combine)[1]>500){
      
      model <- glmmTMB(psu_out ~ environ + age + sex + Tanner_stage + INR + race + (1 | familyID)+ (1 | site), data = data_combine, family = binomial(link = "logit"))
      
      coeff = summary(model)$coefficients$cond
      
      # extract each term, and interaction term
      out1 = c(round(coeff["environ", c("Estimate", "Std. Error", "z value")], 3), sprintf("%.2e", coeff["environ", "Pr(>|z|)"]))
      
      # count number of cases
      total_num = length(unique(data_combine$IID))
      group0_num = length(unique(data_combine$IID[which(data_combine$psu_out==0)]))
      group1_num = length(unique(data_combine$IID[which(data_combine$psu_out==1)]))
      
      out1 = c(out1, total_num, group0_num, group1_num)
      
    }else{
      out1 = NA
    }
    out = rbind(out, out1)
    rm(model, coeff, out1)
  }
  
  colnames(out) = c("Beta", "std", "Z-value", "P_value", "Total", "Group0", "Group1")
  
  out = data.frame(out)
  
  out$P_FDR = p.adjust(as.numeric(out$P_value), method = "fdr")
  
  rownames(out) = colnames(data_fea)[num_environ]
  
  write.csv(out, file = paste0(proj_dir, 'experiment/Table.Step.5.psu.environment.association.', type_psu, '.03.17.2026.csv',sep=""), row.names = TRUE)
}


