
rm(list=ls());

library(readr)
library(coxme)
library(dplyr)

proj_dir = "/local_mount/space/arya/3/users/dmzhi/project/PRS.PSU.environment/"
datadir = "/local_mount/space/arya/3/users/dmzhi/project/public.data/abcd6.0.tsv/"

# ========================================================
# load basic information
# ========================================================
data_fea_all = read.csv(paste0(proj_dir, "data/abcd.6.0.pgs.environment.psu.cross.sectional.NA.03.05.2025.csv")) 
data_fea_all = data_fea_all[which(data_fea_all$analysis_sample_y0_6_5_survive==1),]

# ========================================================
# PRS imaging association
# ========================================================
for (time_point in c("ses-00A", "ses-02A", "ses-04A", "ses-06A")){ # , "ses-04A", "ses-06A"
  
  time_point0 = gsub("-", ".", time_point)
  
  for (img_file in c("mr_y_rsfmri__corr__gpnet", "mr_y_rsfmri__corr__gpnet__aseg")){#"mr_y_rsfmri__corr__gpnet", "mr_y_rsfmri__corr__gpnet__aseg"
    
    img_fea =  read_tsv(paste0(datadir, img_file, ".tsv"))
    
    data_fea <- data_fea_all %>%
      left_join(img_fea %>% filter(session_id == time_point), by = c("participant_id"))
    
    #------------------------ logistic regression with psu as outcome
    for (type_psu in c("SSU2PSU", "NSU2SSU", "NSU2PSU")){ 
      # ========================================================
      # quality control
      # ========================================================
      if (grepl("mr_y_rsfmri__corr", img_file)){
        data_fea$motion_ICV = data_fea[, paste0('mr_y_qc__mot__rsfmri__mot_mean_', time_point0)]
        sel_sub = which(data_fea[, paste0('mr_y_qc__incl__rsfmri_indicator_', time_point0)]==1 & data_fea$motion_ICV<0.2)
      }
      
      data_fea = data_fea[sel_sub, ]
      
      # ========================================================
      # number of image feature
      # ========================================================
      num_image = c(790:(dim(data_fea)[2]-1))  
      
      out = c()
      for (num_fea in num_image){    
        
        print(paste0(num_fea, colnames(data_fea)[num_fea]))
        
        data_fea[which(data_fea[, num_fea] == "n/a"), num_fea] <- NA 
        
        if (type_psu != "SSU2PSU"){
          data_combine = data.frame(IID = data_fea$participant_id, image_fea = data_fea[, num_fea], survived_time = data_fea[, type_psu], 
                                    flag.transition = data_fea$psu2_0_6_5_survive, 
                                    age = data_fea[, paste0('Age_', time_point0)], sex = data_fea$Sex_ses.00A, race = data_fea$genetic_ancestry,
                                    site=data_fea[, paste0('Site_', time_point0)], familyID = data_fea$FamilyID_ses.00A, 
                                    Tanner_stage = data_fea[, paste0('Pubertal_mean_', time_point0)], INR = data_fea$INR.BL,
                                    scanner = data_fea[, paste0('Scanner_', time_point0)], motion = data_fea$motion_ICV)
          
        } else {
          data_combine = data.frame(IID = data_fea$participant_id, image_fea = data_fea[, num_fea], survived_time = data_fea[, type_psu], 
                                    flag.transition = data_fea$SSU2PSU_flag, 
                                    age = data_fea[, paste0('Age_', time_point0)], sex = data_fea$Sex_ses.00A, race = data_fea$genetic_ancestry,
                                    site=data_fea[, paste0('Site_', time_point0)], familyID = data_fea$FamilyID_ses.00A, 
                                    Tanner_stage = data_fea[, paste0('Pubertal_mean_', time_point0)], INR = data_fea$INR.BL,
                                    scanner = data_fea[, paste0('Scanner_', time_point0)], motion = data_fea$motion_ICV)
        }
        
        data_combine = na.omit(data_combine[, c('IID', 'survived_time', 'image_fea', 'motion', 'age', 'sex', 'race', 'scanner', 'site', 
                                                'familyID', 'INR', 'Tanner_stage', 'flag.transition')])
        
        data_combine$survived_time <- as.numeric(data_combine$survived_time)
        data_combine$sex <- as.factor(data_combine$sex)
        data_combine$race <- as.factor(data_combine$race)
        data_combine$site <- as.factor(data_combine$site)
        data_combine$familyID <- as.factor(data_combine$familyID)
        data_combine$scanner <- as.factor(data_combine$scanner)
        
        data_combine$image_fea <- scale(as.numeric(data_combine$image_fea))
        data_combine$age <- as.numeric(scale(data_combine$age))
        data_combine$Tanner_stage <- as.numeric(scale(data_combine$Tanner_stage))
        data_combine$motion <- scale(as.numeric(data_combine$motion))
        data_combine$INR <- scale(as.numeric(data_combine$INR))
        
        # convert the flag.transition
        if (type_psu == "NSU2SSU"){
          data_combine$flag.transition[which(data_combine$flag.transition==2)] = NA
        }else if(type_psu == "NSU2PSU"){
          data_combine$flag.transition[which(data_combine$flag.transition==1)] = NA
          data_combine$flag.transition[which(data_combine$flag.transition==2)] = 1
        }
        
        if (dim(data_combine)[1]>500){
          
          hr <- coxme(Surv(survived_time, flag.transition) ~ image_fea + age + scanner + sex + motion + INR +  Tanner_stage + race + 
                        (1 | familyID) + (1 | site), data = data_combine)
          
          coeff = summary(hr)$coefficients
          
          # hr_reduced = coxme(Surv(survived_time, flag.transition) ~ age + sex + scanner + motion + Tanner_stage + race + (1 | familyID) + (1 | site), data = data_combine)
          
          # extract each term, and interaction term 
          out1 = c(round(coeff["image_fea", c("exp(coef)", "se(coef)", "z", "coef")], 3), sprintf("%.2e", coeff["image_fea", "p"]))
          
          # count number of cases
          total_num = length(unique(data_combine$IID))
          group0_num = length(unique(data_combine$IID[which(data_combine$flag.transition==0)]))
          group1_num = length(unique(data_combine$IID[which(data_combine$flag.transition==1)]))
          
          out1 = c(out1, total_num, group0_num, group1_num)
          
        }else{
          out1 = NA
        }
        out = rbind(out, out1)
        rm(model, coeff, out1)
      }
      
      colnames(out) = c("HR", "SE", "Z-value", "beta", "P_value", "Total", "Group0", "Group1")
      
      rownames(out) = colnames(data_fea)[num_image]
      out = data.frame(out)
      out$P_FDR = p.adjust(as.numeric(out$P_value), method = "fdr")
      
      output_line = paste0(img_file, '_', type_psu, ' ', time_point, ' ', total_num, ' FDR significant association:   ', length(which(out$P_FDR<0.05)))
      output_line2 = paste0(img_file, '_', type_psu, ' ', time_point, ' ', total_num, ' not corrected association:   ', length(which(as.numeric(out$P_value)<0.05)))
      
      print(output_line)
      print(output_line2)
      write(output_line, file = paste0(proj_dir, "experiment/psu.image.y0_5.txt"), append = TRUE)
      write(output_line2, file = paste0(proj_dir, "experiment/psu.image.y0_5.txt"), append = TRUE)
      
      write.csv(out, file = paste0(proj_dir, 'experiment/Table.Step.8.3.psu.survival.time.imaging.association.', type_psu, '.', time_point, '.', img_file, '.meanFD_0.2mm.03.25.2026.csv',sep=""), row.names = TRUE)
    }
  }
}





