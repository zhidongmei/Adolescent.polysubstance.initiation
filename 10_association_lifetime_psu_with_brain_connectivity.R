
rm(list=ls());

library(readr)
library(glmmTMB)
library(dplyr)

proj_dir = "/local_mount/space/arya/3/users/dmzhi/project/PRS.PSU.environment/"
datadir = "/local_mount/space/arya/3/users/dmzhi/project/public.data/abcd6.0.tsv/"

# ========================================================
# load basic information
# ========================================================
data_fea_all = read.csv(paste0(proj_dir, "data/abcd.6.0.pgs.environment.psu.cross.sectional.NA.03.05.2025.csv")) 
data_fea_all = data_fea_all[which(data_fea_all$analysis_sample_y0_6_5_all==1),]
data_fea_all[, "psu_ssu_0_6_5"] = data_fea_all$psu_0_6_5
data_fea_all$psu_ssu_0_6_5 = ifelse(is.na(data_fea_all$psu_ssu_0_6_5),  0,  
                             ifelse(data_fea_all$psu_ssu_0_6_5 == 0, NA, data_fea_all$psu_ssu_0_6_5))

# ========================================================
# PRS imaging association
# ========================================================
for (time_point in c("ses-00A", "ses-02A", "ses-04A", "ses-06A")){
  
  time_point0 = gsub("-", ".", time_point)
  
  for (img_file in c("mr_y_rsfmri__corr__gpnet", "mr_y_rsfmri__corr__gpnet__aseg")){#
    
    img_fea =  read_tsv(paste0(datadir, img_file, ".tsv"))
    
    data_fea <- data_fea_all %>%
      left_join(img_fea %>% filter(session_id == time_point), by = c("participant_id"))
    
    #------------------------ logistic regression with psu as outcome
    for (type_psu in c("psu_ssu_0_6_5")){ #"psu_0_6_5", "ssu_0_6_5", "psu_ssu_0_6_5", "psu_sum_0_6_5"
      
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
      if (type_psu=="psu_ssu_0_6_5"){
        num_image = c(791:(dim(data_fea)[2]-1))
      }else{
        num_image = c(791:(dim(data_fea)[2]-1))
      }
      
      out = c()
      for (num_fea in num_image){    
        
        print(num_fea)
        
        data_fea[which(data_fea[, num_fea] == "n/a"), num_fea] <- NA 
        
        data_combine = data.frame(IID = data_fea$participant_id, psu_out = data_fea[, type_psu], image_fea = data_fea[, num_fea], 
                                  age = data_fea[, paste0('Age_', time_point0)], sex = data_fea$Sex_ses.00A, race = data_fea$genetic_ancestry,
                                  site=data_fea[, paste0('Site_', time_point0)], familyID = data_fea$FamilyID_ses.00A, 
                                  Tanner_stage = data_fea[, paste0('Pubertal_mean_', time_point0)], INR = data_fea$INR.BL, 
                                  scanner = data_fea[, paste0('Scanner_', time_point0)], motion = data_fea$motion_ICV) 
        
        data_combine = na.omit(data_combine[, c('IID', 'psu_out', 'image_fea', 'motion', 'age', 'sex', 'Tanner_stage', 'INR', 'race', 'site', 'scanner', 'familyID')])
        
        data_combine$psu_out <- as.numeric(data_combine$psu_out)
        data_combine$sex <- as.factor(data_combine$sex)
        data_combine$race <- as.factor(data_combine$race)
        data_combine$site <- as.factor(data_combine$site)
        data_combine$familyID <- as.factor(data_combine$familyID)
        data_combine$scanner <- as.factor(data_combine$scanner)
        
        data_combine$image_fea <- scale(as.numeric(data_combine$image_fea))
        data_combine$age <- as.numeric(scale(data_combine$age))
        data_combine$Tanner_stage <- as.numeric(scale(data_combine$Tanner_stage))
        data_combine$motion <- scale(as.numeric(data_combine$motion))
        data_combine$INR <- as.numeric(scale(data_combine$INR))
        
        if (dim(data_combine)[1]>500){
          
          model <- glmmTMB(psu_out ~ image_fea + age + sex + motion + Tanner_stage + INR + race + (1 | scanner) + (1 | familyID) + (1 | site), data = data_combine, family = binomial(link = "logit"))
          
          coeff = summary(model)$coefficients$cond
          
          # extract each term, and interaction term 
          out1 = c(round(coeff["image_fea", c("Estimate", "Std. Error", "z value")], 3), sprintf("%.2e", coeff["image_fea", "Pr(>|z|)"]))
          
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
      
      rownames(out) = colnames(data_fea)[num_image]
      out = data.frame(out)
      out$P_FDR = p.adjust(as.numeric(out$P_value), method = "fdr")
      
      output_line = paste0(img_file, '_', type_psu, ' ', time_point, ' ', total_num, ' FDR significant association:   ', length(which(out$P_FDR<0.05)))
      output_line2 = paste0(img_file, '_', type_psu, ' ', time_point, ' ', total_num, ' not corrected association:   ', length(which(as.numeric(out$P_value)<0.05)))
      
      print(output_line)
      print(output_line2)
      write(output_line, file = paste0(proj_dir, "experiment/psu.image.y0_5.txt"), append = TRUE)
      write(output_line2, file = paste0(proj_dir, "experiment/psu.image.y0_5.txt"), append = TRUE)
      
      write.csv(out, file = paste0(proj_dir, 'experiment/Table.8.1.psu.imaging.association.', type_psu, '.', time_point, '.', img_file, '.meanFD_0.2mm.03.17.2026.csv',sep=""), row.names = TRUE)
    }
  }
}







