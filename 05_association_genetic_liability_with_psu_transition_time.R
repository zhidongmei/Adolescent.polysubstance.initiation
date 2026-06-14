
rm(list=ls())

library(lmerTest)
library(lmtest)
library(coxme)

# ------------------------ 
# plot figure S2 and Table SS correlation among different PGS for substance use
# ------------------------ 
proj_dir = "/local_mount/space/arya/3/users/dmzhi/project/PRS.PSU.environment/"

df = read.csv(paste0(proj_dir, "data/abcd.6.0.pgs.environment.psu.cross.sectional.NA.03.05.2025.csv")) 
df =  df[which(df$analysis_sample_y0_6_5_survive==1),]
df$height_prs2 = df$height_prs;

# ------------------------ 
# Table S_all correlation between different PGS and substance use, OR
# ------------------------ 
for (count in 1){
  
  # year0_3 data 
  analysis_sample = "analysis_sample_y0_6_5_survive";
  phenodata_name_list = c("");
  unrelated_ext = "";
  
  source(paste0(proj_dir, "scripts/util.R"))
  
  # -----------------------------------------
  # load cov data
  # -----------------------------------------
  cov_list = c("Age_ses.00A", "Sex_ses.00A", "FamilyID_ses.00A", "Site_ses.00A", "INR.BL", "Pubertal_mean_ses.00A"); 
  pop_list = c("EUR", "AFR", "AMR");
  
  
  pheno_names = c("NSU2SSU",  "NSU2PSU", "SSU2PSU");  
  prs_names = c("sud", "aud", "cod", "cud", "nid", "oud", 
                "aud_uniq", "cod_uniq", "cud_uniq", "nid_uniq", "oud_uniq" );  
  
  
  for(phenodata_name in phenodata_name_list) {  
    out = c();
    outfile = paste0(proj_dir,  "experiment/psu.survived.time.one.prs.association.", phenodata_name, unrelated_ext, ".r2.out.year0_6_5.03.25.2026.txt"); 
    
    for(pop in pop_list) {   
      
      cat("\n\n\n")
      print(pop);
      
      # for each prs
      for(prs in prs_names) {  
        
        prs_name = paste0(prs, "_prs2");   
        
        # for each pheno 
        for(pheno_name in pheno_names) { 
          
          pheno_name = paste0(pheno_name, phenodata_name);
          
          data = as.data.frame(df[which(df[[analysis_sample]]==1 & df$genetic_ancestry == pop), c("participant_id", prs_name, pheno_name, cov_list, 'psu2_0_6_5_survive', 'SSU2PSU_flag')]); 
          colnames(data)[c(2:10)] = c('SCORE', 'survived_time', 'age', 'sex', 'familyID', 'site', 'INR', 'Tanner_stage', 'flag.transition');
          
          # covariates  - age, sex, site
          data$site = as.factor(data$site);
          data$sex = as.factor(data$sex);
          data$familyID = as.factor(data$familyID);
          
          data$age = as.numeric(scale(data$age));
          data$SCORE = as.numeric(scale(data$SCORE)); 
          data$INR = as.numeric(scale(data$INR)); 
          data$Tanner_stage = as.numeric(scale(data$Tanner_stage)); 
          
          # flag for case
          if (pheno_name == "NSU2SSU"){
            data$flag.transition[which(data$flag.transition==2)] = NA
          }else if (pheno_name == "NSU2PSU"){
            data$flag.transition[which(data$flag.transition==1)] = NA
            data$flag.transition[which(data$flag.transition==2)] = 1
          }else if (pheno_name == "SSU2PSU"){
            data$flag.transition = data$SSU2PSU_flag
          }
          
          data = na.omit(data[, c('survived_time', 'flag.transition', 'SCORE', 'age', 'sex',  'site', 'Tanner_stage', 'familyID', 'INR')])
          
          hr <- coxme(Surv(survived_time, flag.transition) ~ SCORE + age + sex + INR + Tanner_stage + (1|familyID) + (1|site), data = data)
          
          coeff <- summary(hr)$coefficients
          
          out1 = c(pop, analysis_sample, pheno_name, prs_name, round(coeff["SCORE", c("exp(coef)", "se(coef)", "z", "coef")], 3), sprintf("%.2e", coeff["SCORE", "p"]),
                   length(which(!is.na(data$SCORE))), length(which(data$flag.transition==1)), length(which(data$flag.transition==0)))
          
          print(c(prs_name, pheno_name, length(which(!is.na(data$SCORE))), length(which(data$flag.transition==1))))
          
          out = rbind(out, out1);
        }
      } 
    }
    colnames(out)= c("pop", "analysis_sample", "pheno_name", "prs", "HR", "SE", "Z_value", "beta", "P_value", "Participants", "Events", "N_control");
    write.table(out, file=outfile, sep="\t", quote=F, row.names=F, col.names=T, append=F); 
  }
}

