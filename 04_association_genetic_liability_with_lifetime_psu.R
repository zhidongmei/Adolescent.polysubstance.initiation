
rm(list=ls())

library(lmerTest)
library(glmmTMB)
library(performance)
library(lmtest)

# ------------------------ 
# plot figure S2 and Table SS correlation among different PGS for substance use
# ------------------------ 
proj_dir = "/local_mount/space/arya/3/users/dmzhi/project/PRS.PSU.environment/"

df = read.csv(paste0(proj_dir, "data/abcd.6.0.pgs.environment.psu.cross.sectional.NA.03.05.2025.csv")) 
df =  df[which(df$analysis_sample_y0_6_5_all==1),]
df$height_prs2 = df$height_prs;

# ------------------------ 
# Table S_all correlation between different PGS and substance use, OR
# ------------------------ 
for (count in 1){
  
  # year0_6_5 data 
  analysis_sample = "analysis_sample_y0_6_5_all";
  phenodata_name_list = c("_0_6_5");
  
  
  unrelated_ext = "";
  
  source(paste0(proj_dir, "scripts/util.R"))
  
  # -----------------------------------------
  # load cov data
  # -----------------------------------------
  cov_list = c("Age_ses.00A", "Sex_ses.00A", "FamilyID_ses.00A", "Site_ses.00A", "INR.BL", "Pubertal_mean_ses.00A"); 
  pop_list = c("EUR", "AFR", "AMR" );
  
  
  pheno_names = c("psu",  "asu", "ssu", "tob", "alc", "can", "oth");  
  prs_names = c("sud", "aud", "cod", "cud", "nid", "oud", 
                "aud_uniq", "cod_uniq", "cud_uniq", "nid_uniq", "oud_uniq" );  
  
  for(phenodata_name in phenodata_name_list) {  
    out = c();
    outfile = paste0(proj_dir,  "experiment/psu.one.prs.lifetime.association.", phenodata_name, unrelated_ext, ".r2.out.year0_6_5.03.17.2026.txt"); 
    
    for(pop in pop_list) {   
      
      cat("\n\n\n")
      print(pop);
      
      # for each prs
      for(prs in prs_names) {  
        
        prs_name = paste0(prs, "_prs2");   
        
        # for each pheno 
        for(pheno_name in pheno_names) { 
          
          pheno_name = paste0(pheno_name, phenodata_name);
          
          data = as.data.frame(df[which(df[[analysis_sample]]==1 & df$genetic_ancestry == pop), c("participant_id", prs_name, pheno_name, cov_list)]); 
          colnames(data)[c(2, 4:9)] = c('SCORE', 'age', 'sex', 'familyID', 'site', 'INR', 'Tanner_stage')
          
          # covariate  - age, sex, site
          data$site = as.factor(data$site);
          data$sex = as.factor(data$sex);
          data$familyID = as.factor(data$familyID);
          
          data$age = as.numeric(scale(data$age));
          data$SCORE = as.numeric(scale(data$SCORE)); 
          data$INR = as.numeric(scale(data$INR)); 
          data$Tanner_stage = as.numeric(scale(data$Tanner_stage)); 
          
          data = na.omit(data[, c(pheno_name, 'SCORE', 'age', 'sex', 'INR', 'Tanner_stage', 'familyID', 'site')])
          
          # pheno - binary
          family_type = "binomial";
          
          # mixed effects model 
          f0 = as.formula(paste0(pheno_name, "~ age + sex + INR + Tanner_stage + (1|familyID) + (1|site)"));
          f1 = as.formula(paste0(pheno_name, "~ age + sex + INR + Tanner_stage + SCORE + (1|familyID) + (1|site)"));
          
          
          # run regression
          if(family_type=="gaussian") {   
            model_null = lmer(f0, data)  
            model_full = lmer(f1, data) 
          }
          if(family_type=="binomial") {   
            model_null = glmmTMB(f0, data, family=family_type)  
            model_full = glmmTMB(f1, data, family=family_type)
          }
          
          # likelihood ratio test 
          lr_test_stat = lrtest(model_null, model_full)
          lr_test_p = sprintf("%.2e", lr_test_stat$`Pr(>Chisq)`[[2]])
          
          # r2 - figure out glmer r2 calculation
          r2_full = performance::r2(model_full)
          r2_null = performance::r2(model_null)
          
          if(family_type=="binomial") {  
            r2_full_marginal = r2_full$R2_marginal
            r2_null_marginal = r2_null$R2_marginal
          } else {
            r2_full_marginal = r2_full[1,1];
            r2_null_marginal = r2_null[1,1]
          }
          
          r2_diff = sprintf("%.2f", 100 * (r2_full_marginal - r2_null_marginal))
          r2_pct = sprintf("%.2f", ((r2_full_marginal - r2_null_marginal) / r2_full_marginal) * 100)
          
          coeff = summary(model_full)$coefficients$cond
          
          odds_ratio <- exp(as.numeric(coeff["SCORE", "Estimate"]))
          lower_CI <- exp(as.numeric(coeff["SCORE", "Estimate"]) - 1.96 * as.numeric(coeff["SCORE", "Std. Error"]))
          upper_CI <- exp(as.numeric(coeff["SCORE", "Estimate"]) + 1.96 * as.numeric(coeff["SCORE", "Std. Error"]))
          
          out1 = cbind(pop, analysis_sample, pheno_name, prs_name, "SCORE", round(coeff["SCORE", c("Estimate")], 3), round(coeff["SCORE", c("Std. Error")], 3), 
                       round(odds_ratio, 3),  round(lower_CI, 3), round(upper_CI, 3), sprintf("%.2e", coeff["SCORE", "Pr(>|z|)"]), 
                       length(which(data[, pheno_name]==1)), length(which(data[, pheno_name]==0)), length(which(!is.na(data[, pheno_name]))), 
                       rep(lr_test_p), rep(r2_diff), rep(r2_pct));
          
          out = rbind(out, out1);
          
          print(out1);
        }
      } 
    }
    colnames(out)= c("pop", "analysis_sample", "pheno_name", "prs", "predictor", "beta", "se", "or","ci_lower", "ci_upper", "p.value", "case_no", "ctrl_no", "all_no", "lr_test_p", "r2_diff_in_pct", "r2_predictor_model_pct");
    write.table(out, 
                file=outfile, 
                sep="\t", quote=F, row.names=F, col.names=T, append=F); # file.exists(outfile) 
  }
}





