
rm(list=ls());

library(readr)

# ========================================================
# Substance Use Interview
# Release 6.0 Data Table: su_y_sui
# ========================================================
datadir = "C:/Users/DZ1011/Partners HealthCare Dropbox/Dongmei Zhi/abcd6.0.tsv/"

var_list= c("alc", "tob", "can", "oth")
 
gen_su <- function(data, alc, tob, can, oth) {   
  su_y_sui = as.data.frame(data[, c(1:5)]);
  colnames(su_y_sui) = c("IID", "alc", "tob", "can", "oth");
  
  i=2;
  var_list = alc;
  su_y_sui[,i] = 0;
  for(var in var_list) {
    if(var %in% colnames(data)) {  su_y_sui[which(data[, var]==1), i] = 1;  } 
  }
  
  i=3;
  var_list = tob;
  su_y_sui[,i] = 0; 
  for(var in var_list) {
    if(var %in% colnames(data)) {  su_y_sui[ which(data[, var]==1), i] = 1;  } 
  } 
  
  i=4;
  var_list = can;
  su_y_sui[,i] = 0;
  for(var in var_list) {
    if(var %in% colnames(data)) {  su_y_sui[ which(data[, var]==1), i] = 1;  } 
  }
  
  i=5;
  var_list = oth;
  su_y_sui[,i] = 0; 
  for(var in var_list) {
    if(var %in% colnames(data)) {  su_y_sui[ which(data[, var]==1), i] = 1;  } 
  }
  
  for(i in colnames(su_y_sui)[-1]) {
    print(i)
    print(table(su_y_sui[, i]))
  }
  su_y_sui;
} 

gen_eachtime <- function(y0, y1, var_list= c("alc", "tob", "can", "oth"), keep_all=F) {
  
  # original version: we keep everyone
  # now only tox test
  if(keep_all) { 
    dt = as.data.frame(merge(y0, y1, by="IID", all=T));
    
  } else {  
    # if newer data exists, we keep them
    dt = as.data.frame(merge(y0, y1, by="IID", all.y=T));
  }
  
  for(var in colnames(dt)[-1]) {
    dt[which(is.na(dt[,var])), var] = 0;
  }
  
  for(var in var_list) {
    dt = cbind(dt, dt[, paste0(var, ".y")]);
    colnames(dt)[dim(dt)[2]] = var;
    dt[which(dt[, paste0(var, ".x")]==1), var] = 1;
    print(var);
    print(table(dt[, var]));
  } 
  
  dt = dt[, c("IID", var_list)];
  
  if(!keep_all & 1==2) {  
      # if no newer data exists, but they have been cases, we keep them 
      dt1 = as.data.frame(merge(y0, y1, by="IID", all.x=T)); 
      if(length(which(is.na(dt1$alc.y)))>0) {  
        col_no = (dim(dt1)[2]+1)/2;
        dt1 = dt1[which(is.na(dt1$alc.y)), c("IID", paste0(var_list, ".x"))];
        colnames(dt1) = c("IID", var_list);
        
        # only keep when cases
        dt1 = dt1[which(rowSums(dt1[, var_list], na.rm = T)>0),]; 
        
        dt = rbind(dt, dt1);
      }
  }
  
  dt$psu_sum = rowSums(dt[, var_list], na.rm = T);
  dt$psu_sum[rowSums(!is.na(dt[, var_list])) == 0] <- NA
  
  dt$psu2 = dt$psu_sum;
  dt$psu2[which(dt$psu2>=2)] = 2;

  dt$ssu = dt$psu_sum;
  dt$ssu[which(dt$ssu==1)] = 1;
  dt$ssu[which(dt$ssu>1)] = NA;

  dt$psu = dt$psu_sum;
  dt$psu[which(dt$psu==1)] = NA;
  dt$psu[which(dt$psu>1)] = 1;

  dt$asu = dt$psu_sum;
  dt$asu[which(dt$asu>=1)] = 1;
  
  dt;
}

esu_lifetime <- function(y0, y1, var_list= c("alc", "tob", "can", "oth"), keep_all=F) {# lifetime use for each substance use
  
  # original version: we keep everyone
  # now only tox test
  if(keep_all) { 
    dt = as.data.frame(merge(y0, y1, by="IID", all=T));
    
  } else {  
    # if newer data exists, we keep them
    dt = as.data.frame(merge(y0, y1, by="IID", all.y=T));
  }
  
  for(var in var_list) {
    dt = cbind(dt, dt[, paste0(var, ".y")]);
    colnames(dt)[dim(dt)[2]] = var;
    dt[which(dt[, paste0(var, ".x")]==1), var] = 1;
    print(var);
    print(table(dt[, var]));
  } 
  
  dt = dt[, c("IID", var_list)];
  
  if(!keep_all & 1==2) {  
    # if no newer data exists, but they have been cases, we keep them 
    dt1 = as.data.frame(merge(y0, y1, by="IID", all.x=T)); 
    if(length(which(is.na(dt1$alc.y)))>0) {  
      col_no = (dim(dt1)[2]+1)/2;
      dt1 = dt1[which(is.na(dt1$alc.y)), c("IID", paste0(var_list, ".x"))];

      colnames(dt1) = c("IID", var_list);
      
      # only keep when cases
      dt1 = dt1[which(rowSums(dt1[, var_list], na.rm = T)>0),]; 
      
      dt = rbind(dt, dt1);
    }
  }
  
  dt;
}

# ============== multiple substance use at each time point
psu_lifetime <- function(y0, y1, var_list= c("psu_sum"), keep_all=F) {
  
  # original version: we keep everyone
  # now only tox test
  if(keep_all) { 
    dt = as.data.frame(merge(y0, y1, by="IID", all=T));
  } else {  
    # if newer data exists, we keep them
    dt = as.data.frame(merge(y0, y1, by="IID", all.y=T));
  }
  
  for(var in var_list) {
    dt = cbind(dt, dt[, paste0(var, ".y")]);
    colnames(dt)[dim(dt)[2]] = var;
    dt[, var] <- pmax(dt[, paste0(var, ".x")], dt[, paste0(var, ".y")], na.rm = TRUE)
    print(var);
    print(table(dt[, var]));
  } 
  
  dt$psu2 = dt$psu_sum;
  dt$psu2[which(dt$psu2>2)] = 2
  dt$psu2[which(is.na(dt$psu_sum.y) & (dt$psu2<2))] = NA
  
  dt$ssu = dt$psu_sum;
  dt$ssu[which(dt$ssu>1)] = NA
  dt$ssu[which(is.na(dt$psu_sum.y) & (dt$ssu<2))] = NA
  
  dt$psu = dt$psu_sum;
  dt$psu[which(is.na(dt$psu_sum.y) & (dt$psu<2))] = NA
  dt$psu[which(dt$psu==1)] = NA;
  dt$psu[which(dt$psu>1)] = 1;
  
  dt$asu = dt$psu_sum;
  dt$asu[which(dt$psu2>0)] = 1
  dt$asu[which(is.na(dt$psu_sum.y) & (dt$asu==0))] = NA
  
  dt = dt[, c("IID", c("psu_sum", "psu2", "psu", "ssu", "asu"))];
  
  dt;
}

gen_su_sui <- function(data) { 
  alc = c("su_y_sui__use__alc_001");
  alc = c(alc, paste0(alc, "__l"));
  
  tob = c("su_y_sui__use__nic__cig_001" , "su_y_sui__use__nic__cigar_001" , "su_y_sui__use__nic__pipe_001" , "su_y_sui__use__nic__rplc_001",
          "su_y_sui__use__nic__vape_001" ,"su_y_sui__use__nic__hookah_001" ,"su_y_sui__use__nic__chew_001");
  tob = c(tob, paste0(tob, "__l"));
  
  can = c("su_y_sui__use__mj__smoke_001","su_y_sui__use__mj__blunt_001","su_y_sui__use__mj__drink_001", "su_y_sui__use__mj__tinc_001",
          "su_y_sui__use__mj__edbl_001","su_y_sui__use__mj__conc_001", "su_y_sui__use__mj__vape_001__l", "su_y_sui__use__mj__conc__vape_001");
  can = c(can, paste0(can, "__l"));
  
  oth = c("su_y_sui__use__mj__synth_001", 
          "su_y_sui__use__coc_001", 
          "su_y_sui__use__cath_001", 
          "su_y_sui__use__meth_001", 
          "su_y_sui__use__mdma_001", 
          "su_y_sui__use__ket_001", 
          "su_y_sui__use__ghb_001", 
          "su_y_sui__use__opi_001", 
          "su_y_sui__use__hall_001", 
          "su_y_sui__use__shroom_001", 
          "su_y_sui__use__salv_001", 
          "su_y_sui__use__roid_001", 
          "su_y_sui__use__inh__sniff_001",
          # "su_y_sui__use__qc_001", 
          "su_y_sui__use__inh_001", 
          "su_y_sui__use__rxstim_001", 
          "su_y_sui__use__rxsed_001", 
          "su_y_sui__use__dxm_001", 
          "su_y_sui__use__rxopi_001", 
          "su_y_sui__use__othdrg_001", 
          "su_y_sui__use__othdrg_002");
  
  oth = c(oth, paste0(oth, "__l"));
  
  # ---------------------------delete participants with all NA
  can <- setdiff(can, c("su_y_sui__use__mj__conc__vape_001", "su_y_sui__use__mj__conc_001__l", "su_y_sui__use__mj__vape_001__l__l"))
  oth <- setdiff(can, c("su_y_sui__use__othdrg_002", "su_y_sui__use__inh__sniff_001__l"))
  
  name_var = c(alc, tob, can, oth)
  rows_all_na <- which(apply(data[, name_var], 1, function(x) all(x == "n/a")))
  print(paste0("su_y_sui: participants with all NA, n = ", length(rows_all_na)))
  if (length(rows_all_na) > 0){
    data = data[-rows_all_na, ]
  }
  
  rm(rows_all_na, name_var)
  
  gen_su(data, alc, tob, can, oth);
}

name = "su_y_sui";
su_y_sui <- read_tsv(paste0(datadir, "su_y_sui.tsv"))
events = c(0:6);

for(event in events) {
  
  print(paste0("Time ", event))
  
  data <- su_y_sui[which(su_y_sui$session_id==paste0("ses-0", event, "A")), ]
  
  dt = gen_su_sui(data); 
  
  assign(paste0(name, ".year", event), dt)
}

# var_list = c("xskipout_alc_l" , "xskipout_tob_l", "xskipout_mj_l", "xskipout_other_l")
# xskipout_alc_l	Did the participant have at least 2 drinking days in the past 12 months where they consumed one full drink of alcohol or more?	1 = Yes; 0 = No	tlfb_alc_use_l =='1'
# xskipout_tob_l	Did the participant use one or more full cigarette(s) or other tobacco products for at least 2 days in the past 12 months?	1 = Yes; 0 = No	tlfb_cig_use_l == '1' || tlfb_ecig_use_l == '1' || tlfb_hookah_use_l == '1' || tlfb_chew_use_l == '1'
# xskipout_mj_l	Did the participant use any form of marijuana for at least 2 days in the past 12 months.	1 = Yes; 0 = No	tlfb_mj_use_l == '1' || tlfb_blunt_use_l == '1' || tlfb_mj_use_l == '1' || tlfb_edible_use_l == '1' || tlfb_mj_conc_use_l == '1'
# xskipout_other_l	Did the participant use drugs other than alcohol, tobacco, or marijuana for at least 2 days in the past 12 months.	0 = Yes; 0 = No	tlfb_mj_synth_use_l == '1' || tlfb_coc_use_l == '1' || tlfb_bsalts_use_l == '1' || tlfb_meth_use_l == '1' || tlfb_mdma_use_l == '1' || tlfb_ket_use_l == '1' || tlfb_ghb_use_l == '1' || tlfb_opi_use_l == '1' || tlfb_lsd_use_l == '1' || tlfb_shrooms_use_l == '1' || tlfb_salvia_use_l == '1' || tlfb_steroids_use_l == '1' || tlfb_bitta_use_l == '1' || tlfb_inhalant_use_l == '1' || tlfb_amp_use_l == '1' || tlfb_tranq_use_l == '1' || tlfb_cough_use_l == '1' || tlfb_vicodin_use_l == '1' || tlfb_other_use_l == '1'


# ========================================================
# Substance Use Phone Interview (Mid Year)
# Release 6.0 Data Table: su_y_mypi
# ======================================================== 

datadir = "C:/Users/DZ1011/Partners HealthCare Dropbox/Dongmei Zhi/abcd6.0.tsv/"

gen_su_mypi <- function(data) { 
  
  alc = c("su_y_mysu__use__alc__1wk_001", "su_y_mysu__use__alc__1mo_001", "su_y_mysu__use__alc__6mo_001", "su_y_mysu__use__alc__6mo_002");
  tob = c("su_y_mysu__use__nic__1mo_001", "su_y_mysu__use__nic__1wk_001", "su_y_mysu__use__nic__6mo_001", "su_y_mysu__use__nic__6mo_002",
          "su_y_mysu__use__nic__cig__1mo_001", "su_y_mysu__use__nic__cig__6mo_001", "su_y_mysu__use__nic__cigar__6mo_001", "su_y_mysu__use__nic__chew__6mo_002",
          "su_y_mysu__use__nic__chew__1mo_001", "su_y_mysu__use__nic__chew__1wk_001", "su_y_mysu__use__nic__chew__6mo_001",
          "su_y_mysu__use__nic__vape__1mo_001", "su_y_mysu__use__nic__vape__6mo_001");        
  

  can = colnames(su_y_mypi)[grep("su_y_mysu__use__mj__", colnames(su_y_mypi))]
  can = can[grep("su_y_mysu__use__mj__synth__", can, invert=T)]; 
  
  oth = colnames(su_y_mypi)[grep("__use__", colnames(su_y_mypi))]
  oth = oth[grep("__alc_", oth, invert=T)]; 
  oth = oth[grep("__nic_", oth, invert=T)]; 
  oth = oth[grep("__mj_", oth, invert=T)]; 
  oth = oth[grep("__flav_", oth, invert=T)]; 
  oth = oth[grep("__cbd_", oth, invert=T)]; 
  oth = oth[grep("__qc_", oth, invert=T)]; 
  oth = oth[grep("su_y_mysu__use__illdrg_001", oth, invert=T)]; 
  
  oth = oth[grep("participant_id", oth, invert=T)];
  oth = oth[grep("session_id", oth, invert=T)];
  
  # add su_y_mysu__use__mj__synth
  oth = c(oth, colnames(su_y_mypi)[grep("su_y_mysu__use__mj__synth__", colnames(su_y_mypi))])
  
  # ---------------------------delete participants with all NA
  name_var = c(alc, tob, can, oth)
  rows_all_na <- which(apply(data[, name_var], 1, function(x) all(x == "n/a")))
  print(paste0("su_y_mypi: participants with all NA, n = ", length(rows_all_na)))
  if (length(rows_all_na) > 0){
    data = data[-rows_all_na, ]
  }
  
  rm(rows_all_na, name_var)
  
  gen_su(data, alc, tob, can, oth);
}

name = "su_y_mysu";
su_y_mypi <- read_tsv(paste0(datadir,  name, ".tsv"));

events = c(0:5)
for(event in events) {
  
  data <- su_y_mypi[which(su_y_mypi$session_id==paste0("ses-0", event, "M")), ]
  data = data[, -2];
  
  dt = gen_su_mypi(data);
  
  assign(paste0(name, ".month", event), dt)
  
  rm(list=c("dt", "data"));
  
}
#rm(su_y_mypi);


# ===========================================================
# hair toxicology test
# ===========================================================
datadir = "C:/Users/DZ1011/Partners HealthCare Dropbox/Dongmei Zhi/abcd6.0.tsv/"

gen_su_hair_tox <- function(data) {
  var_list = colnames(su_y_hairtox)[grep("_qnt", colnames(su_y_hairtox), invert=T)]
  var_list = var_list[grep("_lab", var_list, invert=T)]; 
  var_list = var_list[grep("__coll_", var_list, invert=T)];
  # var_list = var_list[grep("_scrn", var_list, invert=T)];  # delete Immunoassay screen substance 
  
  alc = var_list[grep("_alc_", var_list)];
  tob = var_list[grep("_nic_", var_list)];
  can = var_list[grep("_mj_", var_list)];
  oth = var_list[!var_list %in% alc];
  oth = oth[!oth %in% tob];
  oth = oth[!oth %in% can];
  
  oth = oth[grep("participant_id", oth, invert=T)];
  oth = oth[grep("session_id", oth, invert=T)];
  
  # ---------------------------delete participants with all NA
  name_var = c(alc, tob, can, oth)
  rows_all_na <- which(apply(data[, name_var], 1, function(x) all(x == "n/a")))
  print(paste0("gen_su_hair: participants with all NA, n = ", length(rows_all_na)))
  if (length(rows_all_na) > 0){
    data = data[-rows_all_na, ]
  }
  
  rm(rows_all_na, name_var)
  
  gen_su(data, alc, tob, can, oth);
}

name = "su_y_hairtox";
su_y_hairtox =  read_tsv(paste0(datadir,  name, ".tsv"));

events = c(0:6);
for(event in events) {
  
  print(event)
  data <- su_y_hairtox[which(su_y_hairtox$session_id==paste0("ses-0", event, "A")), ]
  
  print(var_list)
  dt = gen_su_hair_tox(data); 
  assign(paste0(name, ".year", event), dt)
}


# ===========================================================
# generate co-occurrence at each time point
# ===========================================================
substance_use.y0 = gen_eachtime(su_y_sui.year0, su_y_hairtox.year0, keep_all=T);

substance_use.y0.5 = gen_eachtime(su_y_mysu.month0, su_y_mysu.month0);

substance_use.y1 = gen_eachtime(su_y_sui.year1, su_y_hairtox.year1, keep_all=T);

substance_use.y1.5 = gen_eachtime(su_y_mysu.month1, su_y_mysu.month1);

substance_use.y2 = gen_eachtime(su_y_sui.year2, su_y_hairtox.year2, keep_all=T);

substance_use.y2.5 = gen_eachtime(su_y_mysu.month2, su_y_mysu.month2);

substance_use.y3 = gen_eachtime(su_y_sui.year3, su_y_hairtox.year3, keep_all=T);

substance_use.y3.5 = gen_eachtime(su_y_mysu.month3, su_y_mysu.month3);

substance_use.y4 = gen_eachtime(su_y_sui.year4, su_y_hairtox.year4, keep_all=T);

substance_use.y4.5 = gen_eachtime(su_y_mysu.month4, su_y_mysu.month4);

substance_use.y5 = gen_eachtime(su_y_sui.year5, su_y_hairtox.year5, keep_all=T);

substance_use.y5.5 = gen_eachtime(su_y_mysu.month5, su_y_mysu.month5);

substance_use.y6 = gen_eachtime(su_y_sui.year6, su_y_hairtox.year6, keep_all=T);


# ===========================================================
# generate co-occurrence across life time for each substance use
# ===========================================================
substance_use.y0_0 = substance_use.y0[, c(1:5)]

substance_use.y0_0.5 = esu_lifetime(substance_use.y0[, c(1:5)], substance_use.y0.5[, c(1:5)], keep_all=T);

substance_use.y0_1 = esu_lifetime(substance_use.y0_0.5[, c(1:5)], substance_use.y1[, c(1:5)], keep_all=T);

substance_use.y0_1.5 = esu_lifetime(substance_use.y0_1[, c(1:5)], substance_use.y1.5[, c(1:5)], keep_all=T);

substance_use.y0_2 = esu_lifetime(substance_use.y0_1.5[, c(1:5)], substance_use.y2[, c(1:5)], keep_all=T);

substance_use.y0_2.5 = esu_lifetime(substance_use.y0_2[, c(1:5)], substance_use.y2.5[, c(1:5)], keep_all=T);

substance_use.y0_3 = esu_lifetime(substance_use.y0_2.5[, c(1:5)], substance_use.y3[, c(1:5)], keep_all=T);

substance_use.y0_3.5 = esu_lifetime(substance_use.y0_3[, c(1:5)], substance_use.y3.5[, c(1:5)], keep_all=T);

substance_use.y0_4 = esu_lifetime(substance_use.y0_3.5[, c(1:5)], substance_use.y4[, c(1:5)], keep_all=T);

substance_use.y0_4.5 = esu_lifetime(substance_use.y0_4[, c(1:5)], substance_use.y4.5[, c(1:5)], keep_all=T);

substance_use.y0_5 = esu_lifetime(substance_use.y0_4.5[, c(1:5)], substance_use.y5[, c(1:5)], keep_all=T);

substance_use.y0_5.5 = esu_lifetime(substance_use.y0_5[, c(1:5)], substance_use.y5.5[, c(1:5)], keep_all=T);

substance_use.y0_6 = esu_lifetime(substance_use.y0_5.5[, c(1:5)], substance_use.y6[, c(1:5)], keep_all=T);

# ===========================================================
# generate co-occurrence across life time for polysubstance use
# ===========================================================
substance_use_psu.y0_0 = substance_use.y0[, c(1, 6:10)]

substance_use_psu.y0_0.5 = psu_lifetime(substance_use.y0[, c(1, 6:10)], substance_use.y0.5[, c(1, 6:10)], keep_all=T);

substance_use_psu.y0_1 = psu_lifetime(substance_use_psu.y0_0.5, substance_use.y1[, c(1, 6:10)], keep_all=T);

substance_use_psu.y0_1.5 = psu_lifetime(substance_use_psu.y0_1, substance_use.y1.5[, c(1, 6:10)], keep_all=T);

substance_use_psu.y0_2 = psu_lifetime(substance_use_psu.y0_1.5, substance_use.y2[, c(1, 6:10)], keep_all=T);

substance_use_psu.y0_2.5 = psu_lifetime(substance_use_psu.y0_2, substance_use.y2.5[, c(1, 6:10)], keep_all=T);

substance_use_psu.y0_3 = psu_lifetime(substance_use_psu.y0_2.5, substance_use.y3[, c(1, 6:10)], keep_all=T);

substance_use_psu.y0_3.5 = psu_lifetime(substance_use_psu.y0_3, substance_use.y3.5[, c(1, 6:10)], keep_all=T);

substance_use_psu.y0_4 = psu_lifetime(substance_use_psu.y0_3.5, substance_use.y4[, c(1, 6:10)], keep_all=T);

substance_use_psu.y0_4.5 = psu_lifetime(substance_use_psu.y0_4, substance_use.y4.5[, c(1, 6:10)], keep_all=T);

substance_use_psu.y0_5 = psu_lifetime(substance_use_psu.y0_4.5, substance_use.y5[, c(1, 6:10)], keep_all=T);

substance_use_psu.y0_5.5 = psu_lifetime(substance_use_psu.y0_5, substance_use.y5.5[, c(1, 6:10)], keep_all=T);

substance_use_psu.y0_6 = psu_lifetime(substance_use_psu.y0_5.5, substance_use.y6[, c(1, 6:10)], keep_all=T);


for(event in c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6)) {
  
  dt = get(paste0("substance_use.y0_", event));
  
  print(event);
  
  for(var in var_list) {  
    # print( table(dt[, var])  )
    print(paste(var, sprintf("%.2f", mean(dt[, var], na.rm=T)*100)))
  }
  cat("\n\n") 
}

for(event in c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6)) {
  if (event == 0){
    dt = substance_use.y0
  }else{
    dt = get(paste0("substance_use_psu.y0_", event));
  }
  
  print(event);
  print(table(dt[, "psu2"])) 
  psu_no = length(which(dt[,"psu2"]>1));
  ssu_no = length(which(dt[,"psu2"]==1));
  nsu_no = length(which(dt[,"psu2"]==0));
  
  all_no = psu_no + ssu_no + nsu_no;
  
  print(table(dt[, "psu2"])/all_no) 
  
  print(paste0(sprintf("%.2f", (psu_no*100)/all_no), "(%)")) 
  
  cat("\n\n")
}
rm(list=c( "dt"));


# renamed original version with everyone, even when no newer data exists, to psu.data.v0.Rdb
if(1==1) {  
outdir = "C:/Users/DZ1011/Dropbox/2024.PSU.lancet/data/psu/" 
outfile = paste0(outdir, "abcd.6.0.substance_use.definition.NA.03.05.2026.Rdb");
target_list = c(ls()[grep ("su_y_hairtox.year", ls())], ls()[grep ("su_y_mysu.month", ls())], ls()[grep ("su_y_sui.year", ls())]);
save(list=target_list, file=outfile)

# ------------------lifetime definition
merge_variable = function(name_var, event_time){
  
  df_list <- list()
  for (event in event_time) {
    
    df <- get(paste0(name_var, event))
    if (name_var == "substance_use.y"){
      colnames(df)[2:dim(df)[2]] <- paste0(colnames(df)[2:dim(df)[2]], "_", event)
    }else if (name_var == "substance_use.y0_"){
      colnames(df)[2:dim(df)[2]] <- paste0(colnames(df)[2:dim(df)[2]], "_0_", event)
    }else if (name_var == "substance_use_psu.y0_"){
      colnames(df)[2:dim(df)[2]] <- paste0(colnames(df)[2:dim(df)[2]], "_0_", event)
    }
    
    # save to df list 
    df_list[[as.character(event)]] <- df
  }
  
  # use Reduce to merge all variable by IID
  df_merged <- Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), df_list)
}

substance_use.y_all = merge_variable("substance_use.y", c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6))
substance_use.y0_all = merge_variable("substance_use.y0_", c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6))
substance_use_psu.y0_all = merge_variable("substance_use_psu.y0_", c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6))

substance_use_lifetime = Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), list(substance_use.y_all, substance_use.y0_all, substance_use_psu.y0_all))

# define all time point at 6.0 
for (type_sub in c("alc", "tob", "can", "oth", "alc", "psu2")){
  substance_use_lifetime[, paste0(type_sub, "_0_6_5")] = substance_use_lifetime[, paste0(type_sub, "_0_6")]
  substance_use_lifetime[which(is.na(substance_use_lifetime[, paste0(type_sub, "_0_6_5")])), paste0(type_sub, "_0_6_5")] = substance_use_lifetime[which(is.na(substance_use_lifetime[, paste0(type_sub, "_0_6_5")])), paste0(type_sub, "_0_5.5")]
  substance_use_lifetime[which(is.na(substance_use_lifetime[, paste0(type_sub, "_0_6_5")])), paste0(type_sub, "_0_6_5")] = substance_use_lifetime[which(is.na(substance_use_lifetime[, paste0(type_sub, "_0_6_5")])), paste0(type_sub, "_0_5")]
}

substance_use_lifetime$ssu_0_6_5 = substance_use_lifetime$psu2_0_6_5;
substance_use_lifetime$ssu_0_6_5[which(substance_use_lifetime$ssu_0_6_5==1)] = 1;
substance_use_lifetime$ssu_0_6_5[which(substance_use_lifetime$ssu_0_6_5>1)] = NA;

substance_use_lifetime$psu_0_6_5 = substance_use_lifetime$psu2_0_6_5;
substance_use_lifetime$psu_0_6_5[which(substance_use_lifetime$psu_0_6_5==1)] = NA;
substance_use_lifetime$psu_0_6_5[which(substance_use_lifetime$psu_0_6_5>1)] = 1;

substance_use_lifetime$asu_0_6_5 = substance_use_lifetime$psu2_0_6_5;
substance_use_lifetime$asu_0_6_5[which(substance_use_lifetime$asu_0_6_5>=1)] = 1;

outfile = paste0(outdir, "abcd.6.0.substance_use.lifetime.NA.03.05.2026.Rdb");
save(substance_use_lifetime, file=outfile) 
}

