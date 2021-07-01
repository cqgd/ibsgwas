library(data.table)
setwd("~/work/ibs")

save_zip = function(object,path){
  fwrite(object,path)
  zip(paste0(path,".zip"),path)
}

grep_start = function(col,code){
  grepl(paste0("^",code),col)
}

check_code = function(code,P,cols,col_prefix="has",starts_with=FALSE){
  if((code=="") | is.na(code)){return(NULL)}
  
  if(starts_with){
    noted = P[,lapply(.SD,grep_start,code),.SDcols=cols]
  }else{
    noted = P[,lapply(.SD,"==",code),.SDcols=cols]
  }
  
  noted = rowSums(noted, na.rm=TRUE)>=1
  
  res_col = paste0(col_prefix,"_",code)
  set(P,i=NULL,j=res_col,noted)
  return(res_col)
}

P = fread("pheno/ukb21760.csv")

### COLS ###
#only here for RAM
diag_cols = grep("41202|41204",colnames(P),value = TRUE)
sr_cols = grep("20002",colnames(P),value = TRUE)

P = P[,.SD,.SDcols = c("eid",diag_cols,sr_cols,"20127-0.0")]
### end COLS ###

### CODING ###
X = fread("resource/IBS GWAS_ SR and ICD-10 codes of interest - diagnoses_selfreported.csv")
### end CODING ###

######### DIAGNOSED #########

#just include & ICD10.code!="" when passing the codes to avoid having to unlist and to combine the new cols at the end.
diag.exclusion_cols = sapply(X[category=="exclusion",unique(ICD10.code)], check_code, P=P, cols=diag_cols, col_prefix="diag.exclusion", starts_with=TRUE)
P[, lapply(.SD, sum), .SDcols=unlist(diag.exclusion_cols)]

diag.surrogate_cols = sapply(X[category=="surrogate",unique(ICD10.code)], check_code, P=P, cols=diag_cols, col_prefix="diag.surrogate", starts_with=TRUE)
P[, lapply(.SD, sum), .SDcols=unlist(diag.surrogate_cols)]

diag.IBS_cols = check_code(X[category=="IBS",ICD10.code], P=P, cols=diag_cols, col_prefix="diag.IBS", starts_with=TRUE)
P[, lapply(.SD, sum), .SDcols=unlist(diag.IBS_cols)]

diag.interest_cols = sapply(X[category=="interest",unique(ICD10.code)], check_code, P=P, cols=diag_cols, col_prefix="diag.interest", starts_with=TRUE) #added for clinical table, was not sapply at first -- all results the same?
P[, lapply(.SD, sum), .SDcols=unlist(diag.interest_cols)]

diag.corr_cols = sapply(X[category=="correlation plot",unique(ICD10.code)], check_code, P=P, cols=diag_cols, col_prefix="diag.correlation", starts_with=TRUE) #added for clinical table
P[, lapply(.SD, sum), .SDcols=unlist(diag.corr_cols)]
######### end DIAGNOSED #########



######### SELF-REPORTED #########

sr.exclusion_cols = sapply(X[category=="exclusion",unique(SR.code)], check_code, P=P, cols=sr_cols, col_prefix="sr.exclusion", starts_with=FALSE)
P[, lapply(.SD, sum), .SDcols=unlist(sr.exclusion_cols)]

sr.surrogate_cols = sapply(X[category=="surrogate",unique(SR.code)], check_code, P=P, cols=sr_cols, col_prefix="sr.surrogate", starts_with=FALSE)
P[, lapply(.SD, sum), .SDcols=unlist(sr.surrogate_cols)]

sr.IBS_cols = check_code(X[category=="IBS",unique(SR.code)], P=P, cols=sr_cols, col_prefix="sr.IBS")
P[, lapply(.SD, sum), .SDcols=unlist(sr.IBS_cols)]

sr.interest_cols = sapply(X[category=="interest",unique(SR.code)], check_code, P=P, cols=sr_cols, col_prefix="sr.interest") #added for clinical table
P[, lapply(.SD, sum), .SDcols=unlist(sr.interest_cols)]

sr.corr_cols = sapply(X[category=="correlation plot",unique(SR.code)], check_code, P=P, cols=sr_cols, col_prefix="sr.correlation", starts_with=TRUE) #added for clinical table
P[, lapply(.SD, sum), .SDcols=unlist(sr.corr_cols)]

#############

#Neuroticism score
setnames(P,old="20127-0.0",new="neuroticism.score")

###########

#filter_cols = c(diag.exclusion_cols, diag.surrogate_cols, diag.IBS_cols, sr.exclusion_cols, sr.IBS_cols)
filter_cols = grep("(diag)|(sr)|(neuroticism.score)",colnames(P),value=TRUE)
P[, lapply(.SD, sum), .SDcols=filter_cols]

save_zip(P[,.SD,.SDcols=c("eid",filter_cols)],"filters/diagnoses_selfreported.csv") 




