library(data.table)
setwd("E:/Work/ibs")

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

P = fread("pheno/ukb24282.csv")

S = fread("resource/IBS_filter_codes - surgeries.csv")

surg_cols = grep("41200|41210",colnames(P),value = TRUE)

s_codes = S[category=="exclusion",OPCS4.code]
print(length(s_codes))
for(s in 1:length(s_codes)){
  print(s_codes[s])
  print(Sys.time())
  check_code(s_codes[s],P,cols=surg_cols,col_prefix="surg.exclusion",starts_with=TRUE)  
}

check_code(S[category=="gall",OPCS4.code],P,cols=surg_cols,col_prefix="surg.gall",starts_with=TRUE)

filter_cols = grep("(surg.)",colnames(P),value=TRUE)
save_zip(P[,.SD,.SDcols=c("eid",filter_cols)],"filters/surgeries.csv")

#fwrite(P,"filters/surg_ALL.csv")
#Should be faster than the loop, but somehow freezes:
#surg.exclusion_cols = sapply(S[category=="exclusion",OPCS4.code], check_code, P=P, cols=surg_cols, col_prefix="surg.exclusion", starts_with=TRUE)
#P[, lapply(.SD, sum), .SDcols=unlist(surg.exclusion_cols)]
