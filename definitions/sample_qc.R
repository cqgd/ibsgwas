library(data.table)
library(ggplot2)
setwd("~/work/ibs")

save_zip = function(object,path){
  fwrite(object,path)
  zip(paste0(path,".zip"),path)
}

#Load sample QC data

sqc = fread("qc/ukb_sqc_v2.txt")

#Find Brits passing initial QC

sqc[, in.PC.outlier.check:= 
      het.missing.outliers==0 &
      excluded.from.kinship.inference == 0 &
      excess.relatives == 0 & 
      used.in.pca.calculation == 1 & 
      in.white.British.ancestry.subset == 1]


#Find PC outliers (across dimensions) based on the centers and SDs among those 
#(Neale: only among those which also have phenotype data)

PCs = paste0("PC",1:6)
sds_threshold = 7

distance_sds_squared = function(x, sub){
  ( abs(mean(x[sub])-x) / sd(x[sub]) )^2
}

dists_uni_sdsq = sqc[, lapply(.SD, distance_sds_squared, sub = in.PC.outlier.check==TRUE), .SD=PCs]
dists_multi_sdsq = rowSums(dists_uni_sdsq)
sqc[, PCoutlier :=  dists_multi_sdsq > sds_threshold^2]

#Prepare link between sample QC and phenotype data, add eid's

link_path = c("fam/ukb17670_cal_chr22_v2_s488295.fam")
link = fread(link_path)

nrow(link)==nrow(sqc)
sqc[,eid:=link[,1]]

#Add phenotype data (includes ethnicity)

pheno = fread("pheno/ukb21760.csv")
included_ethnicities = c(1,1001,1002,1003) #See https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001
missing_ethnicities = c(-3,-1,NA)
sqcpheno = merge(sqc, pheno[,.(eid,`21000-0.0`)], by="eid", all.x=TRUE, all.y=FALSE)

#Obtain overview of filters

conditions = list(everyone = expression(rep(TRUE,.N)),
                  het.missing.outliers = expression(het.missing.outliers==0),
                  excluded.from.kinship.inference = expression(excluded.from.kinship.inference == 0),
                  excess.relatives = expression(excess.relatives == 0),
                  used.in.pca = expression(used.in.pca.calculation == 1), 
                  pc.outlier = expression(PCoutlier == FALSE),
                  ethnicity.included.or.missing = expression(`21000-0.0` %in% c(included_ethnicities, missing_ethnicities)),
                  consent = expression(eid>0))


condition_count = function(data,conditions){
  pass_alone = numeric(length(conditions))
  pass_cumulative = numeric(length(conditions))
  alone = matrix(nrow=nrow(data),ncol=length(conditions))
  colnames(alone) <- paste0("sqc.",names(conditions))
  
  cumulative = TRUE
  
  for(i in seq_along(conditions)){
    alone[,i] = data[,eval(conditions[[i]])]
    pass_alone[i] = sum(alone[,i])
    
    cumulative = cumulative & alone[,i]
    pass_cumulative[i] = sum(cumulative)
  }
  
  overview = data.table(filter=names(conditions), pass_alone, pass_cumulative)
  
  passing = data.table(alone, sqc.cumulative=cumulative)
  return(list(overview=overview,passing=passing))
}


qc = condition_count(sqcpheno, conditions)
print(qc$overview)
output = data.table(eid=sqcpheno$eid, qc$passing)

fwrite(qc$overview,"filters/sample_qc_overview.csv")
save_zip(output,"filters/sample_qc.csv")

#Overview of filters for supplement, noting that excess relatives is ultimately ignored:
conditions = list(everyone = expression(rep(TRUE,.N)),
                  het.missing.outliers = expression(het.missing.outliers==0),
                  excluded.from.kinship.inference = expression(excluded.from.kinship.inference == 0),
                  #used.in.pca = expression(used.in.pca.calculation == 1), 
                  pc.outlier = expression(PCoutlier == FALSE),
                  ethnicity.included.or.missing = expression(`21000-0.0` %in% c(included_ethnicities, missing_ethnicities)),
                  consent = expression(eid>0))
qc = condition_count(sqcpheno, conditions)
fwrite(qc$overview,"filters/sample_qc_overview_without_excess_relatives_or_used_in_pca.csv")

