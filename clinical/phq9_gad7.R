library(data.table)

#Load phenotype data and questionnaire field mappings
D = fread("/home/cq/work/ibs/pheno/ukb39562.csv")
#M = fread("/home/cq/work/ibs/resource/GAD7_PHD9_UKBBmaps.csv")
M = fread("/home/cq/work/ibs/resource/IBS GWAS_ PHQ-9 , PHQ-12 and GAD-7 fields - Sheet1.csv")

PHQ9_fields = M[questionnaire=="PHQ-9",field]
GAD7_fields = M[questionnaire=="GAD-7",field]

#Subset to PHQ-9/GAD-7 info
all(c(PHQ9_fields,GAD7_fields)%in%colnames(D))

#Add anxiety diagnosis for anxiety without IBS and IBS without anxiety re-analysis (Lancet reviewer 1)


S = D[,.SD,.SDcols=c("eid",c(PHQ9_fields,GAD7_fields))]; rm(D)

#Subtract 1 from coding 504, so one can sum for scores (originally response 1 == not at all == 0 pts)


S[,(PHQ9_fields):=(.SD-1),.SDcols=PHQ9_fields]
S[,(GAD7_fields):=(.SD-1),.SDcols=GAD7_fields]

#Assess completion and sum scores
S[,oneplus_PHQ9 := rowSums(S[,lapply(.SD,function(col){col%in%(0:3)}),.SDcols=PHQ9_fields])>=1]
S[,all_PHQ9 := rowSums(S[,lapply(.SD,function(col){col%in%(0:3)}),.SDcols=PHQ9_fields])==length(PHQ9_fields)]
S[all_PHQ9==TRUE,PHQ9:=rowSums(S[S$all_PHQ9==TRUE,PHQ9_fields,with=FALSE])]

S[,oneplus_GAD7 := rowSums(S[,lapply(.SD,function(col){col%in%(0:3)}),.SDcols=GAD7_fields])>=1]
S[,all_GAD7 := rowSums(S[,lapply(.SD,function(col){col%in%(0:3)}),.SDcols=GAD7_fields])==length(GAD7_fields)]
S[all_GAD7==TRUE,GAD7:=rowSums(S[S$all_GAD7==TRUE,GAD7_fields,with=FALSE])]

setnames(S,old=M$field,new=paste0(M$questionnaire,".",M$label),skip_absent = TRUE)

S[all_GAD7==TRUE,GAD7.anx:=GAD7>=10]
#save subset
fwrite(S,"~/work/ibs/pheno/ukb39562_GAD7_PHQ9.csv")



