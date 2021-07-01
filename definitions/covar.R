setwd("~/work/ibs")
library(data.table)

# read in the phenotype data (for PCs etc)
phe2 <- fread('pheno/ukb21760.csv')
PCs = paste0("22009-0.",1:20)

#read in the sam file from the imputed data
samfile <- read.table('sample/ukb17670_imp_chr22_v3_s487327.sample')

#read in the filters/IBS definitions
disease <- fread('filters/all_integrated.csv')

#make maps between different files
eid <- as.character(phe2$"eid")
temp <- 1:length(eid)
names(temp) <- as.character(eid)
mapper <- temp[as.character(samfile[-(1:2),1])]

temp <- 1:length(disease$eid)
names(temp) <- as.character(disease$eid)
mapper2 <- temp[as.character(samfile[-(1:2),1])]

#make all the covariates
PCs <- as.matrix(phe2[mapper,paste0("22009-0.",1:20)])
colnames(PCs) <- paste0("PC",1:20)
sex <- as.numeric(samfile[-c(1:2),4])
age <- 2018-(phe2$"34-0.0")[mapper] #Should be 2017 for DHQ competion if further intrepreted
sexage <- sex*age
sexage <- sexage - mean(sexage,na.rm=T)
age <- age - mean(age,na.rm=T)
sex <- sex - mean(sex,na.rm=T)
agesq <- age^2
sexagesq <- sexage^2

#Quickly add subtypes, original 

disease[,final_case.sub.C:= Q_case.C & disease$final_case.any]
disease[,final_case.sub.U:= Q_case.U & disease$final_case.any]
disease[,final_case.sub.M:= Q_case.M & disease$final_case.any]
disease[,final_case.sub.D:= Q_case.D & disease$final_case.any]

#different cc scenarios

colnames(disease) <- gsub("^final_post.infective","final_case.post.infective",colnames(disease)) #won't be necessary after re-running integrate_filters.R

final_case_cols = c(grep("^final_case",colnames(disease),value=TRUE),"final_func.C","final_func.D")

cc_gen = function(case_col,cont_col,data){
  cc = rep(NA,nrow(data))
  cc[data[[case_col]]==TRUE]=1
  cc[data[[cont_col]]==TRUE]=0
  return(cc)
}

#find corresponding control columns, default to final_cont where not possible
final_cont_cols = gsub("final_case","final_cont",final_case_cols)
#final_cont_cols = gsub("final_anxcase","final_anxcont",final_cont_cols)
final_cont_cols[!final_cont_cols%in%colnames(disease)] = "final_cont"

ccs = mapply(cc_gen, final_case_cols, final_cont_cols, MoreArgs=list(data=disease))

#ccs = sapply(final_case_cols, cc_gen, cont_col="final_cont", data=disease) #old, would use final_cont for respondents GWAS, which is incorrect, introduced matching cont to case cols

#re-run all with relaxed controls
final_cont_cols.relaxed = rep("final_cont.relaxed",length(final_case_cols))
final_cont_cols.relaxed[final_case_cols=="final_case.no.anx"]="final_cont.no.anx.relaxed" #except for IBS w/o anxiety, which has its own relaxed control set without anxiety
ccs.relaxed = mapply(cc_gen, final_case_cols, final_cont_cols.relaxed, MoreArgs=list(data=disease))

#ccs.relaxed = sapply(final_case_cols, cc_gen, cont_col="final_cont.relaxed", data=disease)


ccs = ccs[mapper2,]
ccs.relaxed = ccs.relaxed[mapper2,]

quant_qc = disease$final_quant[mapper2]
quant_qc.mauro = disease$final_quant.mauro[mapper2]

colnames(ccs) = gsub("final_case","final_cc",colnames(ccs))
colnames(ccs.relaxed) = paste0(gsub("final_case","final_cc",colnames(ccs)),".conts.relaxed")

#make the output file
out <- data.table(FID=as.character(samfile[-c(1:2),1]),
                  IID=as.character(samfile[-c(1:2),2]),
                  missing=as.character(samfile[-c(1:2),3]),
                  sex.sample=as.character(samfile[-c(1:2),4]), #this is fine, Bolt-LMM only needs FID and IID
                  ccs, ccs.relaxed,
                  quant_qc, quant_qc.mauro,
                  age,sex,agesq,sexage,sexagesq,
                  PCs)

######
#Quickly add in sex-specific phenotype, this should be moved to the integrate_filters code
######
out[,final_cc.male:=final_cc.any][final_cc.any==1 & sex.sample!=1,final_cc.male:=NA] #eliminate non-males from cases only, keeping controls the same
out[,final_cc.female:=final_cc.any][final_cc.any==1 & sex.sample!=2,final_cc.female:=NA]

######

#Add SSS and bowel movement frequency data
Q = fread("filters/questionnaire.csv")
Q[,eid:=as.character(eid)]
out = merge(out, Q[,grepl("(^Q_fq)|(eid)|(Q_SSS$)|(Q_participant)|(phq)|(PHQ)|(type)|(post.infective)",colnames(Q)),with=FALSE],
      by.x="FID", by.y="eid",
      all.x=TRUE,
      sort=FALSE)

quantile.normalize = function(v){
  qnorm((frank(v,ties.method="average",na.last="keep")-0.5)/length(v[!is.na(v)]))
}


#Visualize these overlaps
out[(quant_qc==TRUE) & (Q_fq.loose %in% (0:4)),final_quant.loose:=quantile.normalize(Q_fq.loose)]
out[(quant_qc==TRUE) & (Q_fq.hard %in% (0:4)),final_quant.hard:=quantile.normalize(Q_fq.hard)]
#out[(quant_qc==TRUE) & (Q_fq.loose %in% (0:20)),final_quant.daily:=quantile.normalize(Q_fq.daily)] BUG, NOW FIXED BELOW, RE-RUN DAILY
out[(quant_qc==TRUE) & (Q_fq.daily %in% (0:20)),final_quant.daily:=quantile.normalize(Q_fq.daily)]
out[(quant_qc==TRUE) & (Q_fq.daily.max %in% (1:40)),final_quant.daily.max:=quantile.normalize(Q_fq.daily.max)]
out[(quant_qc==TRUE) & (Q_fq.weekly.min %in% (0:120)),final_quant.weekly.min:=quantile.normalize(Q_fq.weekly.min)]
out[(quant_qc==TRUE) & (is.finite(Q_SSS)),final_quant.SSS:=quantile.normalize(Q_SSS)]

#Mauro-like QC:
out[(quant_qc.mauro==TRUE) & (Q_fq.daily %in% (0:20)),final_quant.daily.mauro:=quantile.normalize(Q_fq.daily)]
out[(quant_qc.mauro==TRUE) & (Q_fq.daily.max %in% (1:40)),final_quant.daily.max.mauro:=quantile.normalize(Q_fq.daily.max)]
out[(quant_qc.mauro==TRUE) & (Q_fq.weekly.min %in% (0:120)),final_quant.weekly.min.mauro:=quantile.normalize(Q_fq.weekly.min)]

#Add in quantile normalization for PHQ-12 tiredness amongst pooled UKBB cases
out[final_cc.any==1 & Q_all_PHQ12==TRUE & `Q_PHQ-12.tired`%in%(0:2),final_quant.phq.tired:=quantile.normalize(`Q_PHQ-12.tired`)]

#Sanity check
out[,lapply(.SD,function(col){sum(!is.na(col))}),.SDcols=grepl("quant",colnames(out))]

#Qonly
#Restricted to Questionnaire participants (so cases and controls differ less wrt filling out Q)
out[Q_participant==TRUE,final_cc.Qonly.any:=final_cc.any]
out[Q_participant==TRUE,final_cc.Qonly.Qbased:=as.numeric(final_cc.Q | final_cc.Q.prev.diag.yes)]
out[Q_participant==TRUE,final_cc.Qonly.no.anx:=final_cc.no.anx]
#out[Q_participant==TRUE,final_cc.Qonly.no.anxbroad:=final_cc.no.anxbroad]

#Non-Q
out[Q_participant==FALSE, final_cc.Qnon.any:=final_cc.any.conts.relaxed]
out[Q_participant==FALSE, final_cc.Qnon.no.anx:=final_cc.no.anx.conts.relaxed]
#out[Q_participant==FALSE, final_cc.Qnon.no.anxbroad:=final_cc.no.anxbroad.conts.relaxed]

#Q_SS within IBS cases to address reviewer comment on whether mild and severe IBS are the same. Are there variants differentiating them?
out[is.finite(Q_SSS)==TRUE & final_cc.any==1, final_quant.caseSSS:=quantile.normalize(Q_SSS)]

#Only cases with family history
all(disease[mapper2,eid] == out$FID)
#out[disease[mapper2,Q_family.history]==TRUE,final_cc.any.fhx:=final_cc.any] #OLD always limited to respondents, as fhx is a DHQ question #Wrong, restricts controls also
#out[disease[mapper2,Q_family.history]==TRUE,final_cc.any.fhx:=final_cc.any] #OLD
#always limited to respondents, for both cases and controls, as fhx is a DHQ question, and cases will be solely respondents
out[final_cc.Qonly.any==1 & disease[mapper2,Q_family.history]==TRUE, final_cc.any.fhx := 1]
out[final_cc.Qonly.any==0, final_cc.any.fhx := 0]
table(out$final_cc.any.fhx)
table(out$final_cc.Qonly.any)

     
#Pure IBS
out[,final_cc.all:=final_cc.diag & final_cc.sr & final_cc.Q & final_cc.Q.prev.diag.yes]


#Add half-ukbb analyses, analyses without sample overlap (so LDSC isn't possibly influenced by this)
set.seed(123)
random_half = sample(1:nrow(out),round(nrow(out)/2))
out[,half_one:=NA][,half_two:=NA] #if you use false, get more controls upon combination
out[random_half,half_one:=TRUE]
out[is.na(half_one),half_two:=TRUE]

add_split_col = function(split_col,data){
  new_split_cols=paste0(split_col,c("_ONE","_TWO"))
  data[half_one==TRUE,(new_split_cols[1]):=get(split_col)]
  data[half_two==TRUE,(new_split_cols[2]):=get(split_col)]
}

split_cols = paste0("final_cc.",c("diag","sr","Q","Q.prev.diag.yes","Qonly.any","Qnon.any"))
sapply(split_cols,add_split_col,data=out)


#cat(paste(grep("(_TWO$)|(_ONE$)",colnames(out),value=TRUE)))
#final_cc.diag_ONE final_cc.diag_TWO final_cc.sr_ONE final_cc.sr_TWO final_cc.Q_ONE final_cc.Q_TWO final_cc.Q.prev.diag.yes_ONE final_cc.Q.prev.diag.yes_TWO final_cc.Qonly.any_ONE final_cc.Qonly.any_TWO final_cc.Qnon.any_ONE final_cc.Qnon.any_TWO

#QC: Check that the above GWASes are split by half  
out[,lapply(.SD,table),.SDcols=sort(grep(paste0(split_cols,collapse="|"),colnames(out),value=TRUE))]

#Questionnaire-derived subtypes applied to anyone who has a diagnosis in any of the sources
#Will be a subset of DHQ-respondents, and a subset of people who have been diagnosed
out[final_cc.any==1 & Q_type.C==1,final_cc.sub.any.C:=1]
out[final_cc.any==1 & Q_type.D==1,final_cc.sub.any.D:=1]
out[final_cc.any==1 & Q_type.M==1,final_cc.sub.any.M:=1]
out[final_cc.any==1 & Q_type.U==1,final_cc.sub.any.U:=1]

out[final_cc.any==0,final_cc.sub.any.C:=0]
out[final_cc.any==0,final_cc.sub.any.D:=0]
out[final_cc.any==0,final_cc.sub.any.M:=0]
out[final_cc.any==0,final_cc.sub.any.U:=0]

#QC: Check that subtype numbers have grown, controls stayed the same
out[,lapply(.SD,table),.SDcols=sort(grep("sub",colnames(out),value=TRUE))]

#D vs C directly
#Subtypes based on IBS but expanded to all IBS cases
out[final_cc.sub.any.D==1,final_cc.any.DvsC:=1]
out[final_cc.sub.any.C==1,final_cc.any.DvsC:=0]
#Subtypes restricted to questionnaire IBS
out[final_cc.sub.D==1,final_cc.DvsC:=1]
out[final_cc.sub.C==1,final_cc.DvsC:=0]

out[,lapply(.SD,table),.SDcols=sort(grep("DvsC",colnames(out),value=TRUE))]

#Within Rome-defintiion, i..e Rome (currently experiencing symptoms) + formal
out[final_cc.Q==1 & final_cc.diag==1,final_cc.Q.formal.diag:=1]
out[final_cc.Q==1 & final_cc.Q.prev.diag.yes==1,final_cc.Q.formal.Q.prev.diag.yes:=1]
out[final_cc.Q==1 & final_cc.sr==1,final_cc.Q.formal.sr:=1]

out[final_cc.Q==0, final_cc.Q.formal.diag:=0]
out[final_cc.Q==0, final_cc.Q.formal.Q.prev.diag.yes:=0]
out[final_cc.Q==0, final_cc.Q.formal.sr:=0]


out[,lapply(.SD,table),.SDcols=c("final_cc.Q",grep("formal",colnames(out),value=TRUE))]

#NATURE GENETICS REVISIONS
#SSS>300, same controls as always, must necessarily be on respondents only or you won't have SSS (so we use Qonly controls as well)
if(max(out$Q_SSS,na.rm=TRUE)>50.1){
  stop("You have 0-500 scale SSS, but are filtering on the scale of 0-50 below.")
}
out[final_cc.Qonly.any==1 & is.finite(Q_SSS) & Q_SSS>=30.0, final_cc.Qonly.any.severe:=1]
out[final_cc.Qonly.any==0, final_cc.Qonly.any.severe:=0]

#Twoplus (this will be a subset of cases we are certain of, same controls as always, meta-analyzed)
any_binary = as.matrix(out[,.SD,.SDcols=c("final_cc.diag", "final_cc.sr", "final_cc.Q", "final_cc.Q.prev.diag.yes")])
twoplus = rowSums(any_binary,na.rm=TRUE)>=2

out[final_cc.Qonly.any==1 & twoplus, final_cc.Qonly.any.twoplus:=1] 
out[final_cc.Qonly.any==0, final_cc.Qonly.any.twoplus:=0]

out[final_cc.Qnon.any==1 & twoplus, final_cc.Qnon.any.twoplus:=1] 
out[final_cc.Qnon.any==0, final_cc.Qnon.any.twoplus:=0]

#Add the dosages for out top HLA hit
#hla_hit, qctools expected dosage format
H = fread("hla/hla_SNPs_dosages.qctools.tsv.hit.transpose",header=FALSE)
Hs = fread("hla/hla_SNPs_dosages_column_eids.qctools.tsv",header=FALSE)
nrow(Hs) == nrow(H)-5
H = data.table(FID=as.character(Hs$V1), hla_hit_rs2736155_C=H[6:nrow(H),V1])

out = merge(out,H,by="FID",all.x=TRUE,sort=FALSE)

#Add the dosages for uc/microscopic colitis/celiac disease SNPs
H = fread("hla/related_traits_SNPs_dosages.qctools.tsv.transpose",header=FALSE)
Hs = fread("hla/related_traits_SNPs_dosages_column_eids.qctools.tsv",header=FALSE)
snp_names = paste0("related_trait_hit_",H[1,])
nrow(Hs) == nrow(H)-5
H = data.table(FID=as.character(Hs$V1), H[6:nrow(H)])
colnames(H)[2:ncol(H)] = snp_names

out = merge(out,H,by="FID",all.x=TRUE,sort=FALSE)

#Add HLA dosages for celiac disease (alternative to HLA allele proxy SNPs above)
#As well as our top HLA association (to check if SNP outside of HLA alleles is still significant)

HLA = fread("hla/ukb_hla_v2.txt")
fam = fread("fam/ukb17670_cal_chr22_v2_s488295.fam")
HLA[,FID:=as.character(fam$V1)]; rm(fam)
DQ2_5 = c("DQA1_501","DQB1_201")
DQ8 = c("DQA1_301","DQB1_302")
ibs_top_hla = "B_801"
HLA_celiac_ibs_top = HLA[,c("FID",DQ2_5,DQ8,ibs_top_hla),with=FALSE]

out = merge(out,HLA_celiac_ibs_top,by="FID",all.x=TRUE,sort=FALSE)

#hla_hit, plink minor count format (best guess, not expected dosage)
#H = fread("hla/ukb_imp_chr6_v3_hla_minorcount.raw")
#hla_hit_col = grep("rs2736155",colnames(H),value=TRUE)
#H = H[,.(as.character(FID), get(hla_hit_col))]
#colnames(H) = c("FID",paste0("hla_hit_",hla_hit_col))




#People not in Mauro's study, IBS-C
#Basically new cases from the DHQ only (therefore, no questionnaire bias meta-analysis is necessary)
#Also the pooled SR/ICD-10 cases, which are the ones mauro would have had
#sum(out$final_cc.any==1 & out$final_cc.sub.any.C==1,na.rm=TRUE)
#sum(out$final_cc.any==1 & out$final_cc.sub.any.C==1 & out$sex>0,na.rm=TRUE)
#sum(out$final_cc.any==1 & out$final_cc.sub.any.C==1 & out$sex>0 & out$final_cc.sr%in%c(NA,0) & out$final_cc.diag%in%c(NA,0),na.rm=TRUE) #check sex coding, 1 is female
#sum(out$final_cc.any==1 & out$final_cc.sub.any.C==1 & out$sex<0 & (out$final_cc.sr==1 & out$final_cc.diag==1),na.rm=TRUE) #check sex coding
out[final_cc.any==0,final_cc.nomauro.sub.any.C.female:=0]
out[final_cc.any==0,final_cc.likemauro.sub.any.C.female:=0]
out[final_cc.any==0,final_cc.nomauro.any.female:=0]
out[final_cc.any==0,final_cc.likemauro.any.female:=0]

out[final_cc.any==1 & final_cc.sub.any.C==1 & sex>0 & final_cc.sr%in%c(NA,0) & final_cc.diag%in%c(NA,0), final_cc.nomauro.sub.any.C.female:=1]
out[final_cc.any==1 & final_cc.sub.any.C==1 & sex>0 & (final_cc.sr==1 | final_cc.diag==1), final_cc.likemauro.sub.any.C.female:=1]
out[final_cc.any==1 & sex>0 & final_cc.sr%in%c(NA,0) & final_cc.diag%in%c(NA,0), final_cc.nomauro.any.female:=1]
out[final_cc.any==1 & sex>0 & (final_cc.sr==1 | final_cc.diag==1), final_cc.likemauro.any.female:=1]


#it appears that bonfiglio actually just used sr, not icd-10 and sr

out[final_cc.any==0,final_cc.likebonfiglio.sub.any.C.female:=0]
out[final_cc.any==0,final_cc.nobonfiglio.sub.any.C.female:=0]
out[final_cc.any==0,final_cc.likebonfiglio.any.female:=0]
out[final_cc.any==0,final_cc.nobonfiglio.any.female:=0]

out[final_cc.any==1 & final_cc.sub.any.C==1 & sex>0 & final_cc.sr%in%c(NA,0), final_cc.nobonfiglio.sub.any.C.female:=1]
out[final_cc.any==1 & final_cc.sub.any.C==1 & sex>0 & final_cc.sr==1, final_cc.likebonfiglio.sub.any.C.female:=1]
out[final_cc.any==1 & sex>0 & final_cc.sr%in%c(NA,0), final_cc.nobonfiglio.any.female:=1]
out[final_cc.any==1 & sex>0 & final_cc.sr==1, final_cc.likebonfiglio.any.female:=1]


table(out$final_cc.nomauro.sub.any.C.female)
table(out$final_cc.likemauro.sub.any.C.female)
table(out$final_cc.nomauro.any.female)
table(out$final_cc.likemauro.any.female)
table(out$final_cc.nobonfiglio.sub.any.C.female)
table(out$final_cc.likebonfiglio.sub.any.C.female)
table(out$final_cc.nobonfiglio.any.female)
table(out$final_cc.likebonfiglio.any.female)

#More Mauro requests:

add_female_mauroreq = function(col,data){
  if(! col %in% colnames(data) ){ 
    stop(cat(col,"not among columns, can't add mauroreq female version")) 
  }
  new_col = paste0(col,".female.mauroreq")
  if( new_col %in% colnames(data)){
    cat(new_col," exists already, removing...\n")
    data[,(new_col):=NULL]
  }
  data[get(col)==1 & sex>0,(new_col):=1]
  data[get(col)==0 & sex>0,(new_col):=0]
  return(new_col)
}

out[final_cc.sr==1 & final_cc.Qonly.any==1, final_cc.Qonly.sr := 1]
out[final_cc.Qonly.any==0, final_cc.Qonly.sr := 0]

out[final_cc.sr==1 & final_cc.Qnon.any==1, final_cc.Qnon.sr := 1]
out[final_cc.Qnon.any==0, final_cc.Qnon.sr := 0]
    
mauroreq = c("final_cc.Q",
"final_cc.sub.C",
"final_cc.sub.D",
"final_cc.sub.M",
"final_cc.sr",
"final_cc.Q.prev.diag.yes",
"final_cc.diag",
"final_cc.Qnon.sr",
"final_cc.Qonly.sr",
#New, the below had men in the controls -- we can rerun
"final_cc.likebonfiglio.any.female",
"final_cc.nobonfiglio.any.female",
"final_cc.likebonfiglio.sub.any.C.female",
"final_cc.nobonfiglio.sub.any.C.female",
"final_cc.female")

mauroreq_cols = character(length(mauroreq))
for(i in seq_along(mauroreq)){ #hopefully less RAM than lapply
  mauroreq_cols[i] = add_female_mauroreq(col=mauroreq[i], data=out)
}
#Anxiety without IBS
#Anxiety with IBS

#Adding GAD-7


#Anxiety single deifntion GWASes
phq9_gad7 = fread("pheno/ukb39562_GAD7_PHQ9.csv",key="eid")
new = colnames(phq9_gad7)[!colnames(phq9_gad7) %in% colnames(out)]
phq9_gad7[,FID:=as.character(eid)]
phq9_gad7 = phq9_gad7[,.SD,.SDcols=c("FID",new)]
out = merge(out,phq9_gad7,all.x=TRUE,by="FID",sort=FALSE)

diag_sr = fread("filters/diagnoses_selfreported.csv", key="eid")
diag_sr[,FID:=as.character(eid)]
out = merge(out,diag_sr,all.x=TRUE,by="FID",sort=FALSE)

quest = fread("filters/questionnaire.csv", key="eid")
quest = quest[,.(FID=as.character(eid),Q_treatment.anx)]
out = merge(out, quest, all.x=TRUE, by = "FID", sort=FALSE)

restrict_cases = function(col,data){
  base_col = "final_cc.anxbroad"
  if(! col %in% colnames(data) ){ 
    stop(cat(col,"not among columns, can't restrict cases based on this.")) 
  }
  new_col = paste0(base_col,"_casesfrom_",col)
  if( new_col %in% colnames(data)){
    cat(new_col," exists already, removing...\n")
    data[,(new_col):=NULL]
  }
  data[get(base_col)==0, (new_col):=0] #shared controls
  data[get(base_col)==1 & get(col)==1,(new_col):=1]
  return(new_col)
}

out[,diag.interest_F40_or_F41:=diag.interest_F40 | diag.interest_F41]
anxdefs = c("diag.interest_F40_or_F41",
"diag.interest_F40", #By itself, F40 has very few people
"diag.interest_F41", 
"sr.interest_1287",
"Q_treatment.anx",
"GAD7.anx")

anxdef_cols = character(length(anxdefs))
for(i in seq_along(anxdef_cols)){ #hopefully less RAM than lapply
  anxdef_cols[i] = restrict_cases(col=anxdefs[i], data=out)
}

#Write some function that takes anxbroad GWAS (so controls restricted to GAD respondents, cases not necessarily).
#Then limits cases to just this

#write the output
write.table(out,file="covar/cc_cov.sample",sep="\t",col.names=T,row.names=F,quote=F)

#delete diag.sr, sr, and relaxed equivalents
#final_cc.Qnon.any final_cc.respondent final_cc.all 
#also re-run diag, sr and relaxed equivalents

#cov_old = fread("covar/cc_cov.sample")
#all(which(out$final_cc.any==1)==which(cov_old$CC==1))
#[1] TRUE

#shared = intersect(colnames(Q),colnames(phe2))
#M = merge(Q[,shared,with=FALSE],phe2[,shared,with=FALSE],by="eid",suffixes=c(".q",".phe2"))

#Sense-checking new nat gen revision gwas numbers
#Severe 
#mean(out[final_cc.Qonly.any==1 & is.finite(Q_SSS),Q_SSS>30.0]) #10% of cases with SSS>30
#table(out$final_cc.Qonly.any) 
#table(out$final_cc.Qonly.any.severe) #cases shrunk to about 10%, good
#twoplus
#table(out$final_cc.Qonly.any.twoplus)
#table(out$final_cc.Qonly.any)
#table(out$final_cc.Qnon.any.twoplus) #A lot more thinning in the non-respondent sample
#table(out$final_cc.Qnon.any)
#sum(twoplus==1) #Cases well split

#Sense-checking mauroreq cols
t(out[,sapply(.SD,table),.SDcols=mauroreq_cols])
paste0(mauroreq_cols,collapse=" ")

#Sense-checking anxdef cols
t(out[,sapply(.SD,table),.SDcols=anxdef_cols])
paste0(anxdef_cols,collapse=" ")
