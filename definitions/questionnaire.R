library(data.table)
setwd("~/work/ibs")
P = fread("pheno/ukb24282.csv")

NAtoF = function(v){
 v[is.na(v)]=FALSE
 return(v)
}

save_zip = function(object,path){
  fwrite(object,path)
  zip(paste0(path,".zip"),path)
}

#Gold-standard previous celiac disease / gluten sensitivity diagnosis
gold = -703:-705
q.cel.glut.gold = NAtoF(P$`21069-0.0`%in%gold)

#Cases
q1.case=P$`21025-0.0` %in% 3:6 #Updated to be more relaxed

q2.case=!(NAtoF(P$`21026-0.0` == 1)) #Cases can have NA here
 
q3.case=P$`21027-0.0` == 1

b_response = (-501):(-504)
Bullet = matrix(FALSE,nrow=nrow(P),ncol=3)
Bullet[,1] = P$`21028-0.0`%in%b_response
Bullet[,2] = NAtoF(P$`21029-0.0`%in%b_response) | NAtoF(P$`21030-0.0`%in%b_response) #Cases can have TRUE and NA
Bullet[,3] = NAtoF(P$`21031-0.0`%in%b_response) | NAtoF(P$`21032-0.0`%in%b_response)

q4to8.case = rowSums(Bullet,na.rm = TRUE)>=2

cases = q1.case & q2.case & q3.case & q4to8.case & !q.cel.glut.gold

#Controls
conts = P$`21025-0.0` %in% 0:1 & !q.cel.glut.gold
conts.strict = conts & (P$`21033-0.0` %in% (-500:-501)) & (P$`21034-0.0` %in% (-500:-501)) & !q.cel.glut.gold
  
#Subtypes
Types = matrix(FALSE,nrow=nrow(P),ncol=4)
colnames(Types) = paste0("type.",c("C","D","M","U"))

#No longer used, now done as per IBSMode.pdf
#rare_response = (-500):(-501) 
#freq_response = (-502):(-504)

never_response = -500
ever_response = (-501):(-504)

Types[,"type.C"] = P$`21033-0.0` %in% ever_response & P$`21034-0.0` %in% never_response
Types[,"type.D"] = P$`21033-0.0` %in% never_response & P$`21034-0.0` %in% ever_response
Types[,"type.M"] = P$`21033-0.0` %in% ever_response & P$`21034-0.0` %in% ever_response
#Types[,"type.U"] = !(Types[,"type.C"]) & !(Types[,"type.D"]) & !(Types[,"type.M"])
Types[,"type.U"] = P$`21033-0.0` %in% never_response & P$`21034-0.0` %in% never_response

#Combinations
case.C = cases & Types[,"type.C"]
case.D = cases & Types[,"type.D"]
case.M = cases & Types[,"type.M"]
case.U = cases & Types[,"type.U"]

#Functional
func = P$`21025-0.0`==0
func.C =  func & Types[,"type.C"] 
func.D = func & Types[,"type.D"]

func.relaxed = P$`21025-0.0`%in%c(0,1)
func.relaxed.C =  func.relaxed & Types[,"type.C"] 
func.relaxed.D = func.relaxed & Types[,"type.D"]

#Previously diagnosed with IBS
prev.diag.yes = NAtoF(P$`21024-0.0`==1)
prev.diag.no = NAtoF(P$`21024-0.0`==0)

#Participant
participant = !is.na(P$`21024-0.0`)


#Severity scoring system (SSS)
#Supplement 0 scores where follow-ups not relevant
P[`21035-0.0`==0,`21036-0.0`:=0]; P[`21035-0.0`==0,`21037-0.0`:=0]
P[`21038-0.0`==0,`21039-0.0`:=0]

SSS_fields = c("21036-0.0", "21037-0.0", "21039-0.0", "21040-0.0", "21041-0.0")
S = data.table( all_SSS = rowSums(P[,lapply(.SD,function(col){col%in%(0:10)}),.SDcols=SSS_fields])==length(SSS_fields) )
S[,oneplus_SSS := rowSums(P[,lapply(.SD,function(col){col%in%(0:10)}),.SDcols=SSS_fields])>=1]
S[,all_SSS := rowSums(P[,lapply(.SD,function(col){col%in%(0:10)}),.SDcols=SSS_fields])==length(SSS_fields)]
S[all_SSS==TRUE,SSS:=rowSums(P[S$all_SSS==TRUE,SSS_fields,with=FALSE])] #P and S in same order


# SSS visualization
# apply(P[,SSS_fields,with=FALSE],2,table)
# apply(P[all_SSS,SSS_fields,with=FALSE],2,table)
# cor(P[S$all_SSS==TRUE,SSS_fields,with=FALSE])
# ggplot(P[S$all_SSS==TRUE,.SD,.SDcols=SSS_fields], aes(x=factor(`21036-0.0`,levels=0:10),y=factor(`21040-0.0`,levels=0:10))) + 
#   geom_bin2d()

#quantile.normalize = function(v){
#  qnorm((frank(v,ties.method="average")-0.5)/length(v))
#}

#S[all_SSS==TRUE, SSS.qn:=quantile.normalize(SSS)]

#Bowel movement frequency
FQ_fields = c("hard"="21033-0.0",
              "loose"="21034-0.0",
              "daily.max"="21042-0.0",
              "weekly.min"="21043-0.0",
              "daily"="21044-0.0")
FQ = P[,FQ_fields,with=FALSE]
setnames(FQ,names(FQ_fields))
FQ[,hard:=abs(hard+500)][,loose:=abs(loose+500)]

#FQ[hard%in%(0:5),hard.qn := quantile.normalize(hard)] #These will stll be filtered by case.any, and then re-quantile normalized.
#FQ[loose%in%(0:5),loose.qn := quantile.normalize(loose)]
#FQ[daily%in%(0:20),daily.qn := quantile.normalize(daily)]
colnames(FQ) <- paste0("fq.",colnames(FQ))


#POST-INFECTIVE IBS
# We need
# Section G
# G1= yes = 01
# Plus
# Either G2=yes= 01
# Or 2 or more of
# G3 a-d = 01

pi_g1 = P$`21070-0.0` == -801 #Sudden onset
pi_g2 = P$`21071-0.0` == 1 #After infectious illness
pi_g3ad_fields = c(fever="21073-0.0", diarrhoea="21074-0.0",bloody_diarrhoea="21075-0.0",vomiting="21076-0.0") 
pi_g3ad = rowSums(P[,pi_g3ad_fields,with=FALSE]==1,na.rm = TRUE)>=2
post.infective = pi_g1 & (pi_g2 | pi_g3ad)

post.infective.C = post.infective & Types[,"type.C"]
post.infective.D = post.infective & Types[,"type.D"]
post.infective.M = post.infective & Types[,"type.M"]
post.infective.U = post.infective & Types[,"type.U"]

#PHQ12-SS is sum of  E1 a-n excluding E1k

#phq_info = fread("~/work/ibs/clinical/phq_fields.csv")
phq_info = fread("~/work/ibs/resource/IBS GWAS_ PHQ-9 , PHQ-12 and GAD-7 fields - Sheet1.csv")[questionnaire=="PHQ-12"]
phq_info[,name:=paste0(questionnaire,".",label)]
phq = as.matrix(P[,phq_info$field,with=FALSE])
colnames(phq) = phq_info$name
#Change N/A for intercourse pain in phq-12 to "not bothered at all", to change completeness
intercourse_na_response = phq[,"PHQ-12.dyspareunia"]==-313
phq[intercourse_na_response,"PHQ-12.dyspareunia"]=-600

phq_response = c(-600,-601,-602)
phq[!(phq%in%phq_response)]=NA
phq = -1*phq-600

oneplus_phq=rowSums(!is.na(phq),na.rm=TRUE)>=1 #A response of -313 to intercourse pain also counts as a start
#adds only 4 people total, who apparently had -313 as their only answer to this anywhere (and weren't previously included)
#all others who had -313 were already considered oneplus, they had at least one other answer

#all_phq=rowSums(!is.na(phq),na.rm=TRUE)==ncol(phq)

no_period_or_nausea = setdiff(colnames(phq),c("PHQ-12.period","PHQ-12.nausea"))
all_PHQ12=rowSums(!is.na(phq[,no_period_or_nausea]),na.rm=TRUE)==length(no_period_or_nausea) #ignore period/nausea questions when considering completeness
PHQ12_sum=rowSums(phq[,no_period_or_nausea],na.rm=TRUE) #also no longer take period/nausea questions into account for sum
phq = data.table(phq,oneplus_phq,all_PHQ12,PHQ12_sum)
phq[!all_PHQ12,PHQ12_sum:=NA]


risk_cols = c("21065-0.0","21067-0.0","21066-0.0","21062-0.0","21063-0.0")
risk_names = c("family.history","childhood.antibiotics","born.by.caesarean","treatment.anx","treatment.depr")
risk_factors = setNames(P[,lapply(.SD,function(x){x[!(x%in%(0:1))]=NA; return(x)}),.SDcols=risk_cols],risk_names) 

q_output = data.table(eid=P$eid,
                      participant,
                      cases, 
                      conts, conts.strict,
                      Types,
                      case.C, case.D, case.M, case.U, 
                      func.C, func.D,
                      func.relaxed.C, func.relaxed.D,
                      prev.diag.yes, prev.diag.no,
                      q.cel.glut.gold,
                      S,
                      FQ,
                      post.infective, post.infective.C, post.infective.D, post.infective.M, post.infective.U,
                      phq,
                      risk_factors)

colnames(q_output)[2:ncol(q_output)] <- paste0("Q_",colnames(q_output)[2:ncol(q_output)])
save_zip(q_output,"filters/questionnaire.csv")





