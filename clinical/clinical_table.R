setwd("~/work/ibs/")
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

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

stats_for_col = function(col, data){
  if(is.numeric(col)){ #row numbers directly
    subset = logical(nrow(data))
    subset[col] = TRUE
  }else if(is.character(col)){ #rows where column==1
    subset = data[[col]]==1
  }else{
    quit("error")
  }
  
  part_one = data[subset, .(
      N=.N,
      age=mean(age_dhq,na.rm=TRUE),
      age_sd=sd(age_dhq,na.rm=TRUE),
      male=100*sum(sex.sample==1,na.rm=TRUE)/.N,
      female=100*sum(sex.sample==2,na.rm=TRUE)/.N,
      Q_participation=100*mean(Q_participant,na.rm=TRUE))]
  
  sss_cols = grep("^Q_SSS.",colnames(data),value=TRUE)
  sss.1 = data[subset,.(`SSS_started` = 100*mean(Q_oneplus_SSS,na.rm=TRUE), `SSS_complete` = 100*mean(Q_all_SSS,na.rm=TRUE))]
  sss.2 = data[subset & Q_all_SSS, lapply(.SD,mean,na.rm=TRUE), .SDcols=sss_cols] #normally na.rm not needed, only for period question
  colnames(sss.2) <- gsub("^Q_","",colnames(sss.2))
  sss.3 = data[subset & Q_all_SSS, .(mean_SSS=mean(Q_SSS,na.rm=TRUE))]
  part_sss = cbind(sss.1,sss.2,sss.3)
  
  extra_surg_cols = grep("surg.extra",colnames(data),value=TRUE)
  part_three = data[subset, lapply(.SD,function(col){100*sum(col==1,na.rm=TRUE)/.N}),.SDcols=extra_surg_cols]
  
  risk_cols = c("21065-0.0","21067-0.0","21066-0.0","21062-0.0","21063-0.0")
  risk_names = c("family.history","childhood.antibiotics","born.by.caesarean","treatment.anx","treatment.depr")
  part_four = data[subset & Q_participant, lapply(.SD,function(col){100*sum(col==1,na.rm=TRUE)/.N}),.SDcols=risk_cols] #only amongst participants
  names(part_four) = risk_names
  
  
  #PHQ-12
  phq_cols = grep("^Q_PHQ-12.",colnames(data),value=TRUE)
  phq.1 = data[subset,.(`PHQ-12_started` = 100*mean(Q_oneplus_phq,na.rm=TRUE), `PHQ-12_complete` = 100*mean(Q_all_PHQ12,na.rm=TRUE))]
  phq.2 = data[subset & Q_all_PHQ12, lapply(.SD,mean,na.rm=TRUE), .SDcols=phq_cols] #normally na.rm not needed, only for period question
  colnames(phq.2) <- gsub("^Q_","",colnames(phq.2))
  phq.3 = data[subset & Q_all_PHQ12, .(mean_PHQ12=mean(Q_PHQ12_sum,na.rm=TRUE))]
  part_phq = cbind(phq.1,phq.2,phq.3)
  
  #GAD-7
  gad_cols = grep("GAD-7\\.",colnames(data),value=TRUE)
  gad.1 = data[subset,.(gad_started = 100*mean(oneplus_GAD7,na.rm=TRUE), gad_complete = 100*mean(all_GAD7,na.rm=TRUE))]
  gad.2 = data[subset & all_GAD7, lapply(.SD,mean,na.rm=TRUE), .SDcols=gad_cols]
  gad.3 = data[subset & all_GAD7, .(mean_gad=mean(GAD7,na.rm=TRUE))]
  part_gad = cbind(gad.1,gad.2,gad.3)
  
  #PHQ-9
  phq9_cols = grep("PHQ-9\\.",colnames(data),value=TRUE)
  phq9.1 = data[subset,.(phq9_started = 100*mean(oneplus_PHQ9,na.rm=TRUE), phq9_complete = 100*mean(all_PHQ9,na.rm=TRUE))]
  phq9.2 = data[subset & all_PHQ9, lapply(.SD,mean,na.rm=TRUE), .SDcols=phq9_cols]
  phq9.3 = data[subset & all_PHQ9, .(mean_phq9=mean(PHQ9,na.rm=TRUE))]
  part_phq9 = cbind(phq9.1,phq9.2,phq9.3)
  
  interest_cols = grep("interest",colnames(data),value=TRUE)
  part_immune = data[subset, lapply(.SD,function(col){100*sum(col==1,na.rm=TRUE)/.N}),.SDcols=interest_cols]
  
  townsend = data[subset,mean(`189-0.0`,na.rm=TRUE)]
  
  #if(nrow(part_two)==0){part_two=data.table(mean_SSS=NA)}
  result = cbind(part_one, part_sss, part_four, part_phq, part_gad, part_phq9, part_three,part_immune,townsend)
  return(result)
}

townsend_only = function(col,data){
    if(is.numeric(col)){ #row numbers directly
      subset = logical(nrow(data))
      subset[col] = TRUE
    }else if(is.character(col)){ #rows where column==1
      subset = data[[col]]==1
    }else{
      quit("error")
    }
  
  townsend = data[subset,mean(`189-0.0`,na.rm=TRUE)]
  return(townsend)
    
}
  
dists_for_col = function(col,data){
  subset = data[[col]]==1
  p1 = data[subset & Q_all_SSS, .(variable="Q_SSS", value=Q_SSS)]
  p2 = data[subset & Q_participant, .(absent.work.weeks=`21045-0.0`, affected.work.weeks=`21047-0.0`)]
  p2 = melt(p2, measure.vars=colnames(p2))
  p2 = p2[value>=0 & value<=52]
  p3 = data[subset & Q_all_phq, .(variable="Q_phq_sum", value=Q_phq_sum)]
  result = data.table(group=col, rbind(p1,p2,p3))
  return(result)
}

#Read data
I = fread("/home/cq/work/ibs/filters/all_integrated.csv",key="eid")
C = fread("/home/cq/work/ibs/covar/cc_cov.sample"); colnames(C)[1] <- "eid"; setkey(C,"eid")
P = fread("pheno/ukb24282.csv",key="eid")

A = fread("pheno/ukb21760.csv",key="eid")[,.(eid,`34-0.0`)]
A[,age_dhq:=2017-`34-0.0`]

E = fread("/home/cq/work/ibs/pheno/ukb39562_GAD7_PHQ9.csv",key="eid")
consecutive_merge = function(a,b){
  new_in_b = c("eid",colnames(b)[!(colnames(b)%in%colnames(a))])
  b = b[,.SD,.SDcols=new_in_b]
  merge(a,b,key="eid",all=TRUE,no.dups=FALSE)
}
X = Reduce(consecutive_merge,list(I,C,P,A,E))
#X[final_cc.respondent==TRUE,tstrsplit(`21023-0.0`,"-")][,table(as.numeric(V1))]# 2017: 153636 & 2018: 2935, only 1 percent completed outside of 2017
#X = I[C,][P,][A,][E,] #Does not retain all=TRUE
rm(I);rm(C);rm(P);rm(A);rm(E);
gc()


Ptwo = fread("pheno/ukb7725.csv",key="eid")
adopted_alive_cols = grep("(^3942-)|(^3912-)",colnames(Ptwo),value=TRUE)
townsend_col = "189-0.0"

Ptwo_mini = Ptwo[,.SD,.SDcols=c("eid",townsend_col,adopted_alive_cols)]
adopted_proxy = rowSums(Ptwo_mini[,lapply(.SD,function(col){!is.na(col)}),.SDcols=adopted_alive_cols])>=1
sum(adopted_proxy)
Ptwo_mini[,adopted_proxy:=adopted_proxy]; rm(adopted_proxy)
#Cases/conts out of adopted.
#Add adopted and townsend to stats, as well as percentage of responses greater than 1
#How does this compare to number of adopted people?
X = merge(X,Ptwo_mini,by="eid",all=TRUE,no.dups=FALSE)
#X = X[Ptwo_mini,]
rm(Ptwo); rm(Ptwo_mini)
gc()

#S = fread("sample/ukb17670_imp_chr22_v3_s487327.sample")


#Add surgical codes of interest
surg_cols = grep("41200|41210",colnames(X),value = TRUE)
surgs = fread("/home/cq/work/ibs/clinical/surgeries.txt")
for(s in 1:nrow(surgs)){
  check_code(surgs$code[s],X,cols=surg_cols,col_prefix=paste0("surg.extra.",surgs$name[s]),starts_with=TRUE)
}

#Pool codes: asthma, appendicits
add_pooled = function(columns, data){
  data[,(columns[3]) := (get(columns[1])==TRUE | get(columns[2])==TRUE)]
}

add_pooled_atopy = function(columns, data){ #just make one function with e.g. expression sep="|"
  data[,(columns[4]) := (get(columns[1])==TRUE | get(columns[2])==TRUE | get(columns[3])==TRUE)]
}

asthma = c("diag.interest_J45", "sr.interest_1111","pooled.interest_asthma")
eczema = c("diag.interest_L20", "sr.interest_1452","pooled.interest_eczema")
hayfever = c("diag.interest_J301", "sr.interest_1387","pooled.interest_hayfever")
appendicitis = c("surg.extra.appendix.excision_H01", "surg.extra.appendix.excision_H02", "surg.extra.appendix.excision.pooled_H01_H02")
hysterectomy = c("surg.extra.hysterectomy_Q07","surg.extra.hysterectomy_Q08","surg.extra.hysterectomy.pooled_Q07_Q08")
atopy = c("pooled.interest_asthma","pooled.interest_eczema","pooled.interest_hayfever","pooled.interest_atopy")

add_pooled(asthma,X)
add_pooled(eczema,X)
add_pooled(hayfever,X)
add_pooled(appendicitis,X)
add_pooled(hysterectomy,X)
add_pooled_atopy(atopy,X)


X[final_cc.any==1 & Q_all_SSS==TRUE & Q_SSS<17.5, final_case.SSS.mild:=1]
X[final_cc.any==1 & Q_all_SSS==TRUE & Q_SSS>=17.5 & Q_SSS<=30.0, final_case.SSS.moderate:=1]
X[final_cc.any==1 & Q_all_SSS==TRUE & Q_SSS>30.0, final_case.SSS.severe:=1]

X[final_cc.Qnon.any==0,final_cont.Qnon:=1]

# cont_cols = c("final_cont","final_cont.relaxed")
# case_cols = grep("^final_cc",colnames(X),value=TRUE)
# case_cols = case_cols[!grepl("(relaxed)|(male)|(Qonly)|(Qnon)|(respondent)|(ONE)|(TWO)|(DvsC)",case_cols)]
# func_cols = grep("^final_func",colnames(X),value=TRUE)
# pi_cols = grep("^final_post.infective",colnames(X),value=TRUE)
# subset_cols = c(cont_cols,func_cols,pi_cols,case_cols)
groups = fread("resource/IBS clinical call table - group_order.csv",header = TRUE)
subset_cols = groups$group

print(subset_cols)
#Then re-do table with "subset = data[[col]]==1 & Q.participant=TRUE"

#in case fo error, debug:
#for(i in seq_along(subset_cols)){print(i);print(subset_cols[i]);stats_for_col(subset_cols[i],data=X)}
#

out = sapply(subset_cols,stats_for_col,data=X)
out = t(out)
out = data.table(out,keep.rownames = TRUE)
colnames(out)[1] <-"group"
out = merge(out,groups,all.x=TRUE,sort=FALSE)
setcolorder(out,c("group","description","main"))

fwrite(X,"~/work/ibs/clinical/clinical_dataset.csv") #save time rebuilding
X = fread("~/work/ibs/clinical/clinical_dataset.csv")

fwrite(out,"~/work/ibs/clinical/clinical_table.csv")
clinical_main = out[main==1,.(group,description,N,male,female,age,age_sd,mean_SSS,mean_PHQ12,mean_gad,family.history,childhood.antibiotics,born.by.caesarean,treatment.anx,treatment.depr,pooled.interest_atopy)]
fwrite(clinical_main,"~/work/ibs/clinical/clinical_table_main.csv")


#Adding asterisks for significant differences to a copy of the clinical table
fields = c("mean_SSS"="Q_SSS",
           "mean_PHQ12"="Q_PHQ12_sum",
           "mean_gad"="GAD7",
           "family.history"="Q_family.history",
           "childhood.antibiotics"="Q_childhood.antibiotics",
           "born.by.caesarean"="Q_born.by.caesarean",
           "treatment.anx"="Q_treatment.anx",
           "treatment.depr"="Q_treatment.depr",
           "pooled.interest_atopy"="pooled.interest_atopy") #Name these based on ln216
ccs = groups[main==1 & group!="final_cont",group]

#if this has no controls, create a manual cc vector with the default controls
#needed for e.g. final_case.SSS.mild

args = data.table(expand.grid(ccs,fields,stringsAsFactors = FALSE))

glm_test = function(cc,field,X){
  formula = as.formula(paste0(cc,"~",field,"+sex+age+agesq+sexage+sexagesq+Q_participant"))
  m = glm(data=X,formula,family="binomial")
  summ = summary(m)$coefficients
  pval = summ[2,"Pr(>|z|)"]
  testp = c(cc=cc,field=field,pval=pval)
  return(testp)
}

R = mapply(glm_test,args[[1]],args[[2]],MoreArgs=list(X=X))
R = data.table(t(R))
R[,pval:=as.numeric(pval)]
sig_threshold=0.05/nrow(R)
R[,sig:=pval<sig_threshold]
field_names = data.table(field=fields,field.name=names(fields))
R = merge(R,field_names,by="field",all=TRUE)
CM = copy(clinical_main)

two_decimal_cols = unique(R$field.name[!grepl("^median_",R$field.name)])
two_decimal_cols = unique(R$field.name)

for(col in two_decimal_cols){
  CM[,(col):=as.character(formatC(unlist(get(col)),1,format="f"))]
}

for(col in unique(R$field.name)){
  CM[,(col):=as.character(unlist(get(col)))]
}

for(i in which(R$sig)){
  current_value = CM[group==R[i,cc], get(R[i,field.name])]
  CM[group==R[i,cc], (R[i,field.name])] = paste0(current_value,"*")
}

fwrite(CM,"~/work/ibs/clinical/clinical_table_main_asterisks.csv")

#END CONSTRUCTION

#PHQ-12 and SSS correlation amongst cases and controls

X[final_case.any==1 & Q_all_SSS==TRUE & Q_all_PHQ12==TRUE,cor.test(Q_SSS,Q_PHQ12_sum)]
X[final_case.any==1 & Q_all_SSS==TRUE & Q_all_PHQ12==TRUE,sd(Q_SSS)]
X[final_case.any==1 & Q_all_SSS==TRUE & Q_all_PHQ12==TRUE,sd(Q_PHQ12_sum)]

X[final_case.any==1 & Q_all_SSS==TRUE & Q_all_PHQ12==TRUE,sd(Q_PHQ12_sum)]


X[Q_all_PHQ12==TRUE & final_case.any==1,quantile(Q_PHQ12_sum)]
# 0%  25%  50%  75% 100% 
# 0    4    6    9   22 
X[Q_all_PHQ12==TRUE & final_case.any==0,quantile(Q_PHQ12_sum)]
# 0%  25%  50%  75% 100% 
# 0    2    4    6   22 

X[Q_all_PHQ12==TRUE,t.test(Q_PHQ12_sum~final_case.any)]

X[final_case.any==1 & Q_all_SSS==TRUE & Q_all_phq==TRUE,.(.N,cor(Q_SSS,Q_phq_sum))]
X[final_case.any==1 & Q_all_SSS==TRUE & Q_all_phq==TRUE,cor.test(Q_SSS,Q_phq_sum)]
#0.3982718 [95% CI: 0.3889246 - 0.4075370] amongst 31402 (pooled) IBS cases who completed all PHQ-12 and SSS questions.
X[final_cont==1 & Q_all_SSS==TRUE & Q_all_phq==TRUE,.(.N, cor(Q_SSS,Q_phq_sum))]
X[final_cont==1 & Q_all_SSS==TRUE & Q_all_phq==TRUE,cor.test(Q_SSS,Q_phq_sum)]
#The PHQ-12 SSS correlation is 0.3982718 [95% CI: 0.3889246 - 0.4075370] amongst 31402 (pooled) IBS cases who completed all PHQ-12 and SSS questions. Amongst 70172 (strict) controls, it's 0.2110893 [95% CI: 0.2040090 - 0.2181475].

#PHQ-12 and IBS-SSS correlation with controls
Z = X[(final_case.any==1 | final_cont==1) & Q_all_SSS==TRUE,.(cc=(final_case.any==1),SSS=Q_SSS,PHQ12=Q_PHQ12_sum)]
Z[,id:=1:.N]
subset = c(Z[cc==TRUE,sample(id,20e3)],
           Z[cc==FALSE,sample(id,20e3)])
Z[cc==TRUE,cc_name:="Case"][cc==FALSE,cc_name:="Control"]
p.SSS_PHQ12 = ggplot(Z[sample(subset)],aes(x=SSS*10,y=PHQ12,col=cc_name)) +
  geom_jitter(alpha=0.05) +
  geom_smooth(data=Z) +
  ylab("PHQ-12") +
  xlab("IBS-SSS") +
  labs(col="IBS")
p.SSS_PHQ12
ggsave("~/work/ibs/plots/SSS_PHQ12_corr.pdf",p.SSS_PHQ12,width=9,height=7)
ggplot(X[(final_case.any==1 | final_cont==1) & Q_all_SSS==TRUE], aes(x=Q_SSS,y=Q_PHQ12_sum,col=final_case.any==1)) + 
         geom_smooth()

#Simplified, without controls
subset = Z[cc==TRUE,sample(id,20e3)]
p.SSS_PHQ12 = ggplot(Z[sample(subset)],aes(x=SSS*10,y=PHQ12,col=cc_name)) +
  geom_jitter(alpha=0.05) +
  geom_smooth(data=Z[cc==TRUE]) +
  ylab("PHQ-12") +
  xlab("IBS-SSS") +
  labs(col="IBS") +
  theme(legend.position="none")
p.SSS_PHQ12
ggsave("~/work/ibs/plots/SSS_PHQ12_corr_simplified.pdf",p.SSS_PHQ12,width=7,height=7)
ggplot(X[(final_case.any==1 | final_cont==1) & Q_all_SSS==TRUE], aes(x=Q_SSS,y=Q_PHQ12_sum,col=final_case.any==1)) + 
  geom_smooth()


##Permutation testing


##Permutation testing

#Group 1 defined by col
#Group 2 defined by col

#Matrix where each col is a group setup, col0 is group definition, AAAABBB, col2 is 

# random_groups = function(col.A,col.B,data,n_nulls=100){
#   subset.A = which(data[[col.A]]==1)
#   subset.B = which(data[[col.B]]==1)
#   #group = as.factor(rep(c(col.A,col.B),c(length(subset.A),length(subset.B)))) #doesn't work with same group twice as test case
#   group = as.factor(rep(LETTERS[1:2],c(length(subset.A),length(subset.B))))
#   obs = c(subset.A,subset.B)
#   null = replicate(n_nulls,sample(obs))
#   setups = data.table(null,obs,group)
#   colnames(setups)[1:n_nulls] = paste0("null.",1:n_nulls)
#   return(setups)
# }
# 
# 
# test_stat = function(indices,group,data){
#   subset.A = indices[as.numeric(group)==1]
#   subset.B = indices[as.numeric(group)==2]
#   # stats.A = stats_for_col(subset.A,data)
#   # stats.B = stats_for_col(subset.B,data)
#   stats.A = townsend_only(subset.A,data)
#   stats.B = townsend_only(subset.B,data)
#   test_stat = stats.A-stats.B
#   return(test_stat)
# }
# 
# n_nulls = 8e3
# setups = random_groups(subset_cols[1],subset_cols[16],data=X,n_nulls)
# null_cols = grep("^null.",colnames(setups),value=TRUE)
# nulls = rbindlist(lapply(setups[,null_cols,with=FALSE],test_stat,group=setups$group,data=X))
# stat_names = colnames(nulls)
# nulls = transpose(nulls)
# 
# # indices = setups$obs
# # group = setups$group
# # data=X
# observed = as.numeric(test_stat(setups$obs,setups$group,data=X))
# names(observed) = stat_names
# ps = rowMeans(abs(nulls)>=abs(observed))
# names(ps) <- stat_names
# table(ps)
# ps = c(rn="p_any_vs_conts",ps)
# if(all(names(ps)==colnames(out))){
#   out = rbind(out,t(ps))
# }
# 
# fwrite(out,"~/work/ibs/clinical/clinical_table.csv")


#Process matrix func
#Calls the group processing for each col
#Considers all stats but the first to be 

#Process group setup func, takes one col of this matrix
#Gets all stats for both groups
#Calculates all differences
#Returns differences



###

dists = rbindlist(lapply(subset_cols,dists_for_col,data=X))

p1=ggplot(dists[variable=="affected.work.weeks"]) + stat_count(mapping = aes(x=value, y=..prop.., group=group, fill=group)) + facet_wrap(~group,ncol=1) + theme(legend.position = "none") + ggtitle("Affected at work weeks") + background_grid()
p2=ggplot(dists[variable=="absent.work.weeks"]) + stat_count(mapping = aes(x=value, y=..prop.., group=group, fill=group)) + facet_wrap(~group,ncol=1) + theme(legend.position = "none") + ggtitle("Absent from work weeks") + background_grid()
p3=ggplot(dists[variable=="Q_SSS"]) + stat_count(mapping = aes(x=value, y=..prop.., group=group, fill=group)) + facet_wrap(~group,ncol=1) + theme(legend.position = "none") + ggtitle("SSS") + background_grid()
p4=ggplot(dists[variable=="Q_phq_sum"]) + stat_count(mapping = aes(x=value, y=..prop.., group=group, fill=group)) + facet_wrap(~group,ncol=1) + theme(legend.position = "none") + ggtitle("PHQ") + background_grid()


plot_grid(p1,p2,p3,ncol=3)


X_ = melt(X,id.vars=c("Q_phq_sum","Q_SSS","eid"),measure.vars=subset_cols)
X_ = X_[value==1][,value:=NULL]
X_ = X_[!is.na(Q_phq_sum) & !is.na(Q_SSS)]
ggplot(X_[,.SD[sample(.N,min(.N,1e3))],by=variable], aes(x=Q_phq_sum, y=Q_SSS)) + 
  geom_point(aes(col=variable),alpha=0.2) + 
  geom_smooth() + 
  facet_wrap(~variable,nrow=3) + labs(col="group",x="PHQ-12",y="SSS") +
  ggtitle("PHQ-12 and SSS correlation\nN=1000 per group")# +
  #theme(legend.position="none")


#X_ = X_[,subset:=sample(rep(0,,by=variable]

##
X[,grp:=final_cse_]
ggplot(X[final_case.any==TRUE]) + geom_violin(aes(x=as.factor(round(Q_SSS/10)*10),y=`189-0.0`))
ggplot(X[final_cc.any%in%c(0:1)]) + geom_violin(aes(x=cut_interval(Q_SSS,5),y=`189-0.0`)) + facet_wrap(~final_cc.any)

ggplot(X[final_cc.any==1]) + geom_hex(aes(x=Q_SSS,y=`189-0.0`))


#SSS distribution with medians
median_grps = c("final_cc.any","final_cc.diag","final_cc.sr","final_cc.Q","final_cc.Q.prev.diag.yes")
medians = out[group%in%median_grps,.(group,description,median_SSS=as.numeric(median_SSS))]
medians[group=="final_cc.Q.prev.diag.yes",description:="Quest. Prev. Diag."]
medians[median_SSS==19,median_SSS:=median_SSS+c(-0.25,0.25)]
medians[group=="final_cc.Q",description:="Quest. Symp."]
ggplot(X[final_cc.any==1 & Q_all_SSS==TRUE]) + geom_histogram(aes(x=Q_SSS),fill="lightgrey",binwidth=1) +
  geom_vline(aes(xintercept=as.numeric(median_SSS),col=description),cex=2,data=medians) +
  geom_point(x=17.5,y=0,shape=6,col="black",fill="black",cex=2) +
  #geom_text(x=17.5,y=0,label="ok") + 
  geom_point(x=30,y=0,shape=6,col="black",fill="black",cex=2) +
  labs(col="Median SSS") +
  xlab("SSS") + ylab("No. of pooled cases") + ggtitle("SSS distribution among pooled IBS cases")


#SSS distribution without medians
S_ = X[final_cc.any==1 & Q_all_SSS==TRUE,.N,by=Q_SSS]
S_[,Q_SSS:=10*Q_SSS]
S_[Q_SSS<175,severity:="mild"]
S_[Q_SSS>=175 & Q_SSS<=300,severity:="moderate"]
S_[Q_SSS>300,severity:="severe"]
library(cowplot)
theme_set(theme_cowplot())
p.SSSdist = ggplot(S_) + 
  geom_bar(aes(x=Q_SSS,y=N,fill=severity),stat="identity") + 
  scale_fill_manual(values=c("severe"="tomato","moderate"="darkorange","mild"="goldenrod")) +
  ylab("No. of pooled cases") + 
  xlab("IBS-SSS") +
  labs(fill="IBS severity")
ggsave("~/work/ibs/plots/SSS_dist.pdf",p.SSSdist,width=7,height=5)

#function(subset_col, data){
#}
#male/female
# C[Q_participant==TRUE, lapply(.SD,function(case_col){sum(case_col==1,na.rm=TRUE)}),.SDcols=cols]
# 
# C[Q_participant==TRUE, lapply(.SD,function(case_col, ){sum(case_col==1,na.rm=TRUE)}),.SDcols=cols]
# 
# a = sapply(cols, function(col, data){data[get(col)==1, sum(sex.sample==2,na.rm=TRUE)/.N]}, data=C)
# b = sapply(cols, function(col, data){data[get(col)==1, sum(sex.sample==1,na.rm=TRUE)/.N]}, data=C)
# rowSums(data.table(a,b))

#Stats for in manuscript text, Wilcoxon tests and GAD-7 and PHQ-9 respondents
#See also https://www.notion.so/cqe/Fisher-tests-for-associated-conditions-and-risk-factors-of-IBS-8a1d86e625f74abdb807f03dca0a80e8
X[final_cc.any%in%c(0,1),sum(!is.na(GAD7))]
#79430
X[final_cc.any%in%c(0,1),sum(!is.na(PHQ9))]
#79087

wilcox.test(X[final_case.any==TRUE,as.numeric(GAD7)], X[final_case.any==FALSE,as.numeric(GAD7)])

# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  X[final_case.any == TRUE, as.numeric(GAD7)] and X[final_case.any == FALSE, as.numeric(GAD7)]
# W = 1873546930, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0


wilcox.test(X[final_case.any==TRUE,as.numeric(PHQ9)], X[final_case.any==FALSE,as.numeric(PHQ9)])
# Wilcoxon rank sum test with continuity correction
# 
# data:  X[final_case.any == TRUE, as.numeric(PHQ9)] and X[final_case.any == FALSE, as.numeric(PHQ9)]
# W = 1890849124, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0


X[final_case.any%in%c(TRUE,FALSE) & Q_all_SSS==TRUE & all_GAD7==TRUE,cor.test(Q_SSS,GAD7)]
# 
# Pearson's product-moment correlation
# 
# data:  Q_SSS and GAD7
# t = 98.526, df = 125473, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2628308 0.2731024
# sample estimates:
#       cor 
# 0.2679742
X[final_case.any%in%c(TRUE,FALSE) & Q_all_SSS==TRUE & all_PHQ9==TRUE,cor.test(Q_SSS,PHQ9)]
# data:  Q_SSS and PHQ9
# t = 111.42, df = 124913, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.2956257 0.3057140
# sample estimates:
#   cor 
# 0.3006782 


X[final_case.any%in%c(TRUE) & Q_all_SSS==TRUE & all_GAD7==TRUE,cor.test(Q_SSS,GAD7)]
#Pearson's product-moment correlation
# 
# data:  Q_SSS and GAD7
# t = 37.63, df = 23825, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2248330 0.2488036
# sample estimates:
#       cor 
# 0.2368543
X[final_case.any%in%c(TRUE) & Q_all_SSS==TRUE & all_PHQ9==TRUE,cor.test(Q_SSS,PHQ9)]
# 
# Pearson's product-moment correlation
# 
# data:  Q_SSS and PHQ9
# t = 42.322, df = 23643, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2534822 0.2771798
# sample estimates:
#       cor 
# 0.2653711  


wilcox.test(X[final_case.Q==TRUE & Q_all_SSS==TRUE,Q_SSS],
            X[final_case.sr | final_case.Q.prev.diag.yes==TRUE | final_case.diag==TRUE & Q_all_SSS==TRUE,Q_SSS])
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  X[final_case.Q == TRUE & Q_all_SSS == TRUE, Q_SSS] and X[final_case.sr | final_case.Q.prev.diag.yes == TRUE | final_case.diag == X[final_case.Q == TRUE & Q_all_SSS == TRUE, Q_SSS] and     TRUE & Q_all_SSS == TRUE, Q_SSS]
# W = 229717462, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0



wilcox.test(X[final_case.any==TRUE,as.numeric(Q_childhood.antibiotics)], X[final_case.any==FALSE,as.numeric(Q_childhood.antibiotics)])

# Wilcoxon rank sum test with continuity correction
# 
# data:  X[final_case.any == TRUE, as.numeric(Q_childhood.antibiotics)] and X[final_case.any == FALSE, as.numeric(Q_childhood.antibiotics)]
# W = 1.0236e+10, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(X[final_case.any==TRUE,as.numeric(Q_family.history)], X[final_case.any==FALSE,as.numeric(Q_family.history)])
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  X[final_case.any == TRUE, as.numeric(Q_family.history)] and X[final_case.any == FALSE, as.numeric(Q_family.history)]
# W = 1.0527e+10, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0


#Family history of IBS between cases and controls...
#... with the pooled IBS cases
table(X$Q_family.history, X$final_cc.any)
fisher.test(X$Q_family.history, X$final_cc.any)


#Childhood antibiotics exposure of IBS between cases and controls...
#... with the pooled IBS cases
table(X$Q_childhood.antibiotics, X$final_cc.any)
fisher.test(X$Q_childhood.antibiotics, X$final_cc.any)


#Anxiety treatment sought/offered between cases and controls...
#... with the pooled IBS cases
table(X$Q_treatment.anx, X$final_cc.any)
fisher.test(X$Q_treatment.anx,X$final_cc.any)
#... in functional D
table(X$Q_treatment.anx, X$final_cc.func.D)
fisher.test(X$Q_treatment.anx,X$final_cc.func.D)
#... in IBS-D
table(X$Q_treatment.anx, X$final_cc.sub.D)
fisher.test(X$Q_treatment.anx, X$final_cc.sub.D)


#Depression treatment sought/offered between cases and controls...
#... with the pooled IBS cases
table(X$Q_treatment.depr, X$final_cc.any)
fisher.test(X$Q_treatment.depr,X$final_cc.any)
#... in functional D
table(X$Q_treatment.depr, X$final_cc.func.D)
fisher.test(X$Q_treatment.depr,X$final_cc.func.D)
#... in IBS-D
table(X$Q_treatment.depr, X$final_cc.sub.D)
fisher.test(X$Q_treatment.depr, X$final_cc.sub.D)


glm(data=X,Q_treatment.anx~final_cc.func.D+sex+age+agesq+sexage+sexagesq,family="binomial")

#glm
X[final_cc.any==TRUE,table(final_cc.func.D)]
X[final_cc.any==TRUE,table(final_cc.func.D,useNA="always")]

X[final_cc.func.D==TRUE,table(final_cc.any==TRUE,useNA="always")]


coef_ci = function(m){
  summary = cbind(estimate=coef(summary(m)),confint.default(m))
  text = paste0("\nDoes ",rownames(summary)[2],
                " help predict ",as.character(formula(m))[2],
                "? OR and 95%CI: ",
                round(exp(summary[2,"Estimate"]),2),
                " [",
                round(exp(summary[2,"2.5 %"]),2),
                "-",
                round(exp(summary[2,"97.5 %"]),2),
                "], p=",
                summary[2,"Pr(>|z|)"],"\n")
  cat(text)
  return(summary)
}

excess_or = function(m1,m2){
  summary1 = cbind(estimate=coef(summary(m1)),confint.default(m1))
  summary2 = cbind(estimate=coef(summary(m2)),confint.default(m2))
  
  covar=2
  
  estimate_diff = summary1[covar,"Estimate"] - summary2[covar,"Estimate"] 
  se_pooled = sqrt(summary1[covar,"Std. Error"]^2 + summary2[covar,"Std. Error"]^2) #adjust for n1, n2?
  
  excess_or = exp(estimate_diff)
  ci_lower = exp(estimate_diff - 1.96*se_pooled)
  ci_upper = exp(estimate_diff + 1.96*se_pooled)
  text = c("excess OR:",
           round(excess_or,2),
           " [",
           round(ci_lower,2),
           "-",
           round(ci_upper,2),
           "]")
  
  cat(text)
}

########
#Depression
########
coef_ci(glm(data=X,final_cc.sub.D~Q_treatment.depr+sex+age+agesq+sexage+sexagesq,family="binomial"))
coef_ci(glm(data=X,final_cc.func.D~Q_treatment.depr+sex+age+agesq+sexage+sexagesq,family="binomial"))
coef_ci(glm(data=X,final_cc.sub.C~Q_treatment.depr+sex+age+agesq+sexage+sexagesq,family="binomial")) #not reported
coef_ci(glm(data=X,final_cc.func.C~Q_treatment.depr+sex+age+agesq+sexage+sexagesq,family="binomial")) #not reported

########
#Anxiety
########
coef_ci(glm(data=X,final_cc.sub.D~Q_treatment.anx+sex+age+agesq+sexage+sexagesq,family="binomial"))
coef_ci(glm(data=X,final_cc.func.D~Q_treatment.anx+sex+age+agesq+sexage+sexagesq,family="binomial"))
coef_ci(glm(data=X,final_cc.sub.C~Q_treatment.anx+sex+age+agesq+sexage+sexagesq,family="binomial")) #not reported
coef_ci(glm(data=X,final_cc.func.C~Q_treatment.anx+sex+age+agesq+sexage+sexagesq,family="binomial")) #not reported

######
#Combined depression/anxiety treatment
######
X[,Q_treatment.anx.or.depr:= Q_treatment.anx | Q_treatment.depr]
excess_or(glm(data=X,final_cc.sub.D~Q_treatment.anx.or.depr+sex+age+agesq+sexage+sexagesq,family="binomial"),
           glm(data=X,final_cc.func.D~Q_treatment.anx.or.depr+sex+age+agesq+sexage+sexagesq,family="binomial"))

excess_or(glm(data=X,final_cc.sub.C~Q_treatment.anx.or.depr+sex+age+agesq+sexage+sexagesq,family="binomial"),
           glm(data=X,final_cc.func.C~Q_treatment.anx.or.depr+sex+age+agesq+sexage+sexagesq,family="binomial"))






#individual
coef_ci(glm(data=X,final_cc.sub.D~Q_treatment.anx.or.depr+sex+age+agesq+sexage+sexagesq,family="binomial"))
coef_ci(glm(data=X,final_cc.func.D~Q_treatment.anx.or.depr+sex+age+agesq+sexage+sexagesq,family="binomial")) #and then the same with C

#direct
X[,func.C.vs.sub.C:=NULL]
X[final_cc.sub.C==1,func.C.vs.sub.C:=1]
X[final_cc.func.C==1,func.C.vs.sub.C:=0]
coef_ci(glm(data=X[Q_all_PHQ12==TRUE], func.C.vs.sub.C~Q_treatment.anx.or.depr+sex+age+agesq+sexage+sexagesq, family="binomial"))


X[,func.D.vs.sub.D:=NULL]
X[final_cc.sub.D==1,func.D.vs.sub.D:=1]
X[final_cc.func.D==1,func.D.vs.sub.D:=0]
coef_ci(glm(data=X[Q_all_PHQ12==TRUE], func.D.vs.sub.D~Q_treatment.anx.or.depr+sex+age+agesq+sexage+sexagesq, family="binomial"))



#######
#PHQ-12
#######
coef_ci(glm(data=X[Q_all_PHQ12==TRUE],final_cc.func.D~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq,family="binomial"))
coef_ci(glm(data=X[Q_all_PHQ12==TRUE],final_cc.sub.D~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq,family="binomial"))
excess_or(glm(data=X[Q_all_PHQ12==TRUE],final_cc.func.D~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq,family="binomial"),
           glm(data=X[Q_all_PHQ12==TRUE],final_cc.sub.D~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq,family="binomial"))

excess_or(glm(data=X[Q_all_PHQ12==TRUE],final_cc.func.C~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq,family="binomial"),
          glm(data=X[Q_all_PHQ12==TRUE],final_cc.sub.C~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq,family="binomial"))
           


#direct
coef_ci_t = function(m){
  summary = cbind(estimate=coef(summary(m)),confint.default(m))
  text = paste0("\nDoes ",rownames(summary)[2],
                " help predict ",as.character(formula(m))[2],
                "? OR and 95%CI: ",
                round(exp(summary[2,"Estimate"]),2),
                " [",
                round(exp(summary[2,"2.5 %"]),2),
                "-",
                round(exp(summary[2,"97.5 %"]),2),
                "], p=",
                summary[2,"Pr(>|t|)"],"\n")
  cat(text)
  return(summary)
}


X[,func.C.vs.sub.C:=NULL]
X[final_cc.sub.C==1,func.C.vs.sub.C:=1]
X[final_cc.func.C==1,func.C.vs.sub.C:=0]
coef_ci(glm(data=X[Q_all_PHQ12==TRUE], func.C.vs.sub.C~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq, family="binomial"))


X[,func.D.vs.sub.D:=NULL]
X[final_cc.sub.D==1,func.D.vs.sub.D:=1]
X[final_cc.func.D==1,func.D.vs.sub.D:=0]
coef_ci(glm(data=X[Q_all_PHQ12==TRUE], func.D.vs.sub.D~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq, family="binomial"))

#X[,func.C.vs.sub.C:=as.integer(func.C.vs.sub.C)]

#X[,func.D.vs.sub.D:=factor(func.D.vs.sub.D)]

typeof(X$func.C.vs.sub.C)
table(X$func.C.vs.sub.C)

typeof(X$final_cc.sub.C)
table(X$final_cc.sub.C)

coef_ci_t(glm(data=X[Q_all_PHQ12==TRUE], func.C.vs.sub.C~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq))


coef_ci(glm(data=X[Q_all_PHQ12==TRUE], func.C.vs.sub.C~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq))

coef_ci(glm(data=X[Q_all_PHQ12==TRUE],final_cc.sub.C~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq,family="binomial"))




X[final_cc.sub.D==1,func.D.vs.sub.D:=1]
X[final_cc.func.D==1,func.D.vs.sub.D:=0]
X[,func.D.vs.sub.D:=factor(func.D.vs.sub.D)]

coef_ci(glm(data=X[Q_all_PHQ12==TRUE], func.D.vs.sub.D~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq))
coef_ci(glm(data=X[Q_all_PHQ12==TRUE], func.C.vs.sub.C~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq))

coef_ci(glm(data=X[Q_all_PHQ12==TRUE],final_cc.func.C~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq,family="binomial"))
coef_ci(glm(data=X[Q_all_PHQ12==TRUE],final_cc.sub.C~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq,family="binomial"))
compare_or(glm(data=X[Q_all_PHQ12==TRUE],final_cc.func.C~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq,family="binomial"),
           glm(data=X[Q_all_PHQ12==TRUE],final_cc.sub.C~Q_PHQ12_sum+sex+age+agesq+sexage+sexagesq,family="binomial"))


#GAD-7
coef_ci(glm(data=X[all_GAD7==TRUE],final_cc.any~GAD7+sex+age+agesq+sexage+sexagesq+Q_participant,family="binomial")) #Q_participant added at NEJM review stage

#PHQ-9
coef_ci(glm(data=X[all_PHQ9==TRUE],final_cc.any~PHQ9+sex+age+agesq+sexage+sexagesq+Q_participant,family="binomial")) #Q_participant added at NEJM review stage

#Formal diagnosis with IBS-SSS

#Define informal diagnosis (0, Rome III and no other diagnoses) vs formal diagnosis (1, Rome III and any of ICD-10, UKB unprompted self-rep, DHQ prompted self-rep.)
X[,final_cc.formality:=NULL]
X[final_cc.Q==1, final_cc.formality:=0] #order is important
table(X$final_cc.formality)
X[final_cc.Q==1 & (final_cc.Q.prev.diag.yes==1 | final_cc.diag==1 | final_cc.sr==1),final_cc.formality:=1] #order is important
table(X$final_cc.formality)

coef_ci(glm(data=X[Q_all_SSS==TRUE],final_cc.formality~Q_SSS+sex+age+agesq+sexage+sexagesq,family="binomial"))


#coef_ci(glm(data=X,Q_childhood.antibiotics~Q_SSS+sex+age+agesq+sexage+sexagesq,family="binomial"))


#Below tests are only among respondents, no need to control for DHQ respondent status
X[!is.na(Q_family.history),all(Q_participant)]
X[!is.na(Q_childhood.antibiotics),all(Q_participant)]

#Family history ~ IBS case/control 
coef_ci(glm(data=X,Q_family.history~final_cc.any+sex+age+agesq+sexage+sexagesq,family="binomial"))
#Caesarian section ~ IBS case/control
coef_ci(glm(data=X,Q_born.by.caesarean~final_cc.any+sex+age+agesq+sexage+sexagesq,family="binomial"))
#Childhood antibiotics ~ IBS case/control 
coef_ci(glm(data=X,Q_childhood.antibiotics~final_cc.any+sex+age+agesq+sexage+sexagesq,family="binomial"))

#Childhood antibiotics ~ IBS-SSS
coef_ci(glm(data=X,Q_childhood.antibiotics~Q_SSS+sex+age+agesq+sexage+sexagesq,family="binomial"))
#Family history ~ IBS-SSS
coef_ci(glm(data=X,Q_family.history~Q_SSS+sex+age+agesq+sexage+sexagesq,family="binomial"))
#Childhood antibiotics ~ Anxiety case/control
coef_ci(glm(data=X,Q_childhood.antibiotics~final_cc.anx+sex+age+agesq+sexage+sexagesq,family="binomial"))



# glm(data=X[Q_all_SSS==TRUE],final_cc.formality~Q_SSS+sex+age+agesq+sexage+sexagesq,family="binomial")
# 
# wilcox.test(X[final_cc.formality==1 & Q_all_SSS==TRUE, Q_SSS],
#             X[final_cc.formality==0 & Q_all_SSS==TRUE, Q_SSS])


# 
# > coef_ci = function(m){
#   +   print(cbind(estimate=coef(summary(m)),confint.default(m)))
#   + }
# > coef_ci(glm(data=X,final_cc.sub.D~Q_treatment.depr+sex+age+agesq+sexage+sexagesq,family="binomial"))
# Estimate   Std. Error    z value      Pr(>|z|)         2.5 %        97.5 %
#   (Intercept)          -2.7924464151 0.0979895922 -28.497378 1.262444e-178 -2.984502e+00 -2.6003903434
# Q_treatment.deprTRUE  0.7001017393 0.0284788193  24.583243 1.908641e-133  6.442843e-01  0.7559191994
# sex                   3.5270531954 1.8929985706   1.863210  6.243278e-02 -1.831558e-01  7.2372622166
# age                   0.0838153098 0.0689639498   1.215350  2.242328e-01 -5.135155e-02  0.2189821676
# agesq                -0.0008734584 0.0007021321  -1.244009  2.134964e-01 -2.249612e-03  0.0005026952
# sexage               -0.0433077161 0.0283740880  -1.526312  1.269321e-01 -9.891991e-02  0.0123044744
# sexagesq              0.0001460063 0.0000922436   1.582834  1.134593e-01 -3.478782e-05  0.0003268004

exp(c(0.7001017393,6.442843e-01,0.7559191994))
#2.013958 1.904623 2.129568

# > coef_ci(glm(data=X,final_cc.sub.D~Q_treatment.anx+sex+age+agesq+sexage+sexagesq,family="binomial"))
# Estimate   Std. Error    z value      Pr(>|z|)         2.5 %        97.5 %
#   (Intercept)         -2.8205110202 9.813911e-02 -28.739930 1.210267e-181 -3.012860e+00 -2.6281619050
# Q_treatment.anxTRUE  0.8118775556 2.874730e-02  28.241872 1.791169e-175  7.555339e-01  0.8682212268
# sex                  3.7295537121 1.895241e+00   1.967852  4.908501e-02  1.495058e-02  7.4441568436
# age                  0.0909019918 6.905028e-02   1.316461  1.880194e-01 -4.443407e-02  0.2262380525
# agesq               -0.0008703208 7.025827e-04  -1.238745  2.154400e-01 -2.247358e-03  0.0005067160
# sexage              -0.0463462281 2.840794e-02  -1.631453  1.027947e-01 -1.020248e-01  0.0093323166
# sexagesq             0.0001541413 9.234262e-05   1.669233  9.507123e-02 -2.684687e-05  0.0003351296
exp(c(0.8118775556, 7.555339e-01,0.8682212268))
#2.252133 2.128748 2.382669

# > coef_ci(glm(data=X,final_cc.func.D~Q_treatment.depr+sex+age+agesq+sexage+sexagesq,family="binomial"))
# Estimate   Std. Error      z value      Pr(>|z|)         2.5 %        97.5 %
#   (Intercept)          -2.739372e+00 1.125344e-01 -24.34253919 6.954178e-131 -2.9599354576 -2.5188088580
# Q_treatment.deprTRUE  3.025606e-01 3.464943e-02   8.73205179  2.500937e-18  0.2346490057  0.3704722891
# sex                   6.972058e-02 2.084618e+00   0.03344525  9.733195e-01 -4.0160563385  4.1554974928
# age                   1.558848e-02 7.627629e-02   0.20436860  8.380655e-01 -0.1339103008  0.1650872573
# agesq                -7.664969e-04 6.356510e-04  -1.20584559  2.278770e-01 -0.0020123500  0.0004793561
# sexage               -7.899735e-03 3.123351e-02  -0.25292497  8.003262e-01 -0.0691162926  0.0533168227
# sexagesq              5.420294e-05 9.952316e-05   0.54462642  5.860105e-01 -0.0001408589  0.0002492648
exp(c(3.025606e-01,0.2346490057,0.3704722891))
#1.353320 1.264465 1.448419

# > coef_ci(glm(data=X,final_cc.func.D~Q_treatment.anx+sex+age+agesq+sexage+sexagesq,family="binomial"))
# Estimate   Std. Error      z value      Pr(>|z|)         2.5 %        97.5 %
#   (Intercept)         -2.743114e+00 1.125517e-01 -24.37203875 3.385855e-131 -2.9637110124 -2.5225165538
# Q_treatment.anxTRUE  3.323975e-01 3.563415e-02   9.32806085  1.078255e-20  0.2625558364  0.4022391205
# sex                  1.338485e-01 2.084879e+00   0.06419965  9.488113e-01 -3.9524386007  4.2201355804
# age                  1.768962e-02 7.628749e-02   0.23188099  8.166305e-01 -0.1318311172  0.1672103549
# agesq               -7.618536e-04 6.356546e-04  -1.19853391  2.307092e-01 -0.0020077137  0.0004840065
# sexage              -8.841107e-03 3.123757e-02  -0.28302803  7.771553e-01 -0.0700656097  0.0523833963
# sexagesq             5.669385e-05 9.953558e-05   0.56958378  5.689600e-01 -0.0001383923  0.0002517800
exp(c(3.323975e-01,0.2625558364,0.4022391205))
# [1] 1.394307 1.300249 1.495169

#Distinct from IBS-D, functional diarrhea was more common in men (60% vs 30% prevalence in men)
#Distinct from IBS-D, functional diarrhea was more common in men (30% vs 60% prevalence in men) but less strongly associated with anxiety (OR and 95%CI: 1.39 [1.30-1.50] vs 2.25 [2.13-2.38]) or depression (OR and 95%CI: 1.35 [1.26-1.45] vs 2.01 [1.90-2.13]).

