setwd("~/work/ibs/")
library(data.table)
library(ggplot2)
library(cowplot)

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
    Q_participation=100*mean(Q_participant,na.rm=TRUE),
    SSS_started = 100*mean(Q_oneplus_SSS,na.rm=TRUE),
    SSS_complete = 100*mean(Q_all_SSS,na.rm=TRUE))]
  
  part_two = data[subset & Q_all_SSS, .(median_SSS=median(Q_SSS, na.rm=TRUE))] #only amongst completed SSS
  
  extra_surg_cols = grep("surg.extra",colnames(data),value=TRUE)
  part_three = data[subset, lapply(.SD,function(col){100*sum(col==1,na.rm=TRUE)/.N}),.SDcols=extra_surg_cols]
  
  risk_cols = c("21065-0.0","21067-0.0","21066-0.0","21062-0.0","21063-0.0")
  risk_names = c("family.history","childhood.antibiotics","born.by.caesarean","treatment.anx","treatment.depr")
  part_four = data[subset & Q_participant, lapply(.SD,function(col){100*sum(col==1,na.rm=TRUE)/.N}),.SDcols=risk_cols] #only amongst participants
  names(part_four) = risk_names
  
  
  phq_cols = grep("^Q_phq_",colnames(data),value=TRUE)
  phq_cols = phq_cols[phq_cols!="Q_phq_sum"]
  
  phq.1 = data[subset,.(phq_started = 100*mean(Q_oneplus_phq,na.rm=TRUE), phq_complete = 100*mean(Q_all_phq,na.rm=TRUE))]
  phq.2 = data[subset & Q_all_phq, lapply(.SD,mean,na.rm=TRUE), .SDcols=phq_cols] #normally na.rm not needed, only for period question
  colnames(phq.2) <- gsub("^Q_","",colnames(phq.2))
  phq.3 = data[subset & Q_all_phq, .(median_phq=median(Q_phq_sum,na.rm=TRUE))]
  part_phq = cbind(phq.1,phq.2,phq.3)
  
  interest_cols = grep("interest",colnames(data),value=TRUE)
  part_immune = data[subset, lapply(.SD,function(col){100*sum(col==1,na.rm=TRUE)/.N}),.SDcols=interest_cols]
  
  townsend = data[subset,mean(`189-0.0`,na.rm=TRUE)]
  
  if(nrow(part_two)==0){part_two=data.table(median_SSS=NA)}
  result = cbind(part_one, part_two, part_four, part_phq, part_three,part_immune,townsend)
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
C = fread("/home/cq/work/ibs/covar/cc_cov.sample",key="FID")
P = fread("pheno/ukb24282.csv",key="eid")

A = fread("pheno/ukb21760.csv",key="eid")[,.(eid,`34-0.0`)]
A[,age_dhq:=2017-`34-0.0`]

E = fread("/home/cq/work/ibs/pheno/ukb39562_GAD7_PHQ9.csv",key="eid")
#X[final_cc.respondent==TRUE,tstrsplit(`21023-0.0`,"-")][,table(as.numeric(V1))], 2017: 153636 & 2018: 2935, only 1 percent completed outside of 2017
X = I[C,][P,][A,][E,]
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
X = X[Ptwo_mini,]
rm(Ptwo); rm(Ptwo_mini)
gc()


groups = fread("resource/IBS clinical call table - group_order.csv",header = TRUE)
subset_cols = groups$group

library(ggradar)
library(cowplot)

#Add the severity groups, mild, moderate, severe
#Thresholds of SSS (<175, 175-300, >300) = (mild, moderate, severe)
X[final_cc.any==1 & Q_all_SSS==TRUE & Q_SSS<17.5, final_case.SSS.mild:=1]
X[final_cc.any==1 & Q_all_SSS==TRUE & Q_SSS>=17.5 & Q_SSS<=30.0, final_case.SSS.moderate:=1]
X[final_cc.any==1 & Q_all_SSS==TRUE & Q_SSS>30.0, final_case.SSS.severe:=1]

sss_cols = paste0("final_case.SSS.",c("severe","moderate","mild"))
subset_cols = unique(c(subset_cols,sss_cols))


phq9_responses_for_col = function(col, data){
  if(is.numeric(col)){ #row numbers directly
    subset = logical(nrow(data))
    subset[col] = TRUE
  }else if(is.character(col)){ #rows where column==1
    subset = data[[col]]==1
  }else{
    quit("error")
  }
  
  phq_cols = grep("^PHQ-9.",colnames(data),value=TRUE)
  
  P = melt(data[subset & all_PHQ9],id.vars = "eid",measure.vars = phq_cols,variable.name = "question",value.name="response")
  P[,total_responses:=.N,by=question]
  P = P[,.(fraction=.N/total_responses[1], N=.N,total_responses=unique(total_responses)),by=c("question","response")]
  P[,group:=col]
  return(P)
}


gad7_responses_for_col = function(col, data){
  if(is.numeric(col)){ #row numbers directly
    subset = logical(nrow(data))
    subset[col] = TRUE
  }else if(is.character(col)){ #rows where column==1
    subset = data[[col]]==1
  }else{
    quit("error")
  }
  
  phq_cols = grep("^GAD-7.",colnames(data),value=TRUE)
  
  P = melt(data[subset & all_GAD7],id.vars = "eid",measure.vars = phq_cols,variable.name = "question",value.name="response")
  P[,total_responses:=.N,by=question]
  P = P[,.(fraction=.N/total_responses[1], N=.N,total_responses=unique(total_responses)),by=c("question","response")]
  P[,group:=col]
  return(P)
}

phq_responses_for_col = function(col, data){ #PHQ-12
  if(is.numeric(col)){ #row numbers directly
    subset = logical(nrow(data))
    subset[col] = TRUE
  }else if(is.character(col)){ #rows where column==1
    subset = data[[col]]==1
  }else{
    quit("error")
  }
  
  phq_cols = grep("^Q_PHQ-12",colnames(data),value=TRUE)
  phq_cols = phq_cols[phq_cols!="Q_PHQ12_sum"]
  
  P = melt(data[subset & Q_all_phq],id.vars = "eid",measure.vars = phq_cols,variable.name = "question",value.name="response")
  P[,total_responses:=.N,by=question]
  P = P[,.(fraction=.N/total_responses[1], N=.N,total_responses=unique(total_responses)),by=c("question","response")]
  P[,group:=col]
  return(P)
}




names = fread("/home/cq/work/ibs/resource/analysis_names.csv")
names = names[,.(group=nice_name, chique_name)]


radar_grps= c("final_cc.any",sss_cols,"final_cont")#"final_cc.sub.C","final_cc.sub.D","final_cc.sub.M","final_cc.sub.U") #"final_func.C","final_func.D"
color_grps = c("darkorchid","tomato","darkorange","goldenrod","limegreen")

#Mean PHQ-9 score
Pr = rbindlist(lapply(subset_cols,phq9_responses_for_col,data=X))

q = Pr[,.(mean_response=sum(response*fraction)),by=c("group","question")]
RD = dcast(q,group~question,value.var="mean_response",fun.aggregate = sum,fill=0)
RD = merge(RD,names,by="group",all.x=TRUE)
setcolorder(RD,c("group","chique_name"))
RD[group=="final_cont",chique_name:="UKBB: Control"]
RD[group=="final_cc.anx",chique_name:="UKBB: Anxiety"]
colnames(RD) <- gsub("^PHQ-9\\.","",colnames(RD))
colnames(RD) <- gsub("__","\n",colnames(RD))
colnames(RD) <- gsub("_"," ",colnames(RD))
colnames(RD)[2] <- "chique_name"
RD = RD[match(radar_grps,group)]
RD[,`chique_name`:=gsub("UKBB: ","",`chique_name`)]


ggradar(RD[,2:ncol(RD)],grid.mid=1,grid.max=2,centre.y=0,values.radar = c("0","1","2"),legend.position="bottom",group.colours=color_grps)
RD[,chique_name:=factor(chique_name,levels=unique(chique_name))] #crashes if ggradar is never run without?
ggradar(RD[,2:ncol(RD)],grid.min=0,grid.mid=1,grid.max=1.52,centre.y=0,values.radar = c("0","1","1.5"),legend.position="bottom",group.colours=color_grps,legend.text.size = 11)
ggsave("~/work/ibs/plots/phq9_radar_new_meanscore_croppedy.pdf",width=7.2,height=7.2)
ggsave("~/work/ibs/plots/phq9_radar_new_meanscore_croppedy.svg",width=7.2,height=7.2)

#Mean GAD-7 score
Pr = rbindlist(lapply(subset_cols,gad7_responses_for_col,data=X))

q = Pr[,.(mean_response=sum(response*fraction)),by=c("group","question")]
RD = dcast(q,group~question,value.var="mean_response",fun.aggregate = sum,fill=0)
RD = merge(RD,names,by="group",all.x=TRUE)
setcolorder(RD,c("group","chique_name"))
RD[group=="final_cont",chique_name:="UKBB: Control"]
RD[group=="final_cc.anx",chique_name:="UKBB: Anxiety"]
colnames(RD) <- gsub("^GAD-7\\.","",colnames(RD))
colnames(RD) <- gsub("__","\n",colnames(RD))
colnames(RD) <- gsub("_"," ",colnames(RD))
colnames(RD)[2] <- "chique_name"
RD = RD[match(radar_grps,group)]
RD[,`chique_name`:=gsub("UKBB: ","",`chique_name`)]

ggradar(RD[,2:ncol(RD)],grid.mid=1,grid.max=2,centre.y=0,values.radar = c("0","1","2"),legend.position="bottom",group.colours=color_grps)
RD[,chique_name:=factor(chique_name,levels=unique(chique_name))] #crashes if ggradar is never run without?
ggradar(RD[,2:ncol(RD)],grid.min=0,grid.mid=0.5,grid.max=1,centre.y=0,values.radar = c("0","0.5","1"),legend.position="bottom",group.colours=color_grps,legend.text.size = 11)
ggsave("~/work/ibs/plots/gad7_radar_new_meanscore_croppedy.pdf",width=7.2,height=7.2)
ggsave("~/work/ibs/plots/gad7_radar_new_meanscore_croppedy.svg",width=7.2,height=7.2)

#Mean PHQ-12 score
Pr = rbindlist(lapply(subset_cols,phq_responses_for_col,data=X))

q = Pr[,.(mean_response=sum(response*fraction)),by=c("group","question")]
RD = dcast(q,group~question,value.var="mean_response",fun.aggregate = sum,fill=0)
RD = merge(RD,names,by="group",all.x=TRUE)
setcolorder(RD,c("group","chique_name"))
RD[group=="final_cont",chique_name:="UKBB: Control"]
RD[group=="final_cc.anx",chique_name:="UKBB: Anxiety"]
colnames(RD) <- gsub("^Q_PHQ-12\\.","",colnames(RD))
colnames(RD) <- gsub("__","\n",colnames(RD))
colnames(RD) <- gsub("_"," ",colnames(RD))
colnames(RD)[2] <- "chique_name"
RD = RD[match(radar_grps,group)]
RD[,`chique_name`:=gsub("UKBB: ","",`chique_name`)]
RD = RD[,.SD,.SDcols=setdiff(colnames(RD),c("period","nausea"))]

ggradar(RD[,2:ncol(RD)],grid.mid=1,grid.max=2,centre.y=0,values.radar = c("0","1","2"),legend.position="bottom",group.colours=color_grps)
RD[,chique_name:=factor(chique_name,levels=unique(chique_name))] #crashes if ggradar is never run without?
ggradar(RD[,2:ncol(RD)],grid.min=0,grid.mid=1,grid.max=1.52,centre.y=0,values.radar = c("0","1","1.5"),legend.position="bottom",group.colours=color_grps,legend.text.size = 11)
ggsave("~/work/ibs/plots/phq12_radar_new_meanscore_croppedy.pdf",width=7.2,height=7.2)
ggsave("~/work/ibs/plots/phq12_radar_new_meanscore_croppedy.svg",width=7.2,height=7.2)
