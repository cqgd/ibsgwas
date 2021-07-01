#Created to satisfy NEJM reviewer comment
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(data.table)
library(mvtnorm)
library(numDeriv)

setwd("~/work/ibs")


#Liability model code
getProbs <- function(T1,T2,rho){
  #function for getting the four probabilities as a function of the liability parameters
  p00 <- pmvnorm(lower=c(-Inf,-Inf),upper=c(T1,T2),mean=c(0,0),sigma=cbind(c(1,rho),c(rho,1)))
  p01 <- pmvnorm(lower=c(-Inf,T2),upper=c(T1,Inf),mean=c(0,0),sigma=cbind(c(1,rho),c(rho,1)))
  p10 <- pmvnorm(lower=c(T1,-Inf),upper=c(Inf,T2),mean=c(0,0),sigma=cbind(c(1,rho),c(rho,1)))
  p11 <- pmvnorm(lower=c(T1,T2),upper=c(Inf,Inf),mean=c(0,0),sigma=cbind(c(1,rho),c(rho,1)))
  c(p00,p01,p10,p11)
}

getL <- function(par,data){
  #function to calculate the log likelihood
  ps <- getProbs(par[1],par[2],par[3])
  temp <- -sum(data*log(ps))
  temp
}

getL_cont <- function(par,x,d){
  #calculate the negative log likelihood
  T1 <- par[1]
  rho <- par[2]
  pd <- 1 - pnorm(T1,mean=rho*x,sd=sqrt(1 - rho^2))
  -sum(d*log(pd)) -sum((1 - d)*log(1 - pd))
}
#

liability_cont = function(d,x){
  initpar <- c(qnorm(1 - mean(d)),cor(d,x))
  fitmd <- optim(initpar,getL_cont,x=x,d=d,method="L-BFGS-B",lower=c(-Inf,-1 + 1e-8),upper=c(Inf,1 - 1e-8))
  est_par <- fitmd$par
  SE_par <- sqrt(diag(solve(hessian(function(p) getL_cont(p,x,d),fitmd$par))))
  mle = c(est_par[1],est_par[2],SE_par[1],SE_par[2])
  names(mle) <- c("mle.Z1t","mle.rho","mle.Z1t_SE","mle.rho_SE")
  return(mle)
}



contingency_table = function(col1,col2,data){
  data = rev(as.numeric(data[,table((as.logical(get(col1))), (as.logical(get(col2))))]))
  names(data) <- paste0("N_",c("both","just_disorder","just_IBS","neither"))
  return(data)
}

extract_counts = function(disorder, P, codes){
  codes = codes[comment==disorder]
  all_corr_cols = intersect(colnames(P),c(paste0("diag.correlation_",codes$ICD10.code),
                                          paste0("sr.correlation_",codes$SR.code)))
  if(length(all_corr_cols)>=1){
    all_corr_cols = all_corr_cols[!is.na(all_corr_cols)]
    corr_binary_col = paste0("corr_binary.",disorder)
    pooled_pheno = as.numeric(rowSums(P[,.SD,.SDcols=all_corr_cols])>=1)
    #pooled_pheno_Qnon[
    #pooled_pheno_Qonly
    set(P,i=NULL,j=corr_binary_col,pooled_pheno)
    set(P,i=which(P$Q_participant==FALSE),j=paste0(corr_binary_col,".Qnon"),pooled_pheno[which(P$Q_participant==FALSE)]) #Used for glm
    set(P,i=which(P$Q_participant==TRUE),j=paste0(corr_binary_col,".Qonly"),pooled_pheno[which(P$Q_participant==TRUE)]) #Used for glm
    set(P,i=which(P$sex<0 & P$Q_participant==FALSE),j=paste0(corr_binary_col,".Qnon.male"),pooled_pheno[which(P$sex<0 & P$Q_participant==FALSE)])
    set(P,i=which(P$sex<0 & P$Q_participant==TRUE),j=paste0(corr_binary_col,".Qonly.male"),pooled_pheno[which(P$sex<0 & P$Q_participant==TRUE)])
    set(P,i=which(P$sex>0 & P$Q_participant==FALSE),j=paste0(corr_binary_col,".Qnon.female"),pooled_pheno[which(P$sex>0 & P$Q_participant==FALSE)])
    set(P,i=which(P$sex>0 & P$Q_participant==TRUE),j=paste0(corr_binary_col,".Qonly.female"),pooled_pheno[which(P$sex>0 & P$Q_participant==TRUE)])
    
    
    counts_Qnon_male = contingency_table("final_cc.Qnon.any.male",paste0(corr_binary_col,".Qnon.male"),P)
    counts_Qonly_male = contingency_table("final_cc.Qonly.any.male",paste0(corr_binary_col,".Qonly.male"),P)
    counts_Qnon_female = contingency_table("final_cc.Qnon.any.female",paste0(corr_binary_col,".Qnon.female"),P)
    counts_Qonly_female = contingency_table("final_cc.Qonly.any.female",paste0(corr_binary_col,".Qonly.female"),P)
    
    data = data.table(rbind(counts_Qnon_male,counts_Qonly_male,counts_Qnon_female,counts_Qonly_female))
    data[,disorder:=disorder][,DHQ:=c("male non-respondents","male respondents","female non-respondents","female respondents")]
    #data = rev(as.numeric(P[,table((as.logical(final_cc.any)), (as.logical(get(corr_binary_col))))])) #disease counts, no NAtoF
  }else{
    data = NULL
  }
  return(data)
}

nice_glm_summary = function(m){
  OR = coef(m)[2]
  OR_ci = confint.default(m)[2,]
  OR_SE = unname(0.5*diff(OR_ci)/abs(qnorm(0.05/2))) #==sqrt(diag(vcov(m.Qonly)))[2]
  OR_p = 2*(1-pnorm(abs(OR/OR_SE)))
  OR_summ = c(OR,OR_SE,OR_ci,OR_p)
  names(OR_summ) = c("OR","OR_SE","OR_ci_lower","OR_ci_upper","OR_p")
  return(OR_summ)
}

OR_from_glm = function(disorder, P){ #Call this as per !!!! below
  corr_binary_col = paste0("corr_binary.",disorder)
  if(corr_binary_col%in%colnames(P)){
    m.Qnon = glm(final_cc.Qnon.any~get(paste0(corr_binary_col,".Qnon"))+sex+age+agesq+sexage+sexagesq,data=P,family="binomial")
    m.Qonly = glm(final_cc.Qonly.any~get(paste0(corr_binary_col,".Qonly"))+sex+age+agesq+sexage+sexagesq,data=P,family="binomial")
    
    data = data.table(rbind(nice_glm_summary(m.Qnon), nice_glm_summary(m.Qonly)))
    data[,disorder:=disorder][,DHQ:=c("non-respondents","respondents")]
  }else{
    warning(paste0("Disorder code ",disorder," has no column in P, call extract_counts() first"))
    data = NULL
  }
  return(data)
}


fish = function(ct){
  ct = as.numeric(ct)
  #if(any(ct[1:4]==0)){return(rep(NA,6))}
  res = fisher.test(cbind(ct[c(1,3)],ct[c(2,4)]))
  
  p = ct[1]/(ct[1]+ct[3])
  odds_disorder_case_given_IBS_case = p/(1-p)
  p = ct[2]/(ct[2]+ct[4])
  odds_disorder_case_given_IBS_cont = p/(1-p)
  summ = c(odds_disorder_case_given_IBS_case=odds_disorder_case_given_IBS_case,  
           odds_disorder_case_given_IBS_cont=odds_disorder_case_given_IBS_cont, 
           OR = log(unname(res$estimate)), 
           OR_SE = 0.5*diff(log(res$conf.int))/abs(qnorm(0.05/2)),
           OR_ci_lower = log(res$conf.int[1]), OR_ci_upper=log(res$conf.int[2]), OR_p=res$p.value)
  #this is all on a log scale, and meta-analysis will also be on a log scale
  return(summ)
}

liability = function(ct){
  data = as.numeric(ct)
  initRho <- (data[4]*data[1] - data[2]*data[3])/sqrt((data[1]+data[3])*(data[2] + data[4])*(data[1] + data[2])*(data[3] + data[4]))
  initpar <- c(qnorm(1 - (data[2] + data[4])/sum(data)),qnorm(1 - (data[3] + data[4])/sum(data)),initRho)
  
  fitmd <- try(optim(initpar,getL,data=data,method="L-BFGS-B",lower=c(-Inf,-Inf,-1 + 1e-8),upper=c(Inf,Inf,1 - 1e-8)))
  if(!class(fitmd)=="list"){
    mle = rep(NA,6)
  }else{
    est_par <- fitmd$par
    SE_par <- sqrt(diag(solve(hessian(getL,fitmd$par,data=data))))
    mle = c(est_par,SE_par)
  }
  
  names(mle) = paste0("mle.",c("Z1t","Z2t","rho","Z1t_SE","Z2t_SE","rho_SE"))
  return(mle)
}


inverse_se_meta = function(betas,ses,prefix=""){
  w = 1/(ses^2)
  se = sqrt(1/rowSums(w))
  beta = rowSums(betas*w)/rowSums(w)
  
  Z = beta/se
  p = pnorm(-abs(Z))
  
  meta = data.table(est=beta,se,Z,p)
  setNames(meta,paste0(prefix,toupper(colnames(meta)),".meta"))
}


#Format genetic correlations
eid_col = "eid"
diag_sr = fread("filters/diagnoses_selfreported.csv", key=eid_col)

rg_paths = c("/home/cq/work/ibs/ldhub/metal_Qonly_Qnon_any_ldhub.sumstats.8eaf0500-db2b-408a-bb30-2e2812c3547e.rg.results.csv",
             "/home/cq/work/ibs/ldhub/metal_ICD_ROME_EURUSA_Qonly_Qnon_ldhub.sumstats.f50840be-6a0b-4a89-a8cb-c12d716e35c6.rg.results.csv",
             "/home/cq/work/ibs/ldhub/metal_MAURO_ldhub.sumstats.b167734c-3cf0-4d0c-bf3a-ebfd7b476c19.rg.results.csv",
             "/home/cq/work/ibs/ldhub/subtypes/final_cc.sub.D_ldhub.sumstats.a6fa4985-ec34-4d3a-b189-67650c94e12d.rg.results.csv",
             "/home/cq/work/ibs/ldhub/subtypes/final_cc.sub.M_ldhub.sumstats.4addc1b2-2594-4d6c-9123-f823c9d3b69c.rg.results.csv",
             "/home/cq/work/ibs/ldhub/subtypes/final_cc.sub.C_ldhub.sumstats.7a618909-cfde-4e39-98de-a8df2456a4c5.rg.results.csv",
             "/home/cq/work/ibs/ldhub/defs/rg/final_cc.sr_ldhub.txt.rg.results.csv",
             "/home/cq/work/ibs/ldhub/defs/rg/final_cc.diag_ldhub.txt.rg.results.csv",
             "/home/cq/work/ibs/ldhub/defs/rg/final_cc.Q.prev.diag.yes_ldhub.txt.rg.results.csv",
             "/home/cq/work/ibs/ldhub/defs/rg/final_cc.Q_ldhub.txt.rg.results.csv",
             "/home/cq/work/ibs/ldhub/natgen/final_cc.Qonly.any.severe_ldhub.sumstats.40288636-bd20-4b65-b909-e49e8ca16343.rg.results.csv", #NatGen revision
             "/home/cq/work/ibs/ldhub/natgen/metal_Qonly_Qnon_any.twoplus_ldhub.sumstats.46ed25b5-d583-427d-bf07-faa3d1977a56.rg.results.csv")

R = lapply(rg_paths,fread)
R = rbindlist(R)
R[,ci_lower:=rg-1.96*se][,ci_upper:=rg+1.96*se]


med_rg = R[,.(median_rg=median(rg)),by=trait2]
R = merge(R,med_rg,by="trait2",all.x=TRUE)

setorderv(R,"median_rg",na.last = TRUE,order=-1)
R[,tname:=paste0(trait2,"_",Category,"_",PMID)]
R[,tname:=factor(tname,levels=rev(unique(tname)))]
R[,t2f:=factor(trait2,levels=unique(trait2))]
R[,worst_p:=max(p),by=trait2]
t1legend = c("metal_Qonly_Qnon"="UK Biobank data only",
             "metal_MAURO"="Bellygenes data only",
             "metal_big_pooled"="Discovery cohort",
             "sub.D"="IBS-D",
             "sub.C"="IBS-C",
             "sub.M"="IBS-M",
             "sr"="Hospital ICD-10", #Must be swapped, as data are from before labeling fix of Feb 2019.
             "diag"="Unprompted self-rep.", #Must be swapped with sr.
             "Q.prev.diag.yes"="DHQ self-rep.",
             "Q"="DHQ Rome III",
             "metal_Qonly_Qnon_any.twoplus"="Higher-specificity",
             "final_cc.Qonly.any.severe"="Severe (IBS-SSS>300)") #For transfer report
R[,trait1:=t1legend[trait1]]
best_defs = R[,tname[which.min(se)],by=trait2]$V1
#best_defs = R[,tname[which.min(p)],by=trait2]$V1 #
R = R[tname%in%best_defs]

#Agnostic selection:
imp = R[worst_p<0.05 & Category!="ukbb",unique(tname)]

#Let's add phenotype data
P = fread("/home/cq/work/ibs/filters/all_integrated.csv")
setnames(P,old="eid",new="FID")
C = fread("~/work/ibs/covar/cc_cov.sample")
P = merge(P,C[,.SD,.SDcol=c("FID",setdiff(colnames(C),colnames(P)))],by="FID",all.x=TRUE,sort=FALSE); rm(C); gc()
corr_cols = grep("correlation",colnames(P),value=TRUE)
P[,lapply(.SD,function(x){sum(x==1,na.rm=TRUE)}),.SDcols=corr_cols]
#Add sex-specific respondent and non-respondent cols
P[sex<0,final_cc.Qnon.any.male:=final_cc.Qnon.any]
P[sex>0,final_cc.Qnon.any.female:=final_cc.Qnon.any]
P[sex<0,final_cc.Qonly.any.male:=final_cc.Qonly.any]
P[sex>0,final_cc.Qonly.any.female:=final_cc.Qonly.any]

codes = fread("/home/cq/work/ibs/resource/IBS GWAS_ SR and ICD-10 codes of interest - diagnoses_selfreported.csv")
codes = codes[category=="correlation plot"]

#Get contingency tables for phenotypic correlations
cts = rbindlist(sapply(unique(codes$comment),extract_counts,P=P,codes=codes))

#Calculate OR and rho under liability model (happens for non-responds and respondents separately), 
OR_info = t(apply(cts[,1:4],1,fish)) #THESE ARE LOG(ODDS) DATA
LIABILITY_info = t(apply(cts[,1:4],1,liability))
ALL = cbind(cts,OR_info,LIABILITY_info)


#Add neuroticism
Qnon_neuro_male = P$sex<0 & is.finite(P$neuroticism.score) & P$final_cc.Qnon.any%in%c(0,1)
Qonly_neuro_male = P$sex<0 & is.finite(P$neuroticism.score) & P$final_cc.Qonly.any%in%c(0,1)
Qnon_neuro_female = P$sex>0 & is.finite(P$neuroticism.score) & P$final_cc.Qnon.any%in%c(0,1)
Qonly_neuro_female = P$sex>0 & is.finite(P$neuroticism.score) & P$final_cc.Qonly.any%in%c(0,1)

Qnon_neuro_mle_male = liability_cont(P[Qnon_neuro_male,final_cc.Qnon.any.male],P[Qnon_neuro_male,neuroticism.score])
Qonly_neuro_mle_male = liability_cont(P[Qonly_neuro_male,final_cc.Qonly.any.male],P[Qonly_neuro_male,neuroticism.score])
Qnon_neuro_mle_female = liability_cont(P[Qnon_neuro_female,final_cc.Qnon.any.female],P[Qnon_neuro_female,neuroticism.score])
Qonly_neuro_mle_female = liability_cont(P[Qonly_neuro_female,final_cc.Qonly.any.female],P[Qonly_neuro_female,neuroticism.score])

neuro_mle = data.table(rbind(Qnon_neuro_mle_male,Qonly_neuro_mle_male,Qnon_neuro_mle_female,Qonly_neuro_mle_female))
neuro_mle[,disorder:="Neuroticism"][,DHQ:=c("male non-respondents","male respondents","female non-respondents","female respondents")]

ALL = rbind(ALL,neuro_mle,use.names=TRUE,fill=TRUE)


#Meta-analyze 

setorder(ALL,disorder,DHQ)
if(!all(table(ALL$disorder)==4)){warning("Watch out with meta-analysis, respondent and non-respondent data needs to be in the same order of diseases")}
#fwrite(ALL,"~/work/ibs/tables/phenotypic_correlation_pre-meta-analysis.csv")

rho_betas = cbind(ALL[DHQ=="male respondents",mle.rho], 
                  ALL[DHQ=="male non-respondents",mle.rho],
                  ALL[DHQ=="female respondents",mle.rho],
                  ALL[DHQ=="female non-respondents",mle.rho])

rho_ses = cbind(ALL[DHQ=="male respondents",mle.rho_SE], 
                ALL[DHQ=="male non-respondents",mle.rho_SE],
                ALL[DHQ=="female respondents",mle.rho_SE],
                ALL[DHQ=="female non-respondents",mle.rho_SE])

rho_meta = inverse_se_meta(rho_betas, rho_ses, prefix="rho.")


OR_betas = cbind(ALL[DHQ=="male respondents",OR], 
                 ALL[DHQ=="male non-respondents",OR],
                 ALL[DHQ=="female respondents",OR],
                 ALL[DHQ=="female non-respondents",OR])

OR_ses = cbind(ALL[DHQ=="male respondents",OR_SE], 
               ALL[DHQ=="male non-respondents",OR_SE],
               ALL[DHQ=="female respondents",OR_SE],
               ALL[DHQ=="female non-respondents",OR_SE])

OR_meta = inverse_se_meta(OR_betas, OR_ses, prefix="OR.")

PC = cbind(disorder=unique(ALL$disorder),rho_meta,OR_meta)

#Add glm results
glms = rbindlist(sapply(unique(ALL$disorder),OR_from_glm,P=P),fill=TRUE)
glm_betas = cbind(glms[DHQ=="non-respondents",OR],
                  glms[DHQ=="respondents",OR])
glm_ses = cbind(glms[DHQ=="non-respondents",OR_SE],
                glms[DHQ=="respondents",OR_SE])
glm_meta = inverse_se_meta(glm_betas,glm_ses,prefix="glm.")
glm_meta[,disorder:=unique(glms$disorder)]
PC = merge(PC,glm_meta,by="disorder",all.x=TRUE,sort=FALSE)

#Merge phenotypic with genetic correlations
setnames(PC,old="disorder",new="trait2")
PC[,trait1:="UK Biobank data only"]

print("Phenotypic correlations with no matching genetic correlation:")
print(PC[!(PC$trait2%in%R$trait2),trait2])

R = merge(R,PC,by=c("trait2","trait1"),all.x=TRUE,sort=FALSE)
R[grepl("(Crohn)|(Ulcerative)",trait2),rho.EST.meta:=NA]

#Add confidence intervals (maybe do before merge)
R[,rho.ci_upper.meta:=rho.EST.meta+qnorm(1-(0.05/2))*rho.SE.meta][,rho.ci_lower.meta:=rho.EST.meta-qnorm(1-(0.05/2))*rho.SE.meta]
R[,rho.ci_lower.meta:=pmax(-1,rho.ci_lower.meta)][,rho.ci_upper.meta:=pmin(1, rho.ci_upper.meta)]

R[,OR.ci_upper.meta:=OR.EST.meta+qnorm(1-(0.05/2))*OR.SE.meta][,OR.ci_lower.meta:=OR.EST.meta-qnorm(1-(0.05/2))*OR.SE.meta]
R[,glm.ci_upper.meta:=glm.EST.meta+qnorm(1-(0.05/2))*glm.SE.meta][,glm.ci_lower.meta:=glm.EST.meta-qnorm(1-(0.05/2))*glm.SE.meta]
#R[,OR.ci_lower.meta:=pmax(-1,OR.ci_lower.meta)][,OR.ci_upper.meta:=pmin(1, OR.ci_upper.meta)]

#Fix order, rename traits, select traits

R[trait2=="Non-cancer illness code_ self-reported: anxiety/panic attacks",trait2:="Anxiety or panic attacks"]
R[trait2=="Doctor diagnosed hayfever or allergic rhinitis",trait2:="Hayfever"]
R[trait2=="Crohns disease",trait2:="Crohn's disease"]

manual_traits = c("Anxiety or panic attacks","Asthma","Ulcerative colitis","Crohn's disease")
auto_then_manual_order = c(setdiff(unique(R$trait2),manual_traits),manual_traits)
thematic_order = c("Anxiety or panic attacks",setdiff(unique(R$trait2),manual_traits),setdiff(manual_traits,"Anxiety or panic attacks"))

R[,t2f:=factor(trait2,levels=thematic_order)]
R[,trait1:=factor(trait1,levels=t1legend)]

selected_traits = unique(c(setdiff(codes[category=="correlation plot",comment],
                                   c("Subjective well being","Eczema","Crohns disease","Doctor diagnosed hayfever or allergic rhinitis","Excessive daytime sleepiness")),
                           "Anxiety or panic attacks","Crohn's disease"))

R[trait2%in%selected_traits,plot:=TRUE]

#Plotting pheno corr
p.rg = ggplot(R[plot==TRUE], aes(x=rg,y=trait1,col=trait1)) +
  geom_point() + 
  geom_errorbarh(aes(xmin=ci_lower, xmax=ci_upper)) + 
  geom_rect(data = R[trait2%in%manual_traits & trait1=="UK Biobank data only"],fill="goldenrod1",color=NA,xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.1, show.legend=FALSE) +
  #geom_text(aes(x=-1,label=p_label)) +
  xlim(-1.09,1.09) +
  geom_vline(xintercept=0,col="lightgrey") +
  background_grid() +
  ylab("Trait compared to IBS") + xlab("Genetic correlation") +
  labs(col="IBS definition") +
  #facet_grid(t2f~.,scales = "free",space="free_y") +
  facet_wrap(~t2f,strip.position="left",ncol=1) +
  theme(strip.text.y = element_text(angle = 180),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        legend.position="none") +
  scale_y_discrete(limits = rev(levels(R$trait1)))

#p.rg

p.mle = ggplot(R[plot==TRUE], aes(x=rho.EST.meta,y=trait1,col=trait1)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin=rho.ci_lower.meta, xmax=rho.ci_upper.meta)) + 
  geom_rect(data = R[trait2%in%manual_traits & trait1=="UK Biobank data only"],fill="goldenrod1",color=NA,xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.1, show.legend=FALSE) +
  #geom_text(aes(x=-1,label=p_label)) +
  #xlim(-1.09,1.09) +
  geom_vline(xintercept=0,col="lightgrey") +
  background_grid() +
  ylab("") + xlab("Phenotypic correlation") +
  labs(col="IBS definition") +
  #facet_grid(t2f~.,scales = "free",space="free_y") +
  facet_wrap(~t2f,strip.position="right",ncol=1) +
  theme(strip.text.y=element_blank(), #strip.text.y = element_text(angle = 0),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        legend.position="right") +
  scale_y_discrete(limits = rev(levels(R$trait1)))

#p.mle

p.OR = ggplot(R[plot==TRUE], aes(x=exp(OR.EST.meta),y=trait1,col=trait1)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin=exp(OR.ci_lower.meta), xmax=exp(OR.ci_upper.meta))) + 
  geom_rect(data = R[trait2%in%manual_traits & trait1=="UK Biobank data only"],fill="goldenrod1",color=NA,xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.1, show.legend=FALSE) +
  #geom_text(aes(x=-1,label=p_label)) +
  #xlim(-1.09,1.09) +
  geom_vline(xintercept=0,col="lightgrey") +
  background_grid() +
  ylab("") + xlab("OR") +
  labs(col="IBS definition") +
  #facet_grid(t2f~.,scales = "free",space="free_y") +
  facet_wrap(~t2f,strip.position="right",ncol=1) +
  theme(strip.text.y=element_blank(), #strip.text.y = element_text(angle = 0),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        legend.position="none") +
  scale_y_discrete(limits = rev(levels(R$trait1))) +
  scale_x_continuous(trans="log10")


#p.OR

p.glm = ggplot(R[plot==TRUE], aes(x=exp(glm.EST.meta),y=trait1,col=trait1)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin=exp(glm.ci_lower.meta), xmax=exp(glm.ci_upper.meta))) + 
  geom_rect(data = R[trait2%in%manual_traits & trait1=="UK Biobank data only"],fill="goldenrod1",color=NA,xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.1, show.legend=FALSE) +
  #geom_text(aes(x=-1,label=p_label)) +
  #xlim(-1.09,1.09) +
  geom_vline(xintercept=0,col="lightgrey") +
  background_grid() +
  ylab("") + xlab("OR") +
  labs(col="IBS definition") +
  #facet_grid(t2f~.,scales = "free",space="free_y") +
  facet_wrap(~t2f,strip.position="right",ncol=1) +
  theme(strip.text.y=element_blank(), #strip.text.y = element_text(angle = 0),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        legend.position="none") +
  scale_y_discrete(limits = rev(levels(R$trait1))) +
  scale_x_continuous(trans="log10",limits = c(1,NA))

#p.glm

corr_plot = plot_grid(p.rg,p.glm,p.mle,rel_widths=c(0.23,0.13,0.25),nrow=1)
print(corr_plot)
ggsave("~/work/ibs/plots/corr_plot_thematic_subtypes_and_defs.pdf",corr_plot,width=16,height=10)
#ggsave("~/work/ibs/plots/corr_plot_glm.pdf",corr_plot,width=16,height=6)
#ggsave("~/work/ibs/plots/corr_plot_new.pdf",corr_plot,width=16,height=6)
#ggsave("~/work/ibs/plots/corr_plot_log.pdf",corr_plot,width=16,height=6)
fwrite(R[plot==TRUE],"~/work/ibs/tables/genetic_phenotypic_correlation_plot_with_subtypes_and_defs.csv")
fwrite(R[,.(trait1, trait2,	rg,	se,	z,	p,	h2_obs,	h2_obs_se,	h2_int,	h2_int_se,	gcov_int,	gcov_int_se,	ci_lower,	ci_upper, PMID,	Category,	ethnicity, note)],
       "~/work/ibs/tables/genetic_phenotypic_correlation_plot_with_subtypes_and_defs_ALL_RG.csv")




#FOR UEG POSTER: see ueg_poster_rg_plot.R

data_combos = c("Discovery cohort","Bellygenes data only","UK Biobank data only")
Rmini = R[plot==TRUE & trait1 %in% data_combos]
Rmini[,trait1:=factor(trait1,levels=data_combos)]
saveRDS(Rmini,"tables/Rmini_for_rgs_big3_ueg_plot.rds")
Rmini = readRDS("tables/Rmini_for_rgs_big3_ueg_plot.rds")

#see ueg_poster_rg_plot.R

#FOR NatGen comment: see natgen_rg_plot.R
saveRDS(R,"tables/R_for_natgen_rg_plot.rds")
