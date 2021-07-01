library(data.table)
setwd("~/work/ibs")

grep_start = function(col,code){
  grepl(paste0("^",code),col)
}

save_zip = function(object,path){
  fwrite(object,path)
  zip(paste0(path,".zip"),path)
}

add_any = function(pattern,M,suffix="_any"){
  cols = grep(pattern, colnames(M))
  new_col = paste0(pattern,suffix)
  M[,(new_col) := rowSums(.SD,na.rm=TRUE)>=1, .SDcols=cols] 
}

condition_count = function(data, FI, initial=TRUE){
  n = nrow(FI)
  skip_totals = numeric(n+1)
  #on the last run, no filters get skipped
  
  for(s in 1:(n+1)){
    pass_alone = numeric(n)
    pass_cumulative = numeric(n)
    cumulative = rep(initial,nrow(data))
    
    for(i in 1:n){
      if(i==s){
        alone = cumulative
      }else{
        alone = data[[ FI[i,col] ]]
        if(FI[i,not]=="not"){
          alone = !alone
        }
      }
      
      cumulative = get(FI[i,logic])(cumulative,alone)
      
      if(s==(n+1)){
        pass_alone[i] = sum(alone, na.rm=TRUE)
        pass_cumulative[i] = sum(cumulative, na.rm=TRUE)
      }
    }
    skip_totals[s] = sum(cumulative, na.rm=TRUE)
    cat(s,"/",n+1,"\n")
  }
  skip_diff = skip_totals[1:n]  - skip_totals[n+1]

  out = data.table(FI, pass_alone, pass_cumulative, skip_diff)
  if(initial==TRUE){
    first=list("","","all","beginning with all",nrow(data),nrow(data),"")
  }else{
    first=list("","","all","beginning without any",nrow(data),0,"")
  }
  out = rbind(first,out,use.names=TRUE)

  return(list(pass=cumulative, overview=out))
}

filter = function(cols,logic,not,start,desc=cols,ignore="_any$"){
  cols = grep(paste0("^",start),fs,value=TRUE)
  ig = grep(ignore,cols)
  if(length(ig)){
    cols = cols[-ig]
  }
  
  n=length(cols)
  data.table(logic=rep(logic,n),
             not=rep(not,n),
             col=cols,
             desc=desc)
}


eid_col = "eid"
diag_sr = fread("filters/diagnoses_selfreported.csv", key=eid_col)
quest = fread("filters/questionnaire.csv", key=eid_col)
surg = fread("filters/surgeries.csv", key=eid_col)
sqc = fread("filters/sample_qc.csv", key=eid_col)

M = diag_sr[quest][surg][sqc]


#New: adding GAD7 scores used for anxiety cases
phq9_gad7 = fread("pheno/ukb39562_GAD7_PHQ9.csv",key=eid_col)
M = merge(M,phq9_gad7,all.x=TRUE,by="eid",sort=FALSE)

summarize_cols = c("diag.exclusion","diag.surrogate",
                   "sr.exclusion","sr.exclusion",
                   "surg.exclusion")

sapply(summarize_cols,add_any,M=M)

fs = colnames(M)


case_base = rbind(filter(fs,"&","not","Q_q.cel.glut.gold","Questionnaire celiac/gluten sensitivity (F8:3,4,5)"),
                  filter(fs,"&","not","diag.exclusion","ICD10 exclusion incl. children"),
                  filter(fs,"&","not","sr.exclusion","Self-reported"),
                  filter(fs,"&","not","surg.exclusion","OPCS4 exclusion"),
                  #filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)"))
                  filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)|(used.in.pca)|(excess.relatives)"))

FI.case.diag = rbind(filter(fs,"|","","diag.IBS", "ICD10 IBS"),
                case_base)

FI.case.sr = rbind(filter(fs,"|","","sr.IBS", "SR IBS"),
                case_base)

FI.case.Q = rbind(filter(fs,"|","","Q_cases$","Questionnaire algorithm IBS"),
                case_base)

FI.case.Q_prev.diag.yes = rbind(filter(fs,"|","","Q_prev.diag.yes","Questionnaire prev. diag. IBS"),
                case_base)

FI.case.any = rbind(filter(fs,"|","","diag.IBS", "ICD10 IBS"),
                    filter(fs,"|","","sr.IBS", "SR IBS"),
                    filter(fs,"|","","Q_cases$","Questionnaire algorithm IBS"),
                    filter(fs,"|","","Q_prev.diag.yes","Questionnaire prev. diag. IBS"),
                    case_base)

FI.cont = rbind(filter(fs,"&","not","diag.IBS", "ICD10 IBS"),
                filter(fs,"&","not","sr.IBS", "Self-reported IBS"),
                filter(fs,"&","not","diag.exclusion","ICD10 exclusion incl. children"),
                filter(fs,"&","not","diag.surrogate","ICD10 surrogate incl. children"),
                filter(fs,"&","not","sr.exclusion","Self-reported exclusion"),
                filter(fs,"&","not","sr.surrogate","Self-reported surrogate"),
                filter(fs,"&","not","surg.exclusion","OPCS4 exclusion"),
                filter(fs,"&","not","surg.gall","OPCS4 gall"), #can be made exclusion, then ignored for cases.
                filter(fs,"&","","Q_prev.diag.no","Questionnaire no prev. diag. IBS"),
                filter(fs,"&","not","Q_q.cel.glut.gold","Questionnaire celiac/gluten sensitivity (F8:3,4,5)"),
                filter(fs,"&","","Q_conts$","Questionnaire control (q1)"),
                filter(fs,"&","","Q_conts.strict$","Questionnaire control (q1,9,10)"),
                #filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)"))
                filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)|(used.in.pca)|(excess.relatives)"))

  


FI.cont.relaxed = rbind(filter(fs,"&","not","diag.IBS", "ICD10 IBS"),
                filter(fs,"&","not","sr.IBS", "Self-reported IBS"),
                filter(fs,"&","not","diag.exclusion","ICD10 exclusion incl. children"),
                filter(fs,"&","not","diag.surrogate","ICD10 surrogate incl. children"),
                filter(fs,"&","not","sr.exclusion","Self-reported exclusion"),
                filter(fs,"&","not","sr.surrogate","Self-reported surrogate"),
                filter(fs,"&","not","surg.exclusion","OPCS4 exclusion"),
                filter(fs,"&","not","surg.gall","OPCS4 gall"), #can be made exclusion, then ignored for cases.
                #filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)"))
                filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)|(used.in.pca)|(excess.relatives)"))


FI.quant = rbind(filter(fs,"&","not","diag.exclusion","ICD10 exclusion incl. children"),
                 filter(fs,"&","not","diag.surrogate","ICD10 surrogate incl. children"),
                 filter(fs,"&","not","sr.exclusion","Self-reported exclusion"),
                 filter(fs,"&","not","sr.surrogate","Self-reported surrogate"),
                 filter(fs,"&","not","surg.exclusion","OPCS4 exclusion"),
                 filter(fs,"&","not","surg.gall","OPCS4 gall"), #can be made exclusion, then ignored for cases.
                 filter(fs,"&","not","Q_q.cel.glut.gold","Questionnaire celiac/gluten sensitivity (F8:3,4,5)"),
                 #filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)"))
                 filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)|(used.in.pca)|(excess.relatives)"))

FI.quant.mauro = rbind(filter(fs,"&","not","Q_q.cel.glut.gold","Questionnaire celiac/gluten sensitivity (F8:3,4,5)"),
                 filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)|(used.in.pca)|(excess.relatives)"))


FI.case.respondent = rbind(filter(fs,"&","","Q.participant$","Digestive health questionnaire respondent"),
                           filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)|(used.in.pca)|(excess.relatives)"))


FI.cont.respondent = rbind(filter(fs,"&","not","Q.participant$","Digestive health questionnaire respondent"),
                           filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)|(used.in.pca)|(excess.relatives)"))



#QC for clinical table:
#Post-infective
FI.post.infective = rbind(filter(fs,"|","","Q_post.infective","Questionnaire post-infective IBS"),
                          case_base)

FI.post.infective.C = rbind(filter(fs,"|","","Q_post.infective.C","Questionnaire post-infective IBS, subtype C"),
                            case_base)

FI.post.infective.U = rbind(filter(fs,"|","","Q_post.infective.U","Questionnaire post-infective IBS, subtype U"),
                            case_base)

FI.post.infective.M = rbind(filter(fs,"|","","Q_post.infective.M","Questionnaire post-infective IBS, subtype M"),
                            case_base)

FI.post.infective.D = rbind(filter(fs,"|","","Q_post.infective.D","Questionnaire post-infective IBS, subtype D"),
                            case_base)
#Functional
FI.func.C = rbind(filter(fs,"|","","Q_func.C","Questionnaire functional IBS, subtype C"),
                  case_base)
FI.func.D = rbind(filter(fs,"|","","Q_func.D","Questionnaire functional IBS, subtype D"),
                  case_base)

FI.func.relaxed.C = rbind(filter(fs,"|","","Q_func.relaxed.C","Questionnaire relaxed functional IBS, subtype C"),
                          case_base)
FI.func.relaxed.D = rbind(filter(fs,"|","","Q_func.relaxed.D","Questionnaire relaxed functional IBS, subtype D"),
                          case_base)
#End of clinical table QC

#Anxiety

#Anxiety GWAS
any_anxiety = rbind(filter(fs,"|","","diag.interest_F40","Phobic anxiety (ICD-10)"),
                    filter(fs,"|","","diag.interest_F41","Other anxiety (ICD-10)"),
                    filter(fs,"|","","sr.interest_1287","Self-reported anxiety/panic attacks"),
                    filter(fs,"|","","Q_treatment.anx","Treatment sought or offered for anxiety (DHQ)"))

no_anxiety = rbind(filter(fs,"&","not","diag.interest_F40","No phobic anxiety (ICD-10)"),
                   filter(fs,"&","not","diag.interest_F41","No other anxiety (ICD-10)"),
                   filter(fs,"&","not","sr.interest_1287","No self-reported anxiety/panic attacks"),
                   filter(fs,"&","not","Q_treatment.anx","No treatment sought or offered for anxiety (DHQ)")) #as.logical(1)==TRUE obviously, as.logical(-818)==TRUE also, so leaves only true 0 responses.


any_anxietybroad = rbind(any_anxiety,
                    filter(fs,"|","","GAD7.anx","GAD-7 based anxiety, score>=10"))

no_anxietybroad = rbind(no_anxiety,
                   filter(fs,"&","not","GAD7.anx","GAD-7 based anxiety, score>=10")) #as.logical(1)==TRUE obviously, as.logical(-818)==TRUE also, so leaves only true 0 responses.

FI.case.anx = rbind(any_anxiety,
                   filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)|(used.in.pca)|(excess.relatives)"))

FI.cont.anx = rbind(no_anxiety,
                   filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)|(used.in.pca)|(excess.relatives)"))

FI.case.anxbroad = rbind(any_anxietybroad,
                    filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)|(used.in.pca)|(excess.relatives)"))

FI.cont.anxbroad = rbind(no_anxietybroad,
                    filter(fs,"&","","sqc.","UKBB Genetic QC",ignore="(cumulative)|(used.in.pca)|(excess.relatives)"))


#Anxiety GWAS without IBS

no_ibs = rbind(filter(fs,"&","not","diag.IBS", "ICD10 IBS"),
               filter(fs,"&","not","sr.IBS", "SR IBS"),
               filter(fs,"&","not","Q_cases$","Questionnaire algorithm IBS"),
               filter(fs,"&","not","Q_prev.diag.yes","Questionnaire prev. diag. IBS"))


FI.case.anx.no.ibs = rbind(FI.case.anx,
                          no_ibs)

FI.cont.anx.no.ibs = rbind(FI.cont.anx,
                          no_ibs)

FI.case.anxbroad.no.ibs = rbind(FI.case.anxbroad,
                           no_ibs,
                           filter(fs,"&","","Q_participant","Must have completed the DHQ (to ensure absence of IBS via Rome III)"), #Could also do Q_all_SSS
                           filter(fs,"&","","all_GAD7","For symmetry with IBS sans anxiety analysis"))

FI.cont.anxbroad.no.ibs = rbind(FI.cont.anxbroad,
                           no_ibs,
                           filter(fs,"&","","Q_participant","Must have completed the DHQ (to ensure absence of IBS via Rome III)"),
                           filter(fs,"&","","all_GAD7","For symmetry with IBS sans anxiety analysis"))

#IBS GWAS without anxiety
FI.case.no.anx = rbind(FI.case.any,
                       no_anxiety)

FI.cont.no.anx = rbind(FI.cont,
                       no_anxiety)

FI.cont.no.anx.relaxed = rbind(FI.cont.relaxed, #the only definition restricted to DHQ, so if we do this in non-respondents separately we'll have to relax to non-DHQ standards
                       no_anxiety)


FI.case.no.anxbroad = rbind(FI.case.any,
                       no_anxietybroad,
                       filter(fs,"&","","all_GAD7","Must have completed GAD-7 (to ensure absence of anxiety via GAD-7)"),
                       filter(fs,"&","","Q_participant","Must have completed the DHQ (to ensure absence of anxiety via treatment sought)"))

FI.cont.no.anxbroad = rbind(FI.cont,
                       no_anxietybroad,
                       filter(fs,"&","","all_GAD7","Must have completed GAD-7 (to match cases)"),
                       filter(fs,"&","","Q_participant","Must have completed the DHQ (to match cases)"))

#FI.cont.no.anxbroad.relaxed = rbind(FI.cont.relaxed, #the only definition restricted to DHQ, so if we do this in non-respondents separately we'll have to relax to non-DHQ standards
#                               no_anxietybroad)
#
  


F.case.sr = condition_count(M,FI.case.sr,initial=FALSE) #swapped, fixed feb 11 
F.case.diag = condition_count(M,FI.case.diag,initial=FALSE) #swapped
F.case.Q = condition_count(M,FI.case.Q,initial=FALSE)
F.case.Q.prev.diag.yes = condition_count(M,FI.case.Q_prev.diag.yes,initial=FALSE)
F.case.any = condition_count(M,FI.case.any,initial=FALSE)
F.cont = condition_count(M,FI.cont,initial=TRUE)
F.cont.relaxed = condition_count(M,FI.cont.relaxed,initial=TRUE)
F.quant = condition_count(M,FI.quant,initial=TRUE)
F.quant.mauro = condition_count(M,FI.quant.mauro,initial=TRUE)
F.case.respondent = condition_count(M,FI.case.respondent,initial=TRUE)
F.cont.respondent = condition_count(M,FI.cont.respondent,initial=TRUE)

F.post.infective = condition_count(M,FI.post.infective,initial=FALSE)
F.post.infective.C = condition_count(M,FI.post.infective.C,initial=FALSE)
F.post.infective.U = condition_count(M,FI.post.infective.U,initial=FALSE)
F.post.infective.M = condition_count(M,FI.post.infective.M,initial=FALSE)
F.post.infective.D = condition_count(M,FI.post.infective.D,initial=FALSE)
F.func.C = condition_count(M,FI.func.C,initial=FALSE) 
F.func.D = condition_count(M,FI.func.D,initial=FALSE) 
F.func.relaxed.C = condition_count(M,FI.func.relaxed.C,initial=FALSE) 
F.func.relaxed.D = condition_count(M,FI.func.relaxed.D,initial=FALSE) 

F.case.anx = condition_count(M,FI.case.anx,initial=FALSE)
F.cont.anx = condition_count(M,FI.cont.anx,initial=TRUE)
F.case.anx.no.ibs = condition_count(M,FI.case.anx.no.ibs,initial=FALSE)
F.cont.anx.no.ibs = condition_count(M,FI.cont.anx.no.ibs,initial=TRUE)
F.case.no.anx = condition_count(M,FI.case.no.anx, initial=FALSE)
F.cont.no.anx = condition_count(M,FI.cont.no.anx, initial=TRUE)
F.cont.no.anx.relaxed = condition_count(M,FI.cont.no.anx.relaxed, initial=TRUE)



F.case.anxbroad = condition_count(M,FI.case.anxbroad,initial=FALSE)
F.cont.anxbroad = condition_count(M,FI.cont.anxbroad,initial=TRUE)
F.case.anxbroad.no.ibs = condition_count(M,FI.case.anxbroad.no.ibs,initial=FALSE)
F.cont.anxbroad.no.ibs = condition_count(M,FI.cont.anxbroad.no.ibs,initial=TRUE)
F.case.no.anxbroad = condition_count(M,FI.case.no.anxbroad, initial=FALSE)
F.cont.no.anxbroad = condition_count(M,FI.cont.no.anxbroad, initial=TRUE)
#F.cont.no.anxbroad.relaxed = condition_count(M,FI.cont.no.anxbroad.relaxed, initial=TRUE)

M[,final_case.diag:=F.case.diag$pass]
M[,final_case.sr:=F.case.sr$pass]
M[,final_case.Q:=F.case.Q$pass]
M[,final_case.Q.prev.diag.yes:=F.case.Q.prev.diag.yes$pass]
M[,final_case.any:=F.case.any$pass]
M[,final_cont:=F.cont$pass]
M[,final_cont.relaxed:=F.cont.relaxed$pass]
M[,final_quant:=F.quant$pass]
M[,final_quant.mauro:=F.quant.mauro$pass]
M[,final_case.respondent:=F.case.respondent$pass]
M[,final_cont.respondent:=F.cont.respondent$pass]

M[,final_case.post.infective:=F.post.infective$pass] #just added "_case" here for GWAS
M[,final_case.post.infective.C:=F.post.infective.C$pass]
M[,final_case.post.infective.U:=F.post.infective.U$pass]
M[,final_case.post.infective.M:=F.post.infective.M$pass]
M[,final_case.post.infective.D:=F.post.infective.D$pass]
M[,final_case.func.C:=F.func.C$pass]
M[,final_case.func.D:=F.func.D$pass]
M[,final_func.relaxed.C:=F.func.relaxed.C$pass]
M[,final_func.relaxed.D:=F.func.relaxed.D$pass]

M[,final_case.anx:=F.case.anx$pass]
M[,final_cont.anx:=F.cont.anx$pass]
M[,final_case.anx.no.ibs:=F.case.anx.no.ibs$pass]
M[,final_cont.anx.no.ibs:=F.cont.anx.no.ibs$pass]
M[,final_case.no.anx:=F.case.no.anx$pass]
M[,final_cont.no.anx:=F.cont.no.anx$pass]
M[,final_cont.no.anx.relaxed:=F.cont.no.anx.relaxed$pass]

M[,final_case.anxbroad:=F.case.anxbroad$pass]
M[,final_cont.anxbroad:=F.cont.anxbroad$pass]
M[,final_case.anxbroad.no.ibs:=F.case.anxbroad.no.ibs$pass]
M[,final_cont.anxbroad.no.ibs:=F.cont.anxbroad.no.ibs$pass]
M[,final_case.no.anxbroad:=F.case.no.anxbroad$pass]
#M[,final_cont.no.anxbroad.relaxed:=F.cont.no.anxbroad.relaxed$pass]


M[,lapply(.SD,sum,na.rm=TRUE),.SDcols=grep("^final_",colnames(M),value=TRUE)]

fwrite(F.case.diag$overview,"groups/F.case.diag.csv")
fwrite(F.case.sr$overview,"groups/F.case.sr.csv")
fwrite(F.case.Q$overview,"groups/F.case.Q.csv")
fwrite(F.case.Q.prev.diag.yes$overview,"groups/F.case.Q.prev.diag.yes.csv")
fwrite(F.case.any$overview,"groups/F.case.any.csv")
fwrite(F.cont$overview,"groups/F.cont.csv")
fwrite(F.cont.relaxed$overview,"groups/F.cont.relaxed.csv")
fwrite(F.quant$overview,"groups/F.quant.csv")
fwrite(F.quant.mauro$overview,"groups/F.quant.mauro.csv")
fwrite(F.case.respondent$overview,"groups/F.case.respondent.csv")
fwrite(F.cont.respondent$overview,"groups/F.cont.respondent.csv")
fwrite(F.post.infective$overview,"groups/F.post.infective.csv")
fwrite(F.post.infective.C$overview,"groups/F.post.infective.C.csv")
fwrite(F.post.infective.U$overview,"groups/F.post.infective.U.csv")
fwrite(F.post.infective.M$overview,"groups/F.post.infective.M.csv")
fwrite(F.post.infective.D$overview,"groups/F.post.infective.D.csv")
fwrite(F.func.C$overview,"groups/F.func.C.csv")
fwrite(F.func.D$overview,"groups/F.func.D.csv")
fwrite(F.func.relaxed.C$overview,"groups/F.func.relaxed.C.csv")
fwrite(F.func.relaxed.D$overview,"groups/F.func.relaxed.D.csv")
fwrite(F.case.anx$overview,"groups/F.case.anx.csv")
fwrite(F.cont.anx$overview,"groups/F.cont.anx.csv")
fwrite(F.case.anx.no.ibs$overview,"groups/F.case.anx.no.ibs.csv")
fwrite(F.cont.anx.no.ibs$overview,"groups/F.cont.anx.no.ibs.csv")
fwrite(F.case.no.anx$overview,"groups/F.case.no.anx.csv")
fwrite(F.cont.no.anx$overview,"groups/F.cont.no.anx.csv")
fwrite(F.cont.no.anx.relaxed$overview,"groups/F.cont.no.anx.relaxed.csv")

fwrite(F.case.anxbroad$overview,"groups/F.case.anxbroad.csv")
fwrite(F.cont.anxbroad$overview,"groups/F.cont.anxbroad.csv")
fwrite(F.case.anxbroad.no.ibs$overview,"groups/F.case.anxbroad.no.ibs.csv")
fwrite(F.cont.anxbroad.no.ibs$overview,"groups/F.cont.anxbroad.no.ibs.csv")
fwrite(F.case.no.anxbroad$overview,"groups/F.case.no.anxbroad.csv")
fwrite(F.cont.no.anxbroad$overview,"groups/F.cont.no.anxbroad.csv")
#fwrite(F.cont.no.anxbroad.relaxed$overview,"groups/F.cont.no.anxbroad.relaxed.csv")

save_zip(M,"filters/all_integrated.csv")

#Quickly add disease description from ICD10/SR/OPCS4 codes

X_diag_sr_path = "resource/IBS_filter_codes - diagnoses_selfreported.csv"
X_diag_sr = fread(X_diag_sr_path) 
X_surg_path = "resource/IBS_filter_codes - surgeries.csv"
X_surg = fread(X_surg_path)

G.case.diag = fread("groups/F.case.diag.csv")
G.case.sr = fread("groups/F.case.sr.csv")
G.case.Q = fread("groups/F.case.Q.csv")
G.case.Q.prev.diag.yes = fread("groups/F.case.Q.prev.diag.yes.csv")
G.case.any = fread("groups/F.case.any.csv")
G.cont = fread("groups/F.cont.csv")
G.cont.relaxed = fread("groups/F.cont.relaxed.csv")
G.quant = fread("groups/F.quant.csv")
G.quant.mauro = fread("groups/F.quant.mauro.csv")
G.case.respondent = fread("groups/F.case.respondent.csv")
G.cont.respondent = fread("groups/F.cont.respondent.csv")


G.post.infective = fread("groups/F.post.infective.csv")
G.post.infective.C = fread("groups/F.post.infective.C.csv")
G.post.infective.U = fread("groups/F.post.infective.U.csv")
G.post.infective.M = fread("groups/F.post.infective.M.csv")
G.post.infective.D = fread("groups/F.post.infective.D.csv")

G.func.C = fread("groups/F.func.C.csv")
G.func.D = fread("groups/F.func.D.csv")
G.func.relaxed.C = fread("groups/F.func.relaxed.C.csv")
G.func.relaxed.D = fread("groups/F.func.relaxed.D.csv")

G.case.anx=fread("groups/F.case.anx.csv")
G.cont.anx=fread("groups/F.cont.anx.csv")
G.case.anx.no.ibs=fread("groups/F.case.anx.no.ibs.csv")
G.cont.anx.no.ibs=fread("groups/F.cont.anx.no.ibs.csv")
G.case.no.anx=fread("groups/F.case.no.anx.csv")
G.cont.no.anx=fread("groups/F.cont.no.anx.csv")
G.cont.no.anx.relaxed=fread("groups/F.cont.no.anx.relaxed.csv")



annotate = function(G, X, code_col, name_col){
  G[,code := tstrsplit(G$col,split="_")[[2]] ]
  X[,(code_col):= as.character(get(code_col))]
  G[,(name_col):= X[match(G$code, X[[code_col]]), get(name_col)]]
}


triple_annotate = function(G,X_diag_sr,X_surg,path){
  annotate(G, X_diag_sr,code_col="ICD10.code",name_col="ICD10.name")
  annotate(G, X_diag_sr,code_col="SR.code",name_col="SR.name")
  annotate(G, X_surg,code_col="OPCS4.code",name_col="OPCS4.name")
  fwrite(G,path)
  return(G)
}


#Make triple annotation filter a separate func, the below was feasible for one case and one control grouping only.

#G=G.case;X=X_diag_sr;code_col="ICD10.code";name_col="ICD10.name"

triple_annotate(G.post.infective, X_diag_sr, X_surg,"groups/G.post.infective.csv")
triple_annotate(G.post.infective.C, X_diag_sr, X_surg,"groups/G.post.infective.C.csv")
triple_annotate(G.post.infective.U, X_diag_sr, X_surg,"groups/G.post.infective.U.csv")
triple_annotate(G.post.infective.M, X_diag_sr, X_surg,"groups/G.post.infective.M.csv")
triple_annotate(G.post.infective.D, X_diag_sr, X_surg,"groups/G.post.infective.D.csv")


triple_annotate(G.func.C, X_diag_sr, X_surg,"groups/G.func.C.csv")
triple_annotate(G.func.D, X_diag_sr, X_surg,"groups/G.func.D.csv")
triple_annotate(G.func.relaxed.C, X_diag_sr, X_surg,"groups/G.func.relaxed.C.csv")
triple_annotate(G.func.relaxed.D, X_diag_sr, X_surg,"groups/G.func.relaxed.D.csv")

triple_annotate(G.case.anx, X_diag_sr, X_surg,"groups/G.case.anx.csv")
triple_annotate(G.cont.anx, X_diag_sr, X_surg,"groups/G.cont.anx.csv")
triple_annotate(G.case.anx.no.ibs, X_diag_sr, X_surg,"groups/G.case.anx.no.ibs.csv")
triple_annotate(G.cont.anx.no.ibs, X_diag_sr, X_surg,"groups/G.cont.anx.no.ibs.csv")
triple_annotate(G.case.no.anx, X_diag_sr, X_surg,"groups/G.case.no.anx.csv")
triple_annotate(G.cont.no.anx, X_diag_sr, X_surg,"groups/G.cont.no.anx.csv")
triple_annotate(G.cont.no.anx.relaxed, X_diag_sr, X_surg,"groups/G.cont.no.anx.relaxed.csv")
#the below triplets and fwrite can be replaced by triple annotate calls also

annotate(G.case.sr, X_diag_sr,code_col="ICD10.code",name_col="ICD10.name")
annotate(G.case.sr, X_diag_sr,code_col="SR.code",name_col="SR.name")
annotate(G.case.sr, X_surg,code_col="OPCS4.code",name_col="OPCS4.name")
fwrite(G.case.sr,"groups/G.case.sr.csv")

annotate(G.case.diag, X_diag_sr,code_col="ICD10.code",name_col="ICD10.name")
annotate(G.case.diag, X_diag_sr,code_col="SR.code",name_col="SR.name")
annotate(G.case.diag, X_surg,code_col="OPCS4.code",name_col="OPCS4.name")
fwrite(G.case.diag,"groups/G.case.diag.csv")

annotate(G.case.Q, X_diag_sr,code_col="ICD10.code",name_col="ICD10.name")
annotate(G.case.Q, X_diag_sr,code_col="SR.code",name_col="SR.name")
annotate(G.case.Q, X_surg,code_col="OPCS4.code",name_col="OPCS4.name")
fwrite(G.case.Q,"groups/G.case.Q.csv")

annotate(G.case.Q, X_diag_sr,code_col="ICD10.code",name_col="ICD10.name")
annotate(G.case.Q, X_diag_sr,code_col="SR.code",name_col="SR.name")
annotate(G.case.Q, X_surg,code_col="OPCS4.code",name_col="OPCS4.name")
fwrite(G.case.Q,"groups/G.case.Q.csv")

annotate(G.case.Q.prev.diag.yes, X_diag_sr,code_col="ICD10.code",name_col="ICD10.name")
annotate(G.case.Q.prev.diag.yes, X_diag_sr,code_col="SR.code",name_col="SR.name")
annotate(G.case.Q.prev.diag.yes, X_surg,code_col="OPCS4.code",name_col="OPCS4.name")
fwrite(G.case.Q.prev.diag.yes,"groups/G.case.Q.prev.diag.yes.csv")


annotate(G.case.any, X_diag_sr,code_col="ICD10.code",name_col="ICD10.name")
annotate(G.case.any, X_diag_sr,code_col="SR.code",name_col="SR.name")
annotate(G.case.any, X_surg,code_col="OPCS4.code",name_col="OPCS4.name")
fwrite(G.case.any,"groups/G.case.any.csv")

annotate(G.cont, X_diag_sr,code_col="ICD10.code",name_col="ICD10.name")
annotate(G.cont, X_diag_sr,code_col="SR.code",name_col="SR.name")
annotate(G.cont, X_surg,code_col="OPCS4.code",name_col="OPCS4.name")
fwrite(G.cont,"groups/G.cont.csv")

annotate(G.cont.relaxed, X_diag_sr,code_col="ICD10.code",name_col="ICD10.name")
annotate(G.cont.relaxed, X_diag_sr,code_col="SR.code",name_col="SR.name")
annotate(G.cont.relaxed, X_surg,code_col="OPCS4.code",name_col="OPCS4.name")
fwrite(G.cont.relaxed,"groups/G.cont.relaxed.csv")

annotate(G.quant, X_diag_sr,code_col="ICD10.code",name_col="ICD10.name")
annotate(G.quant, X_diag_sr,code_col="SR.code",name_col="SR.name")
annotate(G.quant, X_surg,code_col="OPCS4.code",name_col="OPCS4.name")
fwrite(G.quant,"groups/G.quant.csv")

annotate(G.quant.mauro, X_diag_sr,code_col="ICD10.code",name_col="ICD10.name")
annotate(G.quant.mauro, X_diag_sr,code_col="SR.code",name_col="SR.name")
annotate(G.quant.mauro, X_surg,code_col="OPCS4.code",name_col="OPCS4.name")
fwrite(G.quant.mauro,"groups/G.quant.mauro.csv")

fwrite(G.case.respondent,"groups/G.case.respondent.csv")
fwrite(G.cont.respondent,"groups/G.cont.respondent.csv")

