library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(meta)

read_sumstats = function(path){
  summ = fread(path)
  
  #CONSIDER TEST_ALLELE AS ALLELE1 for MAURO
  #See mauro_ldhub_format, because ldsc doens't understand ALLEL0/ALLEL1
  #This is different from freq_table, where it came from.
  bolt_cols = c("SNP","ALLELE1","ALLELE0","A1FREQ","P_BOLT_LMM","BETA","SE")
  mauro_cols = c("SNP","Test_Allele","Other_Allele","Test_Allele_Freq","P-value","BETA","SE")
  #these are the annotated .pos ones, SNP is rsid from bim. Should be markername if you're doing by by SNP, 'cause then there's no annotation step where the SNP column is added.
  metal_cols = c("SNP","Allele1","Allele2","Freq1","P-value","Effect","StdErr")
  withN_cols = c("SNP","ALLELE1","ALLELE0","A1FREQ","P_BOLT_LMM","BETA","SE")
  standard_cols = c("SNP","ALLELE1","ALLELE0","A1FREQ","P","BETA","SE")
  
  if(any(grepl("^MU_TRANSFORM$",colnames(summ)))){ #Metal, can also do based on path
    setnames(summ,old=withN_cols,new=standard_cols)
  }
  
  else if(any(grepl("BOLT",colnames(summ)))){ #Bolt, but if they come form sumstats_with_N then these already have transformed betas
    summ[,Zscore:=BETA/SE]
    setnames(summ,old=bolt_cols,new=standard_cols)
  }
  
  else if(any(grepl("HetPVal",colnames(summ)))){ #Mauro, can also do based on path
    setnames(summ,old=mauro_cols,new=standard_cols)
  }
  
  else if(any(grepl("^Freq1$",colnames(summ)))){ #Metal, can also do based on path
    setnames(summ,old=metal_cols,new=standard_cols)
  }
  
  else if(any(grepl("im.num.0",colnames(summ)))){
    setnames(summ,
             old=c("assay.name","pvalue","effect","stderr","position","freq.b"),
             new=c("SNP","P","BETA","SE","BP","A1FREQ"))
    summ[,CHR:=as.numeric(gsub("chr","",scaffold))]
    summ[,c("ALLELE0","ALLELE1"):=tstrsplit(alleles,"/")]
  }
  
  #Determine how many of these columns match
  
  summ[,P:=as.numeric(P)]
  summ[,ALLELE0:=toupper(ALLELE0)][,ALLELE1:=toupper(ALLELE1)]
  if(!any(grepl("allele_lex",colnames(summ)))){
    summ[,allele_lex_min:=pmin(ALLELE0,ALLELE1)]
    summ[,allele_lex_max:=pmax(ALLELE0,ALLELE1)]
  }
  
  summ[,CHR:=as.numeric(CHR)]
  summ = summ[CHR%in%(1:22)] #CAREFUL
  
  pattern = paste(paste0("(^",c("CHR","BP","allele_lex_min","allele_lex_max",standard_cols),"$)"),collapse="|")
  subset_cols = grep(pattern,colnames(summ),value=TRUE)
  
  return(summ[,subset_cols,with=FALSE])
}



setwd("~/work/ibs/")


diagnosis_types = c("final_cc.sr",
                    "final_cc.Q.prev.diag.yes",
                    "final_cc.diag",
                    "final_cc.Q")

summ_files = c("meta/runs/metal_ICD_ROME_EURUSA_Qonly_Qnon/METAANALYSIS1.TBL.pos",
               "meta/runs/metal_Qonly_Qnon_any/METAANALYSIS1.TBL.pos",
               paste0("meta/sumstats_with_N/",diagnosis_types,".sumstats"),
               "resource/ibs_replication_summary_stats.csv",
               "meta/runs/metal_MAURO/METAANALYSIS1.TBL.pos",
               "meta/runs/metal_ROME_Q/METAANALYSIS1.TBL.pos")

names(summ_files) <- c("metal_ICD_ROME_EURUSA_Qonly_Qnon",
                       "metal_Qonly_Qnon_any",
                       diagnosis_types,
                       "ibs_replication_summary_stats.csv",
                       "metal_MAURO",
                       "metal_ROME_Q")

GWAS_legend = c("ibs_replication_summary_stats.csv"="23andMe replication",
                "metal_MAURO"="Bellygenes cohorts",
                "metal_ROME_Q"="UKB DHQ Rome III + Bellygenes Rome",
                "final_cc.Q"="UKB: DHQ Rome III",
                "final_cc.Q.prev.diag.yes"="UKB: DHQ self-report",
                "final_cc.sr"="UKB: Unprompted self-rep.",
                "final_cc.diag"="UKB: Hospital ICD-10",
                "metal_Qonly_Qnon_any"="UKB cases in the discovery cohort",
                "metal_ICD_ROME_EURUSA_Qonly_Qnon"="Discovery cohort")

#Extracting lead SNPs from replication table
R = fread("/home/cq/work/ibs/clump/replication/natgen_replication_for_all_hits_extra_cols.csv")
R[,replicated_incl_by_proxy:=any(replicated_1),by=c("analysis","lead_SNP")]
R[,replicated_incl_by_proxy_any_analysis:=any(replicated_1),by=c("lead_SNP")]
relevant_analyses = names(summ_files) #c("metal_ICD_ROME_EURUSA_Qonly_Qnon","final_cc.sr","final_cc.Q.prev.diag.yes","final_cc.diag","final_cc.Q")
R = R[analysis%in%relevant_analyses & new_gwsig_expanded_to_clump==TRUE & SNP_status=="lead"]
setorder(R,CHR,BP)

all_hits = R[,unique(SNP)]


low_mem = function(path,all_hits){
  S = read_sumstats(path)
  S = S[SNP%in%all_hits,]
  return(S)
}

#SLOW, UNCOMMENT FOR RE-BUILD####################################
# summ = list()
# for(f in 1:length(summ_files)){
#   summ[[f]] = low_mem(summ_files[f],all_hits)
# }
# 
# Shs = summ
# 
# names(Shs) = names(summ_files)
# for(s in 1:length(Shs)){
#   Shs[[s]][,GWAS:=names(summ_files)[s]]
# }
# 
# saveRDS(Shs,"~/work/ibs/tables/natgen_forest_data_bigger.RDS")
################################################################


Shs = readRDS("~/work/ibs/tables/natgen_forest_data_bigger.RDS")
X = rbindlist(Shs,use.names=TRUE)[,.(SNP,CHR,BP,GWAS,P,BETA,SE,allele_lex_min,allele_lex_max,ALLELE0,ALLELE1)]
setorder(X,CHR,BP)

X[SNP!="rs622871"] #We ignore this because it's in moderate LD with another chr6 hit, and otherwise forest plot makes it seem as if the discovery yields 7 independent hits.


#Align all with discovery
X[,A0.discovery:=ALLELE0[GWAS=="metal_ICD_ROME_EURUSA_Qonly_Qnon"],by="SNP"]
X[,A1.discovery:=ALLELE1[GWAS=="metal_ICD_ROME_EURUSA_Qonly_Qnon"],by="SNP"]
X[ALLELE1==A0.discovery & ALLELE0==A1.discovery,flipped.wrt.discovery:=TRUE]
X[flipped.wrt.discovery==TRUE,c("ALLELE0","ALLELE1"):=rev(.(ALLELE0,ALLELE1))]
X[flipped.wrt.discovery==TRUE, c("BETA"):=-BETA]

print("Failed to align, check manually:")
X[!(ALLELE1==A1.discovery & ALLELE0==A0.discovery)] #Looks good, insertion==CT

#Flip all where discovery is negative
X[,flipped.because.negative:=any(GWAS=="metal_ICD_ROME_EURUSA_Qonly_Qnon" & BETA<0),by=SNP]
X[flipped.because.negative==TRUE, c("ALLELE0","ALLELE1"):=rev(.(ALLELE0,ALLELE1))]
X[flipped.because.negative==TRUE, c("BETA"):=-BETA]



X[,detail:=paste0(SNP," bp ",BP," chr ",CHR)]
X[,detail:=factor(detail,levels=rev(unique(detail)))]
X[,Zscore:=BETA/SE]

X[,BP.bin:=round(BP/250e3)*250e3]
X[,signal.no:=frank(BP.bin,ties.method = "dense"),by=CHR]
#Heterogeneity test:
GWAS_compared = diagnosis_types
for(hit in all_hits){
  het = metagen(X[SNP==hit & GWAS%in%GWAS_compared,BETA], X[SNP==hit & GWAS%in%GWAS_compared,SE])
  X[SNP==hit & GWAS%in%GWAS_compared,pval.Q:=het$pval.Q]
  X[SNP==hit & GWAS%in%GWAS_compared,Q:=het$Q]
}

#Add replication data from above

X = merge(X,
          unique(R[,.(SNP=lead_SNP,replicated=replicated_incl_by_proxy_any_analysis)]),
          by="SNP")



X[is.na(replicated),replicated:=FALSE]

X[,diff(range(BP)),by=CHR]
setorder(X,CHR,BP)
X[,signal:=paste("CHR",CHR)]
X[,multiple:=max(signal.no)>1,by=CHR]

X[,signal:=factor(signal,levels=unique(signal))]

X[,GWAS_name:=factor(GWAS_legend[GWAS],levels=GWAS_legend)]

X[,OR_raw:=exp(BETA)]
X[,CI.lower:=exp(BETA-1.96*SE)]
X[,CI.upper:=exp(BETA+1.96*SE)]

important_GWAS = c("Discovery cohort","23andMe replication")
X[,importance:=as.numeric(ifelse(GWAS_name%in%important_GWAS,2,1))]
X[,label:=paste0(X$SNP," on ",signal)]

X[,label:=factor(label,levels=unique(label))]
X[,pvaltext:=paste0("\np=",gsub("e-0","e-",formatC(P, format = "e", digits = 1, zero.print = FALSE)))]
theme_set(theme_cowplot())

X[,sig_heterogeneity:=NULL]
X[pval.Q<0.05,sig_heterogeneity:="sighet"]

X[P>5e-8,gw_sig:="Not genome-wide significant"]
X[P<=5e-8,gw_sig:="Genome-wide significant"]

X[,plot:=TRUE]


X[SNP%in%ignore,plot:=FALSE]



X[plot==TRUE,plot_at:=ifelse(P[GWAS_name=="Discovery cohort"]<5e-8,"top","bottom"),by="SNP"]

p.f1 = ggplot(X[plot_at=="top"],aes(y=GWAS_name)) + 
  geom_vline(xintercept=1,lty="solid",color="grey60") + 
  geom_errorbarh(aes(xmin=CI.lower,xmax=CI.upper),show.legend=FALSE) +
  geom_point(aes(x=OR_raw,shape=gw_sig),fill="white",size=3) + 
  facet_wrap(~label,ncol=3,scales="free_x") +
  scale_linetype_manual(labels=c("TRUE"="Yes","FALSE"="No"),
                        values=c("TRUE"="solid","FALSE"="dashed")) +
  #background_grid() +
  labs(y="", x=" ",col="GWAS",fill="",shape="") +
  geom_text(aes(x=OR_raw,label=pvaltext),col="black",alpha=0.5,cex=3) +
  geom_rect(data=unique(X[plot_at=="top",sig_heterogeneity,by=label]),
            aes(xmin=-Inf,xmax=Inf,
                ymax=max(match(GWAS_compared,names(GWAS_legend))+0.5),
                ymin=min(match(GWAS_compared,names(GWAS_legend)))-0.5, 
                fill=sig_heterogeneity, y=NULL),
            col=FALSE,alpha=0.12) +
  geom_rect(data=unique(X[plot_at=="top",replicated,by=label]), 
            aes(xmin=-Inf,xmax=Inf,
                ymax=max(match("23andMe replication",GWAS_legend)+0.5),
                ymin=min(match("23andMe replication",GWAS_legend))-0.5,
                fill=replicated, y=NULL),
            col=FALSE,alpha=0.08) + 
  theme(legend.position="top") +
  guides(
    shape = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  ) +
  scale_shape_manual(breaks=c("Genome-wide significant","Not genome-wide significant"),
                     values=c("Genome-wide significant"=16,"Not genome-wide significant"=21)) + 
  scale_fill_manual(breaks = c("TRUE", "FALSE","sighet"),
                    labels=c("TRUE"="Replicated","FALSE"="Not replicated","sighet"="Significant heterogeneity within UKB"),
                    values=c("TRUE"="green","FALSE"="red","sighet"="skyblue"))

p.f2 = ggplot(X[plot_at=="bottom"],aes(y=GWAS_name)) + 
  geom_vline(xintercept=1,lty="solid",color="grey60") + 
  geom_errorbarh(aes(xmin=CI.lower,xmax=CI.upper),show.legend=FALSE) +
  geom_point(aes(x=OR_raw,shape=gw_sig),fill="white",size=3) + 
  facet_wrap(~label,ncol=4,scales="free_x") +
  scale_linetype_manual(labels=c("TRUE"="Yes","FALSE"="No"),
                        values=c("TRUE"="solid","FALSE"="dashed")) +
  #background_grid() +
  labs(y="", x="OR",col="GWAS",fill="",shape="") +
  geom_text(aes(x=OR_raw,label=pvaltext),col="black",alpha=0.5,cex=3) +
  geom_rect(data=unique(X[plot_at=="bottom",sig_heterogeneity,by=label]),
            aes(xmin=-Inf,xmax=Inf,
                ymax=max(match(GWAS_compared,names(GWAS_legend))+0.5),
                ymin=min(match(GWAS_compared,names(GWAS_legend)))-0.5, 
                fill=sig_heterogeneity, y=NULL),
            col=FALSE,alpha=0.12) +
  geom_rect(data=unique(X[plot_at=="bottom",replicated,by=label]), 
            aes(xmin=-Inf,xmax=Inf,
                ymax=max(match("23andMe replication",GWAS_legend)+0.5),
                ymin=min(match("23andMe replication",GWAS_legend))-0.5,
                fill=replicated, y=NULL),
            col=FALSE,alpha=0.08) + 
  theme(legend.position="bottom") +
  guides(
    shape = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  ) +
  scale_shape_manual(breaks=c("Genome-wide significant","Not genome-wide significant"),
                     values=c("Genome-wide significant"=16,"Not genome-wide significant"=21)) + 
  scale_fill_manual(breaks = c("TRUE", "FALSE","sighet"),
                    labels=c("TRUE"="Replicated","FALSE"="Not replicated","sighet"="Significant heterogeneity within UKB"),
                    values=c("TRUE"="green","FALSE"="red","sighet"="skyblue"))


p.text.discovery = ggplot() + 
  annotate("text", x = 0, y = 0, size=8, angle=270, label = "Discovery cohort loci") + 
  theme_void()

p.text.other = ggplot() + 
  annotate("text", x = 0, y = 0, size=8, angle=270, label = "Other loci") + 
  theme_void()


p.forest.all = plot_grid(get_legend(p.f2 + theme(legend.justification="center", legend.text=element_text(size=16))), NULL,
                         p.f1 + theme(legend.position="none"), p.text.discovery,
                         p.f2 + theme(legend.position="none"), p.text.other,
                         ncol=2,
                         rel_heights = c(1,10,10),
                         rel_widths= c(20,1))

p.forest.all = plot_grid(get_legend(p.f2 + theme(legend.justification="center", legend.text=element_text(size=16))),
                         plot_grid(p.f1 + theme(legend.position="none"), p.text.discovery,
                                   p.f2 + theme(legend.position="none"), p.text.other,
                                   ncol=2,
                                   rel_widths=c(30,1)),
                         ncol=1,
                         rel_heights = c(1,15))

#p.forest.all
ggsave("~/work/ibs/plots/natgen_forest_discoveryfocus_text.pdf",p.forest.all,width=800/47,height=565/47)


