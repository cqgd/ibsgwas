setwd("~/work/ibs")
library(data.table)
library(ggplot2)
library(coloc)
library(cowplot)
theme_set(theme_cowplot())
library(EnsDb.Hsapiens.v75)
library(ggbio)
ensdb <- EnsDb.Hsapiens.v75

#Let's try knownGene, see https://www.biostars.org/p/275286/
#see also https://www.biostars.org/p/367011/
#TxDb plotting crashes, example doesn't, run and doesn't contain gene names,<3
#Only upside it's a good subset
#Maybe bette rto just subset ensdb based on  https://bioconductor.org/packages/devel/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
#And /home/cq/work/ibs/resource/knownCanonical.txt, which has a subset of transcript IDs

#Load meta-analysis summary stats
meta_path = "meta/runs/metal_ICD_ROME_EURUSA_Qonly_Qnon/METAANALYSIS1.TBL.pos"

M = fread(meta_path)
s = median(M$N_CASE)/median(M$N_TOTAL)
M = M[nchar(SNP)>=1]
M = M[,.(snp=SNP, pvalues=`P-value`, beta=Effect, varbeta=StdErr^2, CHR=as.character(CHR),BP_START=BP-0,BP_END=BP+0)]
setkey(M,CHR,BP_START,BP_END)

#Find regions around our hits
H = fread("clump/replication/replication_for_all_hits.csv")
H = H[analysis=="metal_ICD_ROME_EURUSA_Qonly_Qnon" & SNP_status=="lead",
      .(clump, lead_SNP=SNP,SNP_status,CHR=as.character(CHR),BP_START=BP-100e3,BP_END=BP+100e3)] #The zeros are necessary, otherwise setkey shuffles rsids and bp
setkey(H,CHR,BP_START,BP_END)

O = foverlaps(M,H,nomatch=NULL)
O[,.N,by=lead_SNP]

M_ = split(O,as.factor(O$lead_SNP))

meta_dataset_around_hit = function(SS,s){
  dataset = list(type="cc",
                 s=s,
                 #pvalues=SS$pvalues, 
                 beta=SS$beta, 
                 varbeta=SS$varbeta,
                 snp=SS$snp)
  return(dataset)
}
D = lapply(M_,meta_dataset_around_hit,s=s)

#Perform finemapping with coloc
FM = lapply(D,finemap.abf)
for(i in 1:length(FM)){
  FM[[i]] = data.table(FM[[i]])
  FM[[i]][,lead_SNP:=names(FM[i])]
}
FM = rbindlist(FM)

#Define credible sets
Q = merge(O,FM,by=c("lead_SNP","snp"),all.x=TRUE)
setorder(Q,lead_SNP,-SNP.PP)
Q[,cum_SNP.PP:=cumsum(SNP.PP),by=lead_SNP]
Q[,credible_set:=as.integer(cum_SNP.PP<=0.95)]
Q[lead_SNP==snp,credible_set:=2]

#Add LD information extracted with plink
ld_dir = "ld"
ld_files = list.files(ld_dir,pattern=".ld$")
ld_paths = paste0(ld_dir,"/",ld_files)
L = rbindlist(lapply(ld_paths,fread))
setnames(L,old=c("SNP_A","SNP_B"),new=c("lead_SNP","snp"))
L = L[,.(snp,lead_SNP,R2)]
L[,lead_SNP:=tstrsplit(L$lead_SNP,"::")[1]]
L[,snp:=tstrsplit(snp,"::")[1]]

Q = merge(Q,L,by=c("lead_SNP","snp"),all.x=TRUE)
Q[,credible_set_factor:=factor(credible_set,levels=c(2,1,0),labels=c("lead","in","out"))]


leads = unique(Q$lead_SNP)
plot_paths = paste0("plots/hits/",leads,".pdf")

canonical = fread("~/work/ibs/resource/knownCanonical.txt")[,tstrsplit(V5,"\\.")[1]][[1]]

l=1
l_range = seq_along(leads)
for(l in l_range){
print(leads[l])

Q_ = Q[lead_SNP==leads[l]]



#MANHATTAN PLOT
P.manhattan = ggplot(Q_,aes(x=i.BP_START)) +
  geom_point(aes(y=-log10(pvalues)),data=Q_[credible_set==2],col="goldenrod3",alpha=0.5,cex=5) +
  geom_text(aes(y=-log10(pvalues),label=paste0(lead_SNP,"\n")),data=Q_[credible_set==2]) +
  geom_hline(aes(yintercept=-log10(5e-8)),linetype=2,col="grey") +
  geom_point(aes(y=-log10(pvalues),col=R2)) + #are these the right p-values? Maybe report on causal variant
  theme(legend.position=c(0.9,0.5),legend.direction="vertical") + 
  ylab(expression(-log[10](p))) +
  labs(col=expression(r^2))
#P.manhattan

#FINEMAPPING PLOT
cr_set_cols = c("out"="black","in"="limegreen","lead"="goldenrod3")

P.finemap = ggplot(Q_, aes(x=i.BP_START)) + 
  geom_col(aes(y=SNP.PP,col=credible_set_factor,fill=credible_set_factor)) +
  theme(legend.position=c(0.9,0.5)) +
  scale_color_manual(values=cr_set_cols) +
  scale_fill_manual(values=cr_set_cols) +
  ylab(expression(PP[causal])) +
  labs(fill=c("95% cr. set"),col=c("95% cr. set"))
#P.finemap

#TRANSCRIPT PLOT
anno_start = unique(Q_$BP_START)
anno_end = unique(Q_$BP_END)
anno_end-anno_start
anno_chr = unique(Q_$CH)
gr <- GRanges(seqnames = anno_chr, IRanges(anno_start, anno_end), strand = "*")
if(l==3){ 
  P.annot = autoplot(ensdb, which=AnnotationFilterList(TxIdFilter(canonical), GRangesFilter(gr)), names.expr = "gene_name",ylab="Transcripts",padding=0) #Canonical only, rs2736155 is very gene-dense
}else{
  P.annot = autoplot(ensdb, which=GRangesFilter(gr), names.expr = "gene_name",ylab="Transcripts") #Original, all transcript isoforms, use for "rs5803650"
}
P.annot
#P.annot@ggplot = P.annot@ggplot + coord_cartesian(clip = "off") #disable gene name clipping
N_TRANSCRIPTS = layer_scales(P.annot@ggplot)$y$range$range[2]-0.3 #number of vertical layers, transcript tracks




#CHROMOSOME OVEVRIEW PLOT
P.ideo = Ideogram(genome = "hg19",aspect.ratio = 1/40,zoom.region=c(anno_start,anno_end)) + xlim(GRanges(paste0("chr",anno_chr), IRanges(anno_start, anno_end))) 

#COMBINE
#SUBPLOT_HEIGHTS = c(1,5,3,N_TRANSCRIPTS*(3/7)) #scale height based on how many transcripts are needed
SUBPLOT_HEIGHTS = c(1,5,3,2.5) #fixed hieght
P.combined = tracks(P.ideo,P.manhattan,P.finemap,P.annot + coord_cartesian(clip = "off"),heights=SUBPLOT_HEIGHTS,padding=unit(-1.4,"lines"))


P.combined + coord_cartesian(clip = "off")

#ggbio::ggsave(plot_paths[l], P.combined,width=11,height=sum(SUBPLOT_HEIGHTS)*(9/12))
ggbio::ggsave(plot_paths[l], P.combined,width=11,height=9) #for use with a fixed height
}

