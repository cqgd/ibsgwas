library(data.table)
setwd("~/work/ibs")
#library(coloc)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

path_to_gene_tissue = function(X){
  split = tstrsplit(X$path,"/")
  n = length(split)
  split[[n]] = gsub(".RDS$","",split[[n]])
  X[,tissue:=split[[n-1]]]
  X[,gene_id:=split[[n]]]
}

TissueNames = fread("resource/gtex/tissue_names.csv") #converts between underscores and hyphenated/bracketed names
tn = TissueNames$nice_name; names(tn) = TissueNames$systematic_name

coloc_dir = "coloc"
coloc_files = list.files("coloc",pattern="^ENSG.*RDS$",recursive=TRUE)
coloc_paths = paste0(coloc_dir,"/",coloc_files)
L = lapply(coloc_paths,readRDS)
L = lapply(L, function(C){C$summary})
C = data.table(do.call(rbind,L))
C[,path:=coloc_paths]
path_to_gene_tissue(C)
C$tissue = tn[C$tissue]

E = fread("/home/cq/work/ibs/resource/gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct")
E_ = melt(E,id.vars=c("gene_id","Description"),variable.name="tissue",value.name = "expression")


X = merge(C,E_,by=c("gene_id","tissue"),all=TRUE)
rm(E);rm(E_);gc()

regions_dir = "resource/gtex/regions"
regions_files = list.files(regions_dir,pattern=".csv$")
regions_paths = paste0(regions_dir,"/",regions_files)
R = rbindlist(lapply(regions_paths,fread))
R = unique(R)
R$tissue = tn[R$tissue]

X = merge(X,R,by=c("gene_id","tissue"),all.x=TRUE)

H = fread("clump/replication/replication_for_all_hits.csv")
H = H[analysis=="metal_ICD_ROME_EURUSA_Qonly_Qnon",.(clump, SNP,SNP_status,CHR=as.character(CHR),BP_START=BP-0,BP_END=BP+0)] #The zeros are necessary, otherwise setkey shuffles rsids and bp
setkey(H,CHR,BP_START,BP_END)
leads = H[SNP_status=="lead",SNP]

lead_overlap = as.character(lapply(strsplit(X$SNPs,","),function(three_snps,leads){intersect(three_snps,leads)},leads=leads))
X[,lead_overlap:=lead_overlap]

X[,keep:=any(PP.H4.abf>=0.5),by=gene_id]
X = X[keep==TRUE]
X[,coloc_signif:=!is.na(PP.H4.abf) & (round(PP.H4.abf,2)>=0.5)]

X = X[lead_overlap!="character(0)"]

#Normalize data

normalize = function(v){
  v = v-mean(v)
  v = v/sd(v)
  return(v)
}

X[,expression_normalized:=normalize(expression),by=gene_id]

fwrite(X[coloc_signif==TRUE],"coloc/coloc_signif_pooled.csv")

border_cols = c("TRUE"="goldenrod3","FALSE"="transparent","NA"="transparent")

ggplot(X, aes(x=tissue,y=Description)) + 
  geom_tile(aes(fill=expression_normalized,color=coloc_signif),width=0.92,height=0.92,size=1.1) + 
  #geom_text(aes(label=round(100*PP.H4.abf),alpha=PP.H4.abf)) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0), strip.text.y = element_text(angle = 0)) +
  labs(fill="Expression\n(Gene-normalized TPM)",x="Tissue",y="Gene",color="Colocalization\n(GWAS vs expression)") +
  scale_fill_gradient(low="white",high="forestgreen") +
  scale_color_manual(values=border_cols) +
  facet_grid(lead_overlap~.,scales = "free",space="free_y")

ggplot(X, aes(x=tissue,y=Description)) + 
  geom_tile(aes(fill=expression_normalized),width=0.92,height=0.92,size=1.1) + 
  geom_text(aes(label=round(100*PP.H4.abf),alpha=PP.H4.abf)) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0), strip.text.y = element_text(angle = 0)) +
  labs(fill="Expression\n(Gene-normalized TPM)",x="Tissue",y="Gene",color="Colocalization\n(GWAS vs expression)") +
  scale_fill_gradient(low="white",high="forestgreen") +
  scale_color_manual(values=border_cols) +
  facet_grid(lead_overlap~.,scales = "free",space="free_y")

ggplot(X, aes(x=tissue,y=Description)) + 
  geom_tile(aes(fill=expression_normalized),width=0.92,height=0.92,size=1.1) + 
  geom_text(aes(label=round(100*PP.H3.abf),alpha=PP.H3.abf)) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0), strip.text.y = element_text(angle = 0)) +
  labs(fill="Expression\n(Gene-normalized TPM)",x="Tissue",y="Gene",color="Colocalization\n(GWAS vs expression)") +
  scale_fill_gradient(low="white",high="forestgreen") +
  scale_color_manual(values=border_cols) +
  facet_grid(lead_overlap~.,scales = "free",space="free_y")


#Coloc expression plot, leaving out HLA region SNP for paper:

expr_range = X[lead_overlap!="rs2736155" & !grepl("^RP11",Description),range(expression_normalized)]
expr_max_non_testis = X[lead_overlap!="rs2736155" & !grepl("^RP11",Description) & tissue!="Testis",max(expression_normalized)]

p.coloc_expr = ggplot(X[lead_overlap!="rs2736155" & !grepl("^RP11",Description)], aes(x=tissue,y=Description)) + 
  geom_tile(aes(fill=expression_normalized,color=coloc_signif),width=0.92,height=0.92,size=1.1) + 
  #geom_text(aes(label=round(100*PP.H4.abf),alpha=PP.H4.abf)) +
  theme(axis.text.x = element_text(angle = 270+45, hjust = 0), strip.text.y = element_text(angle = 0)) +
  labs(fill="Expression\n(Gene-normalized TPM)",x="Tissue",y="Gene",color="Colocalization\n(GWAS vs expression)") +
  #scale_fill_gradient(low="white",high="forestgreen", limits=c(NA,expr_max_non_testis),oob=scales::squish) +
  scale_color_manual(values=border_cols) +
  facet_grid(lead_overlap~.,scales = "free",space="free_y") +
  scale_fill_gradientn (
    colours = colorRampPalette(c("white", "forestgreen"))(20),
    values = c(seq(0, expr_max_non_testis/max(expr_range), length.out = 19), 1))

#for non-linear color scale, see also https://stackoverflow.com/questions/12834802/non-linear-color-distribution-over-the-range-of-values-in-a-geom-raster

p.coloc_expr

  #scale_fill_gradientn (
   # colours = colorRampPalette(c("white", "forestgreen"))(20),
    #values = c(0, seq(qn01[1], qn01[2], length.out = 18), 1))


ggsave("~/work/ibs/plots/coloc_expr.pdf",p.coloc_expr,width=15,height=7)
ggsave("~/work/ibs/plots/coloc_expr.svg",p.coloc_expr,width=15,height=7)


C_ = melt(C,measure.vars=grep("PP",colnames(C),value=TRUE),variable.name = "H")

C_[,H:=gsub("(^PP.)|(.abf$)","",H)]
C_[,path:=gsub("^coloc/","",path)]
ggplot(C_, aes(x=H,y=path)) + 
  geom_tile(aes(fill=100*value)) + 
  geom_text(aes(label=round(value*100))) + 
  scale_fill_gradient(low="white",high="forestgreen") + 
  labs(fill="PP*100",x="Hypothesis",y="Region within 1 Mbp of gene\nMust be around eGene with qval<0.05, and overlap with GWAS hits") + 
  ggtitle("Colocalization between signals for IBS risk and gene expression\ndataset 1: Pooled IBS meta-analysis\ndataset 2: GTEx eQTLs")

library(cowplot)
ggsave("plots/coloc.pdf",width = 15,height=20)


#X[Description=="LY6G6D",.(expression,PP.H4.abf,tissue)]

