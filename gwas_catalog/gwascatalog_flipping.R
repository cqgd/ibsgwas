library(data.table)
library(googlesheets4)
#Formats FUMAS GWAS Catalog data, matching it to IBS sumstats

#Load GWAS catalog data
X = fread("/home/cq/work/ibs/fuma/results/FUMA_job45052/gwascatalog_dashfix.csv")
X[,gwascat.effectallele:=tstrsplit(Strongest,"-")[[2]]]

#Load IBS sumstats
S = fread("/home/cq/work/ibs/meta/runs/metal_ICD_ROME_EURUSA_Qonly_Qnon/METAANALYSIS1.TBL.pos")

#Flip IBS alleles such that effects are positive
S[,Allele1:=toupper(Allele1)]
S[,Allele2:=toupper(Allele2)]
S[Effect<0,flip:=TRUE]
S[flip==TRUE,Effect:=-Effect]
S[flip==TRUE,c("Allele1","Allele2"):=.(Allele2,Allele1)]
S[flip==TRUE,Freq1:=1-Freq1]

colnames(S) = paste0("IBS.",colnames(S))
setnames(S,old="IBS.SNP",new="snp")

#Merge sumstats and GWAS catalog data
Xm = merge(X,S,by="snp",all.x=TRUE)
Xo = Xm[,.(`IBS Lead SNP`=IndSigSNP, chr, bp, snp, FirstAuth, Link, Study, Trait, `IBS Reference Allele`=IBS.Allele2, `IBS Effect Allele`=IBS.Allele1, `Trait Effect Allele`=gwascat.effectallele, `Trait AF`=RiskAF, `IBS AF`=IBS.Freq1, `IBS Effect`=IBS.Effect, `IBS SE`=IBS.StdErr, `Trait Effect`=OrBeta, `Trait Units and CI`=`95CI`, MappedGene, InitialN, ReplicationN, ReportedGene, Region, Date, Journal)]

#Make the GWAS catalog strand match the IBS sumstats strand
Xo[ ,same_strand:=TRUE]
#...only if the trait effect allele is neither of the IBS alleles (because then flipping clearly isn't sufficient), otherwise we assume the strand is correct
Xo[`Trait Effect Allele`!=`IBS Effect Allele` & `Trait Effect Allele`!=`IBS Reference Allele`, same_strand:=FALSE]
other_strand_from = c("A"="T", "T"="A", "C"="G", "G"="C","?"="?")
Xo[same_strand==FALSE, `Trait Effect Allele`:=other_strand_from[`Trait Effect Allele`]]
Xo[same_strand==FALSE & `Trait Effect Allele`==`IBS Effect Allele` | `Trait Effect Allele`==`IBS Reference Allele`,same_strand:=TRUE]

#Where strands match, add a reference allele for the GWAS catalog trait based on the IBS data
Xo[same_strand==TRUE & `Trait Effect Allele`==`IBS Effect Allele`, `Trait Reference Allele`:=`IBS Reference Allele`]
Xo[same_strand==TRUE & `Trait Effect Allele`==`IBS Reference Allele`, `Trait Reference Allele`:=`IBS Effect Allele`]
Xo[same_strand==FALSE, `Trait Reference Allele`:="?"]

#Where alleles are flipped, change the GWAS catalog data to match the IBS data
Xo[(`Trait Effect Allele`==`IBS Reference Allele` & `Trait Reference Allele` == `IBS Effect Allele`), flipped:=TRUE]
Xo[flipped==TRUE,c("Trait Effect Allele","Trait Reference Allele"):=.(`Trait Reference Allele`, `Trait Effect Allele`)]
Xo[,`Trait AF`:=as.numeric(`Trait AF`)]
Xo[flipped==TRUE,`Trait AF` := 1 -`Trait AF`]
Xo[flipped==TRUE,`Trait Units and CI`:=gsub("decrease","placeholder",`Trait Units and CI`)]
Xo[flipped==TRUE,`Trait Units and CI`:=gsub("increase","decrease",`Trait Units and CI`)]
Xo[flipped==TRUE,`Trait Units and CI`:=gsub("placeholder","increase",`Trait Units and CI`)]

#Compare directions of effect
Xo[,Directions:="?"]
Xo[`Trait Effect Allele`==`IBS Effect Allele` & `IBS Reference Allele` == `Trait Reference Allele`, Directions:="opposite"]
Xo[`Trait Effect Allele`==`IBS Effect Allele` & `IBS Reference Allele` == `Trait Reference Allele` & 
     grepl("increase",`Trait Units and CI`)*2-1==sign(`IBS Effect`), Directions:="identical"]

#Formatting
setorder(Xo,chr,bp)

Xout = Xo[,.(`IBS Lead SNP`,
      Trait, snp, chr, bp,
     `IBS Reference Allele`, `IBS Effect Allele`,
     `Trait Reference Allele`, `Trait Effect Allele`,
     `IBS AF`, `Trait AF`,  
     `IBS Effect`, `IBS SE`, 
     `Trait Effect`, `Trait Units and CI`, `Directions`,
      MappedGene, InitialN, ReplicationN, ReportedGene, Region, Date, Journal)]

#setorder(Xout,Trait)
gs4_create("IBS GWAS: GWAS Catalog Results",sheets=list(Xout))
#I did some massaging, renaming columns, adding e.g. links and author lists back to the final table

#Double check allele frequencies match for C/G variants:
library(ggplot2)
library(ggrepel)
Xout[,alleles:=paste0(`IBS Reference Allele`,"/",`IBS Effect Allele`)]
ggplot(Xout, aes(x=`IBS AF`,y=`Trait AF`)) +
  geom_abline(slope=1) +
  geom_label_repel(aes(label=paste0(snp,"\n",alleles), fill=alleles)) +
  labs(title="Allele frequencies between IBS sumstats and GWAS Catalog\nAfter flipping IBS for positive effect, reversing GWAS catalog, adding GWAS catalog reference allele, flipping GWAS catalog")
  
