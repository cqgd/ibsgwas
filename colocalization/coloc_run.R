library(data.table)
setwd("~/work/ibs")
library(coloc)

meta_dataset = function(meta_path){
  SS = fread(meta_path)
  SS = SS[nchar(SNP)>=1]
  s = median(SS$N_CASE)/median(SS$N_TOTAL)
  SS = SS[,.(snp=SNP, pvalues=`P-value`, beta=Effect, varbeta=StdErr^2)] #Should MAF be min(MAF,1-MAF)?
  dataset = list(type="cc",
                 s=s,
                 #pvalues=SS$pvalues,
                 beta=SS$beta,
                 varbeta=SS$varbeta,
                 snp=SS$snp)
  return(dataset)
}

gtex_gene_dataset = function(gtex_gene_path){
  G = fread(gtex_gene_path)
  G = G[,.(variant_id, beta=slope, varbeta=slope_se^2, rsid, maf)]
  G = G[!is.na(beta) & !is.na(varbeta)]

  dataset = list(type="quant",
                 #pvalues=SS$pvalues,
                 beta=G$beta,
                 varbeta=G$varbeta,
                 snp=G$rsid,
                 #MAF=G$maf, N=500)
                 sdY=1) #GTex uses normalized data
  return(dataset)
}

#meta_path = "meta/runs/metal_ICD_ROME_EURUSA_Qonly_Qnon/METAANALYSIS1.TBL.pos"
meta_path = "/gfs/work/avoda/ibs/METAANALYSIS1.TBL.pos" #Alex's copy
META = meta_dataset(meta_path)

regions_dir = "resource/gtex/regions"
regions_files = list.files(regions_dir,pattern=".csv$")
regions_paths = paste0(regions_dir,"/",regions_files)

tissues = gsub("*.csv","",regions_files)

t = grep("Colon_Sigmoid",tissues)

tissue_regions_dir = "resource/gtex/tissue_all_relevant_regions/"
tissues = intersect(tissues,list.dirs(tissue_regions_dir,full.names=FALSE)) #check if you don't lose any tissues you have regions for
print(tissues)

excluded_tissues = ""
#excluded_tissues = c("Colon_Sigmoid","Colon_Transverse")
tissues = tissues[!tissues%in%excluded_tissues]

regions_paths = paste0(regions_dir,"/",tissues,".csv") #Subset the above via table, don't regenrate this


t_range = seq_along(tissues)

#t_range = 23:24
print(tissues[t_range])

for(t in t_range){
  print(tissues[t])
  TR = fread(regions_paths[t])

  TR = TR[,gtex_gene_path:=paste0(tissue_regions_dir,"/",tissue,"/",gene_id,".csv")]
  for(r in 1:nrow(TR)){
    print(TR[r])
    if(file.exists(TR[r,gtex_gene_path])){
      G = gtex_gene_dataset(TR[r,gtex_gene_path])
      coloc = coloc.abf(META,G)
      coloc_tissue_dir = paste0("coloc/",TR[r,tissue])
      coloc_path = paste0(coloc_tissue_dir,"/",TR[r,gene_id],".RDS")
      dir.create(coloc_tissue_dir,showWarnings = FALSE)
      saveRDS(coloc,coloc_path)
      print(coloc_path)
    }else{
      print(paste0(TR[r,gtex_gene_path]," doesn't exist."))
    }
  }
}


