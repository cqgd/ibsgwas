library(data.table)
setwd("~/work/ibs")


tissues_dir = "~/archive/gtex/GTEx_Analysis_v7_eQTL_all_associations"

tissues = gsub(".allpairs.txt$","",list.files(tissues_dir,pattern=".allpairs.txt$"))

exclude=""
tissues = tissues[!(tissues %in% exclude)]

print(tissues)

save_gene = function(X,path_prefix){
  path = paste0(path_prefix,"/",X$gene_id[1],".csv")
  fwrite(X,path)
}

for(t in seq_along(tissues)){
  regions = fread(paste0("resource/gtex/regions/",tissues[t],".csv"))

  G = fread(paste0(tissues_dir,"/",tissues[t],".allpairs.txt"))
  G = G[gene_id%in%unique(regions$gene_id)]

  R = fread("resource/gtex/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt")

  GR = merge(G,R[,.(variant_id,rsid=rs_id_dbSNP147_GRCh37p13)],by="variant_id",all.x=TRUE,all.y=FALSE)

  path_prefix = paste0("resource/gtex/tissue_all_relevant_regions/",tissues[t],"/")
  dir.create(path_prefix,showWarnings = FALSE, recursive=TRUE)

  returns = by(GR,GR$gene_id,save_gene,path_prefix)

  print(tissues[t])
  print(timestamp())
}



