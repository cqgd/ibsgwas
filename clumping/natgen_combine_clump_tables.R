library(data.table)
setwd("~/work/ibs")
replication_dir = "clump/replication"
files = list.files(replication_dir,pattern="^(replication_)")
exclude = grepl("(^replication_for)|(^replication_counts)",files)
files = files[!exclude]
print(files)
paths = paste0(replication_dir,"/",files)
all = rbindlist(lapply(paths,fread))
setcolorder(all,"trait")
overview_path = paste0(replication_dir,"/overview_natgen.csv")
fwrite(all,overview_path)
print(overview_path)

#Same as before, but adding the severe IBS hit