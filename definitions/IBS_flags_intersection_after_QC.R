#install.packages("UpSetR")
#install.packages("nVennR")
library(UpSetR)
library(nVennR)
library(data.table)

setwd("~/work/ibs")

NAtoF = function(v){
  v[is.na(v)]=FALSE
  return(v)
}

#M = fread("/home/cq/work/ibs/filters/all_integrated.csv") #Old, from raw case definitions
M = fread("~/work/ibs/clinical/clinical_dataset.csv") #New, cc

IBS_flags = c("Hospital ICD-10"="final_case.diag",
              "Unprompted self-rep." = "final_case.sr",
              "DHQ Rome III"="final_case.Q",
              "DHQ self-report"="final_case.Q.prev.diag.yes"
              )
K = M[,lapply(.SD, NAtoF), .SDcols = c("Q_participant","Q_SSS",IBS_flags)]
K = K[,lapply(.SD, as.numeric)]

##Pairwise overlap for abstract
overlap = function(cols,x){
  both = x[get(cols[1]) & get(cols[2]),.N]
  union = x[get(cols[1]) | get(cols[2]),.N]
  return(round(100*both/union,1))
}

combos = expand.grid(IBS_flags,IBS_flags)

overlaps = cbind(combos, apply(combos,1,overlap,x=K[Q_participant==1]))
print(overlaps)
##

#UpSetR
pdf(file="plots/IBS_gwasdef_flag_overlap_UpSet.pdf",onefile=FALSE)
upset(K[,.SD,.SDcols=IBS_flags],text.scale=1.1)
dev.off()

pdf(file="plots/IBS_gwasdef_flag_overlap_UpSet_Q_participant_TRUE.pdf",onefile=FALSE)
upset(K[Q_participant==TRUE,.SD,.SDcols=IBS_flags],text.scale=1.1)
dev.off()

pdf(file="plots/IBS_gwasdef_flag_overlap_UpSet_Q_participant_FALSE.pdf",onefile=FALSE)
upset(K[Q_participant==FALSE,.SD,.SDcols=IBS_flags],text.scale=1.1)
dev.off()




#nVennR
plot_path = "plots/IBS_gwasdef_flag_overlap_nVennR_Q_participant.svg"
sets = lapply(IBS_flags,function(col,M){M[,which(get(col))]},M=M)
names(sets) <- names(IBS_flags)
plotVenn(sets,outFile=plot_path,labelRegions=FALSE)
png_cmd = paste("inkscape -z -e", paste0(plot_path,".png"), "-w 1400 -h 1000", plot_path)
system(png_cmd)

plot_path = "plots/IBS_gwasdef_flag_overlap_nVennR_Q_participant_TRUE.svg"
sets = lapply(IBS_flags,function(col,M){M[,which(get(col))]},M=M[Q_participant==TRUE])
names(sets) <- names(IBS_flags)
plotVenn(sets,outFile=plot_path,labelRegions=FALSE)
png_cmd = paste("inkscape -z -e", paste0(plot_path,".png"), "-w 1400 -h 1000", plot_path)
system(png_cmd)

plot_path = "plots/IBS_gwasdef_flag_overlap_nVennR_Q_participant_FALSE.svg"
sets = lapply(IBS_flags,function(col,M){M[,which(get(col))]},M=M[Q_participant==FALSE])
names(sets) <- names(IBS_flags)
plotVenn(sets,outFile=plot_path,labelRegions=FALSE)
png_cmd = paste("inkscape -z -e", paste0(plot_path,".png"), "-w 1400 -h 1000", plot_path)
system(png_cmd)



#SSS Distribution plot

#Z = M[Q_cases==1,.SD,.SDcols=c(IBS_flags,"Q_SSS")]
Z = M[final_case.Q==1,.SD,.SDcols=c(IBS_flags,"Q_SSS")]
Z = Z[,lapply(.SD,as.numeric)]
setcolorder(Z,IBS_flags)
setnames(Z,old="Q_SSS",new="SSS")
colnames(Z)[1:length(IBS_flags)]=names(IBS_flags)
Z[,SSS:=SSS*10]

trace(UpSetR:::BoxPlotsPlot,edit=TRUE)



function (bdat, att, att_color)
{
  yaxis <- as.character(att)
  col <- match(att, colnames(bdat))
  colnames(bdat)[col] <- "attribute"
  upper_xlim <- as.numeric((max(bdat$x) + 1))
  plot_lims <- factor(as.numeric(0:upper_xlim))
  print(plot_lims)
  bdat$x <- as.factor(bdat$x)
  bdat$severity <- as.factor(rowSums(bdat[, 1:4], na.rm = TRUE))
  print(head(bdat))
  fill_colors = c(`1` = "lightgoldenrod1", `2` = "lightgoldenrod2",
                  `3` = "lightgoldenrod3", `4` = "lightgoldenrod4")
  boxplots <- ggplotGrob(ggplot(bdat, aes(x = x, y = attribute,
                                          fill = severity)) + theme_bw() + ylab(yaxis) + scale_x_discrete(limits = plot_lims,
                                                                                                          expand = c(0, 0)) + theme(plot.margin = unit(c(-0.7,
                                                                                                                                                         0, 0, 0), "cm"), axis.title.y = element_text(vjust = -0.8),
                                                                                                                                    axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                                                                                                                                    panel.border = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                    panel.grid.major = element_blank(), axis.title.x = element_blank(),
                                                                                                                                    legend.position = "none") + geom_violin(draw_quantiles = c(0.25,
                                                                                                                                                                                               0.5, 0.75), colour = "grey80") + stat_summary(geom = "label",
                                                                                                                                                                                                                                             fun.y = median, aes(label = sprintf("%1.0f", ..y..)),
                                                                                                                                                                                                                                             color = "black", fill = "grey80", cex = 2.5) + scale_fill_manual(values = fill_colors))
  return(boxplots)
}


upset(Z,boxplot.summary = c("SSS"), text.scale=1.3, mb.ratio = c(1,1))
#The export is bugged so I just manually saved the above

pdf(file="plots/IBS_gwasdef_flag_overlap_UpSet_Rome_SSS2.pdf",onefile=FALSE)
x = upset(Z,boxplot.summary = c("SSS"), text.scale=1.3, mb.ratio = c(1,1))
dev.off()

upset(Z)




#Pairwise t-tests
Z[,set:=NULL]
Z[,set:=paste(`ICD-10`,`Self-reported`,`Quest. Prev. Diag.`,`Quest. Symptoms`,sep=",")] #determines order
ordered_levels = levels=c("0,0,0,1","0,0,1,1","0,1,0,1","1,0,0,1","0,1,1,1","1,0,1,1","1,1,0,1","1,1,1,1") #to match UpSet
ordered_names = c("QS","QS+QPD","QS+SR","QS+ICD","QS+QPD+SR","QS+QPD+ICD","QS+SR+ICD","QS+QPD+SR+ICD")
names(ordered_names) <- ordered_levels
Z[,set:=ordered_names[set]]
Z[,set:=factor(set,levels=ordered_names)]

lm(SSS~set,data=Z)
formal = pairwise.t.test(Z$SSS,Z$set,p.adjust.method = "bonferroni")
out = data.table(formal$p.value)
out[,set:=rownames(formal$p.value)]
setcolorder(out,"set")
fwrite(out,"tables/SSS_pairwise_T-tests.csv")
out

formal$p.value

##NOT USED
# 
# #UpSetR from Covar
# C = fread("/home/cq/work/ibs/covar/cc_cov.sample")
# ccs = grep("final_cc.",colnames(C),value=TRUE)
# upset(C[,.SD,.SDcols=ccs],text.scale=1.1)
# 
# #Pairwise
# G = data.table(expand.grid(colnames(K),colnames(K)))
# #G[,Var1:=as.character(Var1)][,Var2:=as.character(Var2)]
# 
# G[,value:=sum( K[[Var1]] & K[[Var2]] ), by=1:nrow(G)]
# dcast(G,"Var1~Var2")
# 
# #Overlapping groups, not used
# Za = melt(Z,id.vars = "SSS",value.name = "member",variable.name = "group")
# Za = Za[member==1]
# Za[,.N,by=group]
# 
# 
# #Trying to get nVennR to plot different labels, not used
# sets <- nVennR:::.flattenInput(sets, sNames = NULL)
# nBits <- nVennR:::.getNBits(sets)
# if (nBits == 0) {
#   stop("You must provide at least one list (and seriously consider providing more than one)")
# }
# lresult <- nVennR:::.processVenn(sets, nBits)
# lresult$def
# lresult$def[lresult$def=="19265"]="whatever"
# 
# myVenn = nVennR:::makeVenn(lresult,7000)
# showSVG(myVenn,outFile="~/work/ibs/plots/test.svg")
# 
# plotVenn
# nVennR:::refineVenn(myVenn)
# nVennR:::_nVennR_refineVenn
# #this is definitely it, how to inspect .call function?
# # function (x) 
# # {
# #   .Call("_nVennR_refineVenn", PACKAGE = "nVennR", x)
# # }
