library(ggplot2)
library(cowplot)
library(data.table)
theme_set(theme_cowplot())
setwd("~/work/ibs")
R = readRDS("tables/R_for_natgen_rg_plot.rds")
data_combos = c("Discovery cohort","Bellygenes data only","UK Biobank data only","Higher-specificity","Severe (IBS-SSS>300)")
R = R[plot==TRUE & trait1 %in% data_combos]
R[,trait1 := factor(trait1,levels=data_combos)]
manual_traits = c("Anxiety or panic attacks","Asthma","Ulcerative colitis","Crohn's disease")
p.rg = ggplot(R, aes(x=rg,y=trait1,col=trait1)) +
  geom_vline(xintercept=0,lty="solid",color="grey60") +
  geom_point() + 
  geom_errorbarh(aes(xmin=ci_lower, xmax=ci_upper)) + 
  geom_rect(data = R[trait2%in%manual_traits & trait1=="UK Biobank data only"],fill="goldenrod1",color=NA,xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.1, show.legend=FALSE) +
  #geom_text(aes(x=-1,label=p_label)) +
  xlim(-1.09,1.09) +
  #geom_vline(xintercept=0,col="lightgrey") +
  background_grid() +
  ylab("Trait compared to IBS") + xlab("Genetic correlation") +
  labs(col="IBS definition") +
  #facet_grid(t2f~.,scales = "free",space="free_y",switch="both") +
  #facet_grid(t2f~.,scales = "free",space="free_y") +
  facet_wrap(~t2f,strip.position="left",ncol=1) +
  theme(strip.text.y.left = element_text(angle = 0), #this is somehow not working #now needs .left https://github.com/tidyverse/ggplot2/issues/3888
        #strip.text = element_text(angle = 70),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        legend.position="right") +
  scale_y_discrete(limits = rev(levels(R$trait1)))
p.rg
ggsave("~/work/ibs/plots/natgen_rgs.pdf",p.rg,width=10,height=7)
fwrite(R[plot==TRUE],"~/work/ibs/tables/natgen_rgs.csv")
