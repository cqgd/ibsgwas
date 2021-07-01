#!/bin/bash
#SBATCH --job-name=ibs_anxdefs
#SBATCH --output=/gfs/work/ceijsbouts/ibs/jobs/stream/job_%A_%a.stdout
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --time=6-0
#SBATCH --mem-per-cpu=0
#SBATCH --cpus-per-task=24
#SBATCH --array=90-95
##SBATCH --exclude=kgen0001,kgen0002

##SBATCH --nodelist=kgen0010,kgen0011,kgen0012,kgen0013

#    1             2           3          4                        5            6             7               8              9              10
CCS=(final_cc.diag final_cc.sr final_cc.Q final_cc.Q.prev.diag.yes final_cc.any final_cc.male final_cc.female final_cc.sub.C final_cc.sub.U final_cc.sub.M
#  11             12               13               14                 15              16                    17                     18                        19                          20
 final_cc.sub.D final_quant.hard final_quant.loose final_quant.daily final_quant.SSS final_quant.daily.max final_quant.weekly.min final_cc.any.conts.relaxed final_cc.diag.conts.relaxed final_cc.sr.conts.relaxed
#  21                    22                                      23                    24                  25                     26                           27
 final_cc.Q.conts.relaxed final_cc.Q.prev.diag.yes.conts.relaxed final_cc.Qonly.Qbased final_cc.Qonly.any final_quant.daily.mauro final_quant.daily.max.mauro final_quant.weekly.min.mauro
# 28                 29                30           31
 final_cc.respondent final_cc.Qnon.any final_cc.all final_quant.phq.tired
# 32              #33               #34             #35              #36           #37            #38                          #39                           #40                    #41                    #42    $
final_cc.diag_ONE final_cc.diag_TWO final_cc.sr_ONE final_cc.sr_TWO final_cc.Q_ONE final_cc.Q_TWO final_cc.Q.prev.diag.yes_ONE final_cc.Q.prev.diag.yes_TWO final_cc.Qonly.any_ONE final_cc.Qonly.any_TWO final_cc$
#44                45                 46                 47                 48                49
final_cc.sub.any.C final_cc.sub.any.D final_cc.sub.any.M final_cc.sub.any.U final_cc.any.DvsC final_cc.DvsC
#50                      #51
final_cc.post.infective final_cc.nomauro.sub.any.C.female
#52                               53                                  54                          55
final_cc.nomauro.sub.any.C.female final_cc.likemauro.sub.any.C.female final_cc.nomauro.any.female final_cc.likemauro.any.female
#56                                    57                                     58                              59
final_cc.nobonfiglio.sub.any.C.female final_cc.likebonfiglio.sub.any.C.female final_cc.nobonfiglio.any.female final_cc.likebonfiglio.any.female
#60          #61                 62                     63
final_cc.anx final_cc.anx.no.ibs final_cc.Qonly.no.anx final_cc.Qnon.no.anx
#64           #65                  #66                 #67
final_cc.func.C final_cc.func.D final_quant.caseSSS final_cc.any.fhx
#68               #69                      #70
final_cc.anxbroad final_cc.anxbroad.no.ibs final_cc.no.anxbroad
#71                       #72                       #73
final_cc.Qonly.any.severe final_cc.Qonly.any.twoplus final_cc.Qnon.any.twoplus
#74                         #75                            #76                           #77
final_cc.Q.female.mauroreq final_cc.sub.C.female.mauroreq final_cc.sub.D.female.mauroreq final_cc.sub.M.female.mauroreq
#78                         #79                                      #80                           #81                              #82
final_cc.sr.female.mauroreq final_cc.Q.prev.diag.yes.female.mauroreq final_cc.diag.female.mauroreq final_cc.Qnon.sr.female.mauroreq final_cc.Qonly.sr.female.mauroreq
#83                                                #84                                             #85
final_cc.likebonfiglio.any.female.female.mauroreq final_cc.nobonfiglio.any.female.female.mauroreq final_cc.likebonfiglio.sub.any.C.female.female.mauroreq
#86                                                   #87                             #88                                                      #89
final_cc.nobonfiglio.sub.any.C.female.female.mauroreq final_cc.female.female.mauroreq final_cc.nobonfiglio.sub.any.C.female.female.mauroreq final_cc.female.female.mauroreq
#90                                                  #91                                           #92
final_cc.anxbroad_casesfrom_diag.interest_F40_or_F41 final_cc.anxbroad_casesfrom_diag.interest_F40 final_cc.anxbroad_casesfrom_diag.interest_F41
#93                                          #94                                         #95
final_cc.anxbroad_casesfrom_sr.interest_1287 final_cc.anxbroad_casesfrom_Q_treatment.anx final_cc.anxbroad_casesfrom_GAD7.anx
)

CC_ID=${CCS[${SLURM_ARRAY_TASK_ID}-1]}

mkdir /gfs/work/ceijsbouts/ibs/gwas/${CC_ID}

/gfs/apps/bio/BOLT-LMM-2.3.2/bolt \
--bed=/gfs/work/ceijsbouts/ibs/plink/ukb_chr{1:22}_geno_qc.bed \
--bim=/gfs/work/ceijsbouts/ibs/plink/ukb_chr{1:22}_geno_qc.bim \
--fam=/gfs/work/ceijsbouts/ibs/plink/ukb_chr22_geno_qc.fam \
--phenoFile=/gfs/work/ceijsbouts/ibs/covar/cc_cov.sample \
--phenoCol=${CC_ID} \
--covarFile=/gfs/work/ceijsbouts/ibs/covar/cc_cov.sample \
--qCovarCol=sex \
--qCovarCol=age \
--qCovarCol=agesq \
--qCovarCol=sexage \
--qCovarCol=sexagesq \
--qCovarCol=PC{1:20} \
--LDscoresFile=/gfs/apps/bio/BOLT-LMM-2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=/gfs/apps/bio/BOLT-LMM-2.3.2/tables/genetic_map_hg19_withX.txt.gz \
--bgenFile=/gfs/archive/jostins/ukbb/v3/imputation/ukb_imp_chr{1:22}_v3.bgen \
--sampleFile=/gfs/work/ceijsbouts/ibs/sample/ukb17670_imp_chr22_v3_s487327.sample \
--bgenMinMAF=0.01 \
--bgenMinINFO=0.3 \
--numThreads=24 \
--statsFile=/gfs/work/ceijsbouts/ibs/gwas/${CC_ID}/blmm_output.stats.gz \
--statsFileBgenSnps=/gfs/work/ceijsbouts/ibs/gwas/${CC_ID}/blmm_output.bgen.stats.gz \
--verboseStats \
--remove=/gfs/work/ceijsbouts/ibs/remove/bolt.in_plink_but_not_imputed.FID_IID.968.txt \
--noBgenIDcheck \
--lmmForceNonInf

