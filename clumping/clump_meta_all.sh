#!/bin/bash
#SBATCH --job-name=LDclump_meta
#SBATCH --output=/gfs/work/ceijsbouts/ibs/jobs/stream/job_%A_%a.stdout
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem-per-cpu=0
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22

##export OMP_NUM_THREADS=4

CHR=${SLURM_ARRAY_TASK_ID}
#      1              2            3                    4           5                                6                         7
METAS=(metal_ICD_diag metal_ROME_Q metal_Qonly_Qnon_any metal_MAURO metal_ICD_ROME_EURUSA_Qonly_Qnon metal_ICD_ROME_EURUSA_any metal_ICD_diag)

module load bio/plink/1.90b6.7

for m in {1..7}
do
META=${METAS[${m}-1]}
echo $META
echo $CHR

mkdir /gfs/work/ceijsbouts/ibs/clump/clumped/${META}

plink \
--bfile /gfs/archive/jostins/ukbb/v3/imputation_bfile_infomaf_filtered/ukb_imp_chr${CHR}_v3 \
--keep /gfs/work/ceijsbouts/ibs/clump/keep_unrelated/10k_passing_sqc.txt \
--clump /gfs/work/ceijsbouts/ibs/clump/sumstats_meta/${META}/chr${CHR}.sumstats \
--clump-field P-value \
--clump-snp-field MarkerName \
--clump-p1 5e-8 \
--clump-p2 0.05 \
--clump-r2 0.05 \
--clump-kb 5000 \
--clump-verbose \
--clump-allow-overlap \
--out /gfs/work/ceijsbouts/ibs/clump/clumped/${META}/chr${CHR}

done

#removed --clump-best

