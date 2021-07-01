#!/bin/bash
#SBATCH --job-name=bgen2bfile_infomaf
#SBATCH --output=/gfs/work/ceijsbouts/ibs/jobs/stream/job_%A_%a.stdout
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=0
#SBATCH --cpus-per-task=20
#SBATCH --array=4-22
##SBATCH --exclude=kgen[0001-00022]
##SBATCH --nodelist=kgen0501,kgen0502

CHR=${SLURM_ARRAY_TASK_ID}
/gfs/work/ceijsbouts/app/plink2 \
--bgen /gfs/archive/jostins/ukbb/v3/imputation/ukb_imp_chr${CHR}_v3.bgen \
--sample /gfs/work/ceijsbouts/ibs/sample/ukb17670_imp_chr22_v3_s487327.sample \
--maf 0.0001 \
--make-bed \
--mach-r2-filter 0.9 2.0 \
--out /gfs/archive/jostins/ukbb/v3/imputation_bfile_infomaf_filtered/ukb_imp_chr${CHR}_v3

