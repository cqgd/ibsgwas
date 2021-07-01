#!/bin/bash
#SBATCH --job-name=cojo_meta
#SBATCH --output=/gfs/work/ceijsbouts/ibs/jobs/stream/job_%A_%a.stdout
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem-per-cpu=0
#SBATCH --cpus-per-task=4
#SBATCH --array=1-5

TRAIT="metal_ICD_ROME_EURUSA_Qonly_Qnon"
extract_dir="/gfs/work/ceijsbouts/ibs/clump/extract"

loci_files=($extract_dir/$TRAIT/chr*)

echo ${loci_files[@]}

locus=${SLURM_ARRAY_TASK_ID}
locus_path=${loci_files[$locus-1]}
locus_name=$(basename "$locus_path")

CHR=$(echo $locus_path | sed 's/.*chr\([0-9]*\).*/\1/')

printf "Working with $TRAIT locus $locus \n from $locus_path \n on chr$CHR"

out_dir="/gfs/work/ceijsbouts/ibs/clump/cojo/$TRAIT"
mkdir ${out_dir}
out_file=${out_dir}/${locus_name}


module load bio/gcta/1.92.0b

gcta64 \
        --bfile /gfs/archive/jostins/ukbb/v3/imputation_bfile_infomaf_filtered/ukb_imp_chr${CHR}_v3 \
        --keep /gfs/work/ceijsbouts/ibs/clump/keep_unrelated/10k_passing_sqc.txt \
        --extract $locus_path \
        --cojo-file /gfs/work/ceijsbouts/ibs/clump/sumstats_meta_cojo/${TRAIT}/chr${CHR}.sumstats \
        --cojo-slct \
        --cojo-p 5e-8 \
        --thread-num 4 \
        --prevalence 0.15 \
        --out ${out_file}

