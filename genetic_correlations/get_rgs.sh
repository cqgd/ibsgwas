mkdir logs
mv ./*.log -t logs

#SSAs=(*ldhub*.gz)
#SSBs=(*bgz*)

#SSAs=(*anx*.gz)
#SSBs=(*anx*.gz)

#SSAs=metal_ICD_ROME_EURUSA_Qonly_Qnon_ldhub.sumstats.sumstats.gz
#SSAs=metal_MAURO_ldhub.sumstats.sumstats.gz
#SSBs=(*ldhub*.gz)

#SSAs=(*_ONE_*)
#SSBs=(*_TWO_*)

#SSAs=(*ROME.IBS.META.fix-effect.for_share*)
#SSBs=(*final_cc.Q_ldhub*)
#SSAs=(*ldhub*.gz)
#SSBs=(*ldhub*.gz)

#SSAs=final_cc.Q_ldhub.sumstats.sumstats.gz
#SSBs=final_cc.respondent_ldhub.sumstats.sumstats.gz
#SSBs=(*bgz*)

#SSAs=(./subtypes_functional/*ldhub*.gz)
#SSBs=(./subtypes_functional/*ldhub*.gz)

SSAs=(*anxbroad*.gz)
SSBs=(*anxbroad*.gz)

echo ${SSAs[@]}
echo "AGAINST"
echo ${SSBs[@]}
echo "COMBOS:"

for SSA in "${SSAs[@]}"
do
        for SSB in "${SSBs[@]}"
        do
                echo ${SSA},${SSB}
        done
done

echo "START"

for SSA in "${SSAs[@]}"
do
        for SSB in "${SSBs[@]}"
        do
                echo ${SSA},${SSB}
                ~/devel/ldsc/ldsc.py \
                --rg ${SSA},${SSB} \
                --ref-ld-chr /gfs/devel/ceijsbouts/ldsc/tutorial/eur_w_ld_chr/ \
                --w-ld-chr /gfs/devel/ceijsbouts/ldsc/tutorial/eur_w_ld_chr/ \
                --out ~/work/ibs/rg/${SSA}___${SSB}
        done
done

echo "DONE."



