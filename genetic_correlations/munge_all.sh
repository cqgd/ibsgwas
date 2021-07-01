#SSes=(*ldhub*)
#SSes=(*ROME.IBS.META.fix-effect.for_share*)
#SSes=(*final_cc.Q_ldhub*)
#SSes=(*_ONE*)
#SSes=(*_TWO*)
#SSes=(*anx*)
#SSes=(*func*)
#SSes=(*anxbroad*)
#SSes=(*casesfrom*)
#SSes=(*20002_1287*)
SSes=(*mddformatted*)

for SS in "${SSes[@]}"
do
        echo ${SS}
        ~/devel/ldsc/munge_sumstats.py \
        --sumstats ${SS} \
        --out ~/work/ibs/munged/${SS} \
        --merge-alleles ~/devel/ldsc/tutorial/w_hm3.snplist
done



