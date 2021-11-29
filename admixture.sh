# make chromosome map for the vcf
grep -v "#" calls50.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt

# run vcftools for the combined vcf
vcftools --vcf calls50.vcf  --plink --chrom-map chrom_map.txt --out total 

# convert  with plink
plink --file total --allow-extra-chr --recode 12 --out total2

# run admixture 
for K in 1 2 3 4 5 6 7; do admixture --cv total200.ped $K  | tee log2_${K}.out; done

# check cv
grep -h CV log2_*.out


# 50kbp
# CV error (K=1): 0.39686
# CV error (K=2): 0.39858
# CV error (K=3): 0.42274
# CV error (K=4): 0.45843
# CV error (K=5): 0.54013
# CV error (K=6): 0.57631
# CV error (K=7): 0.63481

# 100kbp
# CV error (K=1): 0.39306
# CV error (K=2): 0.39250
# CV error (K=3): 0.43408
# CV error (K=4): 0.46698
# CV error (K=5): 0.52042
# CV error (K=6): 0.58240
# CV error (K=7): 0.61936