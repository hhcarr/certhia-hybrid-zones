#move to working directory
cd /scratch/summit/hhcarr@xsede.org/certhia21/admix50kbpNTotal/

mkdir genomeScan

#run bcftools query to see sample names so they can be separated
bcftools query -l calls.vcf

# extract sample names for Black Hills
bcftools query -l calls.vcf | grep "BlackHills" > BH
# extract sample names for East
bcftools query -l calls.vcf | grep 'NCarolina\|WVirginia' > East
# extract sample names for Oregon
bcftools query -l calls.vcf | grep "Oregon" > Ore
# extract sample names for California
bcftools query -l calls.vcf | grep "California" > Cal
# extract sample names for Colorado
bcftools query -l calls.vcf | grep "Colorado" > Colo
# extract sample names for Nevada
bcftools query -l calls.vcf | grep "Nevada" > Nev
# extract sample names for Utah
bcftools query -l calls.vcf | grep "Utah" > Utah
# extract sample names for Pacific
bcftools query -l calls.vcf | grep 'Oregon\|California' > Pac
# extract sample names for Rocky Mountain
bcftools query -l calls.vcf | grep 'Colorado\|Utah' > RM

#start interactive session
sinteractive --partition=shas --time=03:00:00 --nodes=1 --ntasks=4
conda activate birbenv

#calculate mean Fst accross all SNPs in the vcf file (whole genome in this case)
vcftools --vcf calls.vcf --weir-fst-pop BH --weir-fst-pop East --out ./BH_East
vcftools --vcf calls.vcf --weir-fst-pop Ore --weir-fst-pop East --out genomeScan/Ore_East
vcftools --vcf calls.vcf --weir-fst-pop Cal --weir-fst-pop East --out genomeScan/Cal_East
vcftools --vcf calls.vcf --weir-fst-pop Utah --weir-fst-pop East --out genomeScan/Utah_East
vcftools --vcf calls.vcf --weir-fst-pop Pac --weir-fst-pop BH --out genomeScan/Pac_BH
vcftools --vcf calls.vcf --weir-fst-pop RM --weir-fst-pop BH --out genomeScan/RM_BH

