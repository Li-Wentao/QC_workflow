#!/bin/bash
### QC ###
	### Step 1: check missingness ###
	plink --bfile $1 --missing   													#check for missing data
	plink --bfile $1 --geno 0.02 --make-bed --out data_out				#setting the threshold 0.02, can be changed accordingly
	plink --bfile $1 --mind 0.02 --make-bed --out data_out				#setting the threshold 0.02, can be changed accordingly

	### Step 2: check sex-discrepancy ###
	plink --bfile data_out --check-sex												#check for sex discrepancies
	grep "PROBLEM" plink.sexcheck | awk '{print$1,$2}' > sex_discrepancy.txt		#locate the problem data
	plink --bfile data_out --remove sex_discrepancy.txt --make-bed --out data_out		#remove the problem data

	### Step 3: check MAF ###
	awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' data_out.bim > snp.txt				#select autosomal SNPs only
	plink --bfile data_out --extract snp.txt --make-bed --out data_out					
	plink --bfile data_out --freq --out MAF_check										#check MAF
	plink --bfile data_out --maf 0.05 --make-bed --out data_out							#remove snps with low MAF

	### Step 4: check Hardy-Weinberg equilibrium (HWE) ###
	plink --bfile data_out --Hardy 													#check Hardy-Weinberg equilibrium
	plink --bfile data_out --hwe 1e-6 --make-bed --out HWE_controls					#HWE threshold for controls
	plink --bfile data_out --hwe 1e-10 --hwe-all --make-bed --out HWE_case			#HWE threshold for case

	### Step 5: check heterozygosity ###
	plink --bfile HWE_case --indep-pairwise 50 5 0.2 --out indep_snp 				#find out the independent snps
	plink --bfile HWE_case --extract indep_snp.prune.in --het --out R_check
	Rscript --no-save check_heterozygosity_rate.R 									#plot the hist-gram of heterozygosity rate 
	Rscript --no-save heterozygosity_outliers_list.R 								#find out the high and low heterozygosity rate
	sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}' > het_fail_ind.txt
	plink --bfile HWE_case --remove het_fail_ind.txt --make-bed --out data_out

	### Step 6: check relatedness ###
	plink --bfile data_out --extract indep_snp.prune.in --genome --min 0.2 --out pihat_min0.2
	awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome > zoom_pihat.genome
	plink --bfile data_out --filter-founders --make-bed --out data_out					#remove data from relatedness

### Population Stratification ###
	# Check whether you have 1000 Genome data to render Population Stratification
	if [ -e ALL.2of4intersection.20100804.genotypes.vcf.gz ]
	then
	    echo "You have already downloaded the 1000 Genome data"
	else
	    echo "You haven't already downloaded the 1000 Genome data, and start donwloading"
		wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz
	fi
	# unzip the vcf file into binary files
	plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --make-bed --out ALL.2of4intersection.20100804.genotypes
	### QC for 1000 Genome data ###
	plink --bfile ALL.2of4intersection.20100804.genotypes --set-missing-var-ids @:#\$1,\$2 --make-bed --out ALL.2of4intersection.20100804.genotypes_no_missing_IDs
	plink --bfile ALL.2of4intersection.20100804.genotypes_no_missing_IDs --geno 0.2 --allow-no-sex --make-bed --out 1KG_MDS
	plink --bfile 1KG_MDS --mind 0.2 --allow-no-sex --make-bed --out 1KG_MDS
	plink --bfile 1KG_MDS --maf 0.05 --allow-no-sex --make-bed --out 1KG_MDS

	### Extract the variants present in HapMap dataset from the 1000 genomes dataset ###
	awk '{print$2}' data_out.bim > HapMap_SNPs.txt
	plink --bfile 1KG_MDS --extract HapMap_SNPs.txt --make-bed --out 1KG_MDS

	### Extract the variants present in 1000 Genomes dataset from the HapMap dataset ###
	awk '{print$2}' 1KG_MDS.bim > 1KG_MDS_SNPs.txt
	plink --bfile data_out --extract 1KG_MDS_SNPs.txt --recode --make-bed --out HapMap_MDS

	### The datasets must have the same build. Change the build 1000 Genomes data build ###
	awk '{print$2,$4}' HapMap_MDS.map > buildhapmap.txt
	plink --bfile 1KG_MDS --update-map buildhapmap.txt --make-bed --out 1KG_MDS

	### Merge the HapMap and 1000 Genomes data sets ###
	# 1) set reference genome
	awk '{print$2,$4}' 1kG_MDS.bim > 1kg_ref-list.txt
	plink --bfile HapMap_MDS --reference-allele 1kg_ref-list.txt --make-bed --out HapMap-adj
	# 2) Resolve strand issues
	awk '{print$2,$5,$6}' 1KG_MDS.bim > 1KGMDS_tmp
	awk '{print$2,$5,$6}' HapMap-adj.bim > HapMap-adj_tmp
	sort 1KGMDS_tmp HapMap-adj_tmp |uniq -u > all_differences.txt
	awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
	plink --bfile HapMap-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_hapmap
	awk '{print$2,$5,$6}' corrected_hapmap.bim > corrected_hapmap_tmp
	sort 1kGMDS_tmp corrected_hapmap_tmp |uniq -u  > uncorresponding_SNPs.txt
	# 3) Remove problematic SNPs from HapMap and 1000 Genomes
	awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
	plink --bfile corrected_hapmap --exclude SNPs_for_exlusion.txt --make-bed --out HapMap_MDS
	plink --bfile 1KG_MDS --exclude SNPs_for_exlusion.txt --make-bed --out 1KG_MDS
	plink --bfile HapMap_MDS --bmerge 1kG_MDS.bed 1KG_MDS.bim 1kG_MDS.fam --allow-no-sex --make-bed --out MDS_merge
	plink --bfile MDS_merge --extract indep_snp.prune.in --genome --out MDS_merge
	plink --bfile MDS_merge --read-genome MDS_merge.genome --cluster --mds-plot 10 --out MDS_merge

	### Exclude ethnic outliers ###
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/20100804.ALL.panel
	# Convert population codes into superpopulation codes
	awk '{print$1,$1,$2}' 20100804.ALL.panel > race_1kG.txt
	sed 's/JPT/ASN/g' race_1kG.txt>race_1kG2.txt
	sed 's/ASW/AFR/g' race_1kG2.txt>race_1kG3.txt
	sed 's/CEU/EUR/g' race_1kG3.txt>race_1kG4.txt
	sed 's/CHB/ASN/g' race_1kG4.txt>race_1kG5.txt
	sed 's/CHD/ASN/g' race_1kG5.txt>race_1kG6.txt
	sed 's/YRI/AFR/g' race_1kG6.txt>race_1kG7.txt
	sed 's/LWK/AFR/g' race_1kG7.txt>race_1kG8.txt
	sed 's/TSI/EUR/g' race_1kG8.txt>race_1kG9.txt
	sed 's/MXL/AMR/g' race_1kG9.txt>race_1kG10.txt
	sed 's/GBR/EUR/g' race_1kG10.txt>race_1kG11.txt
	sed 's/FIN/EUR/g' race_1kG11.txt>race_1kG12.txt
	sed 's/CHS/ASN/g' race_1kG12.txt>race_1kG13.txt
	sed 's/PUR/AMR/g' race_1kG13.txt>race_1kG14.txt
	awk '{print$1,$2,"OWN"}' HapMap_MDS.fam>racefile_own.txt
	cat race_1kG14.txt racefile_own.txt | gsed -e '1i\FID IID race'  > racefile.txt
	Rscript MDS_merged.R 
legend("topright", pch=c(1,1,1,1,3),c("EUR","ASN","AMR","AFR","OWN"),col=c("green","red",470,"blue","black"),bty="o",cex=1)'
	### Exclude ethnic outliers ###
	awk '{ if ($4 <-0.03 && $5 >0.03) print $1,$2 }' MDS_merge.mds > EUR_MDS_merge
	plink --bfile data_out --keep EUR_MDS_merge --make-bed --out HapMap_out
	plink --bfile HapMap_out --extract indep_snp.prune.in --genome --out HapMap_out
	plink --bfile HapMap_out --read-genome HapMap_out.genome --cluster --mds-plot 10 --out HapMap_out_mds
	awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}'  HapMap_out_mds.mds > covar_mds.txt

### Association analysis ###
	# chi-squared test
	plink --bfile HapMap_out --assoc --out assoc_result
	# logistic
	plink --bfile HapMap_out --covar covar_mds.txt --logistic --hide-covar --out logistic_result 
	awk '!/'NA'/' logistic_result.assoc.logistic > logistic_result2.assoc.logistic




















