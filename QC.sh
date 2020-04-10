#!/bin/bash
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

######################################################################
# plink.bed      ( binary file, genotype information )				 #
# plink.fam      ( first six columns of mydata.ped ) 				 #
# plink.bim      ( extended MAP file: two extra cols = allele names) #
######################################################################