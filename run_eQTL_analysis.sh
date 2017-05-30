#!/bin/bash
#
#
#
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00


CD207_snps="./TwinsUK/credible_SNP_list.txt"
genfile="./TwinsUK/Eurobats_genotypes/chr2_Eurobats_Public.gen"
sample_file="./TwinsUK/Eurobats_genotypes/chr2_Eurobats_Public_alt_pheno.sample"
outfile="./TwinsUK/Eurobats_genotypes/chr2_Eurobats_Public_subset_samples"
exp_sample_id="./TwinsUK/Sample_ID.txt"

number_of_samples=672
mean_geno="./TwinsUK/Eurobats_genotypes/chr2_genotype.bimbam"

new_sample_file="./TwinsUK/Eurobats_genotypes/chr2_Eurobats_Public_exp_pheno.sample"
pheno_file="./TwinsUK/Eurobats_genotypes/chr2_pheno.bimbam"
covar_file="./TwinsUK/Eurobats_genotypes/covar.txt"

snp_file="./TwinsUK/Eurobats_genotypes/SNP_annotation_file.txt"

eQTL_out="CD207_eQTL_analysis_test"

relatedness_matrix="./TwinsUK/Eurobats_genotypes/relatedness_all_geno.cXX.txt"

#******************************************************************************************************************

#to get list of SNPs for region (eg, CD207 at chr2p13.3= chr2:71025000-71150000)

#gawk '$3 >= 71025000 && $3 <= 71150000 {print $2}' < ${genfile} > ${CD207_snps}



#extract transcript data for CD207 (transcript id = ENSG00000116031)
#sed -n -e '1p' -e '/ENSG00000116031/p' \
#< /projects/TwinsUK_eczema/TwinsUK_PublicRNAseqData/Skin/EUROBATS.S.annonymousIDs.rpkm \
#> ./TwinsUK/CD207_transcript_data.txt

#read into R and merge with sample data, and create list of ID's with expression data

#module load languages/R-3.3.3-ATLAS

#Rscript ./scripts/get_expression_sample_id.R

#use sample id's to subset .gen and .sample imputation files

module load apps/gtool-0.7.5

gtool -S --g ${genfile} --s ${sample_file} --og ${outfile}.gen --os ${outfile}.sample \
--sample_id ${exp_sample_id} \
--inclusion ${CD207_snps}



#******************************************************************************************************************************************

#convert imputation .gen file to BIMBAM format (mean_genotype file):


cat ${outfile}.gen | awk -v s=${number_of_samples} \
'{ printf $2 "," $4 "," $5; for(i=1; i<=s; i++) \
printf "," $(i*3+3)*2+$(i*3+4); \
printf "\n" }' > ${mean_geno}

#**************************************************************************************************************************************

#create phenotype and covariate file: first merge the expression data


module load languages/R-3.3.3-ATLAS


Rscript ./scripts/merge_exp_data.R


#create phenotype file - extract the phenotype column from the sample file
#first column = case/control status (1 = case, 0 = controls)
#second column = expression data for CD207

tail -n +3 ${new_sample_file} | awk '{print $6, $7}' > ${pheno_file}

#create covariate file
#first column = intercept ("1"), second column = age

tail -n +3 ${new_sample_file} | awk '{print $9, $8}' > ${covar_file}

#**********************************************************************************************************************************

#create snp annotation file in R

module load languages/R-3.3.3-ATLAS

Rscript ./scripts/GEMMA_snp_file.R


#********************************************************************************************************************************

#create relatedness matrix using all genotype data with ./scripts/relatedness.sh
#will need to change sample_id's and no. of samples

#****************************************************************************************************************

#run eQTL analysis "-n 2" will tell it to look at the second column of the pheno file that has the expression data
# use "-n 1 2" to look at both expression and case-control status

./gemma -g ${mean_geno} -p ${pheno_file} -n 2 -a ${snp_file} -c ${covar_file} \
-lmm 4 -o ${eQTL_out} \
-maf 0 \
-k ${relatedness_matrix}

#copy output files to "TwinsUK directory"
cp ./output/${eQTL_out}.assoc.txt ./TwinsUK

#sort results file according to p-value
#may need to modify name of results file in script

module load languages/R-3.3.3-ATLAS
Rscript ./scripts/sort_eQTL_results.R

#***************************************************************************************************************


