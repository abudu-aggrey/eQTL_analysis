#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00


# create genetic relatedness matrix


exp_sample_id="./TwinsUK/ecz_cases.txt"
number_of_samples=86
merged_geno="./TwinsUK/Eurobats_genotypes/All_chr/All_chr_genotype.bimbam"
temp_pheno="./TwinsUK/Eurobats_genotypes/All_chr/temp_pheno.txt"


#first subset individuals with exp data

module load apps/gtool-0.7.5

for i in {1..22}

do

gtool -S --g ./TwinsUK/Eurobats_genotypes/All_chr/gen/chr${i}_Eurobats_Public.gen \
--s ./TwinsUK/Eurobats_genotypes/All_chr/sample/chr${i}_Eurobats_Public.sample \
--og ./TwinsUK/Eurobats_genotypes/All_chr/gen/chr${i}_Eurobats_Public_subset.gen \
--os ./TwinsUK/Eurobats_genotypes/All_chr/sample/chr${i}_Eurobats_Public_subset.sample \
--sample_id ${exp_sample_id} \

done



#convert imputation .gen file to BIMBAM format (mean_genotype file):

for i in {1..22}

do

cat ./TwinsUK/Eurobats_genotypes/All_chr/gen/chr${i}_Eurobats_Public_subset.gen | awk -v \
s=${number_of_samples} \
'{ printf $2 "," $4 "," $5; for(i=1; i<=s; i++) \
printf "," $(i*3+3)*2+$(i*3+4); \
printf "\n" }' > ./TwinsUK/Eurobats_genotypes/All_chr/chr${i}_genotype.bimbam

done



#merge into one file

cat ./TwinsUK/Eurobats_genotypes/All_chr/chr*.bimbam > ${merged_geno}


#create fake pheno_file where there is no missing phenotype data (else those samples will not be included when computing the relatedness matrix)


for i in $(seq 1 $number_of_samples); do echo 1; done > ${temp_pheno}


#create relatedness matrix

gemma -g ${merged_geno} \
-p ${temp_pheno} \
-gk 1 -o relatedness_all_geno_cases

#copy to working directory

cp ./output/relatedness* ./TwinsUK/Eurobats_genotypes/

