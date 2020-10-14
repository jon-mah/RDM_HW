#!/bin/bash
#$ -cwd
#$ -V
#$ -l h_data=1G
#$ -l h_rt=0:10:00
#$ -l h_vmem=50G

# INPUT ARGUMENTS
seed=1
num_ind=8 # Number of samples from single population.
easySFS_proj=16 # Number of chromosomes that sample is projected down into.
s_value=s_0_001
sample_size=10
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..100}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done

# DERIVED ARGUMENTS
inputvcf="${prefix}.vcf"
syn_inputvcf="${prefix}_syn.vcf"
nonsyn_inputvcf="${prefix}_nonsyn.vcf"
syn_output_pops_file="${prefix}_syn_pops_file.txt"
nonsyn_output_pops_file="${prefix}_nonsyn_pops_file.txt"
syn_easySFS_outdir="${prefix}_easySFS_output_syn/"
nonsyn_easySFS_outdir="${prefix}_easySFS_output_nonsyn/"
num_samples=$(($num_ind * 2))

# Separate input `.vcf` file into synonymous and nonsynonymous.
# grep "#\|;MT=1" ${prefix}.vcf > $syn_inputvcf
# grep "#\|;MT=2" ${prefix}.vcf > $nonsyn_inputvcf
# Write pops_file.txt for easySFS.py.
# python write_pops_file.py $num_ind 0 ${prefix}_syn
# python write_pops_file.py $num_ind 0 ${prefix}_nonsyn
# Use easySFS.py to compute synonymous and nonsynonymous sfs.
# python easySFS.py -i $syn_inputvcf -p $syn_output_pops_file -o $syn_easySFS_outdir -f -a --proj $easySFS_proj
# python easySFS.py -i $nonsyn_inputvcf -p $nonsyn_output_pops_file -o $nonsyn_easySFS_outdir -f -a --proj $easySFS_proj

# INPUT ARGUMENTS
sample_size=50
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..100}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done

# INPUT ARGUMENTS
sample_size=100
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..100}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done

# INPUT ARGUMENTS
sample_size=500
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..100}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done

# INPUT ARGUMENTS
sample_size=1000
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..100}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done
