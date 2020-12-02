#!/bin/bash
#$ -cwd
#$ -V
#$ -l h_data=1G
#$ -l h_rt=0:10:00
#$ -l h_vmem=50G

# INPUT ARGUMENTS
s_value=s_0_001_rec
sample_size=10
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..400}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done

# INPUT ARGUMENTS
sample_size=50
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..400}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done

# INPUT ARGUMENTS
sample_size=100
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..400}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done

# INPUT ARGUMENTS
sample_size=500
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..400}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done

# INPUT ARGUMENTS
sample_size=1000
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..400}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done

# INPUT ARGUMENTS
sample_size=5000
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..400}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done

# INPUT ARGUMENTS
sample_size=10000
prefix="../Data/${s_value}/sample_${sample_size}" # Output prefix, and input prefix of given vcf

# Merge all `.vcf` files.

cp ${prefix}_chrom_1.vcf ${prefix}.vcf

for i in {2..400}
do
  sed -i -r "s/^1\t/${i}\t/g" ${prefix}_chrom_${i}.vcf
  grep ";MT=1\|;MT=2" ${prefix}_chrom_${i}.vcf >> ${prefix}.vcf
done
