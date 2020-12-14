#!/bin/bash
#$ -cwd
#$ -t 2-1000:1
#$ -V
#$ -l h_data=10G
#$ -l h_rt=00:05:00

# INPUT ARGUMENTS
s_value=s_0_rec

# s_0_rec
# python compute_HW_departure.py ../Data/${s_value}/sample_10.vcf
# python compute_HW_departure.py ../Data/${s_value}/sample_50.vcf
# python compute_HW_departure.py ../Data/${s_value}/sample_100.vcf
# python compute_HW_departure.py ../Data/${s_value}/sample_500.vcf
# python compute_HW_departure.py ../Data/${s_value}/sample_1000.vcf
# python compute_HW_departure.py ../Data/${s_value}/sample_5000.vcf
# python compute_HW_departure.py ../Data/${s_value}/sample_10000.vcf

python3 compute_HW_expectation.py 10000 $SGE_TASK_ID

# Hardy Weinberg Expectation
# for i in {2..200}
# do
#   popsize=${i/2+1}
#   python compute_HW_expectation.py ${popsize} ${i}
# done

rm -rf *.e*
rm -rf *.o*
