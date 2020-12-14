#!/bin/bash
#$ -cwd
#$ -V
#$ -m bea
#$ -l h_data=10G
#$ -l h_rt=01:00:00

# INPUT ARGUMENTS
s_value=s_0_rec

# given s value
python3 compute_statistic.py ../Data/${s_value}/sample_10.vcf
python3 compute_statistic.py ../Data/${s_value}/sample_50.vcf
python3 compute_statistic.py ../Data/${s_value}/sample_100.vcf
python3 compute_statistic.py ../Data/${s_value}/sample_500.vcf
python3 compute_statistic.py ../Data/${s_value}/sample_1000.vcf
python3 compute_statistic.py ../Data/${s_value}/sample_5000.vcf
python3 compute_statistic.py ../Data/${s_value}/sample_10000.vcf

# python3 compute_HW_expectation.py 10000 $SGE_TASK_ID

# Hardy Weinberg Expectation
# for i in {2..200}
# do
#   popsize=${i/2+1}
#   python compute_HW_expectation.py ${popsize} ${i}
# done

rm -rf *.e*
rm -rf *.o*
