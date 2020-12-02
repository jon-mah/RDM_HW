#!/bin/bash
#$ -cwd
#$ -V
#$ -m bea
#$ -l h_data=10G
#$ -l h_rt=1:00:00

# INPUT ARGUMENTS
s_value=s_0_001_rec

# s_0_rec
python compute_HW_departure.py ../Data/${s_value}/sample_10.vcf
python compute_HW_departure.py ../Data/${s_value}/sample_50.vcf
python compute_HW_departure.py ../Data/${s_value}/sample_100.vcf
python compute_HW_departure.py ../Data/${s_value}/sample_500.vcf
python compute_HW_departure.py ../Data/${s_value}/sample_1000.vcf
# python compute_HW_departure.py ../Data/${s_value}/sample_5000.vcf
# python compute_HW_departure.py ../Data/${s_value}/sample_10000.vcf
