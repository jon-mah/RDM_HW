#!/bin/bash
#$ -t 1:1
#$ -cwd
#$ -V
#$ -m a
#$ -l h_data=10G
#$ -l h_rt=0:30:00

# INPUT ARGUMENTS
seed=1

# If `SLiM` is not an executable, then comment out the next line.
slim -d chrom=$SGE_TASK_ID -d init_seed=$seed recessive_deleterious_mutations.slim

# If `SLiM` is not an executable, then uncomment the next line and
# provide a path to `SLiM`.
# ./slim -d chrom=$SGE_TASK_ID -d init_seed=$seed AW_to_neutral_simulation.slim
