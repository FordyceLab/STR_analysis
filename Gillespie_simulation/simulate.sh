#!/bin/bash
#
#SBATCH --job-name=simulation
#SBATCH --output=txt_output/%j.txt
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

echo $RUN_NUM
echo $TARGET
cd $TARGET
ml python/3.6.1
python3 ../simulate.py $RUN_NUM $TARGET
echo 'Job finished'
