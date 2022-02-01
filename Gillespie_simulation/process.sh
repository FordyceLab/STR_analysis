#!/bin/bash
#
#SBATCH --job-name=process
#SBATCH --output=txt_output/%j.txt
#SBATCH --time=10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

ml python/3.6.1
python3 process.py 100 n_flanks -log True
python3 process.py 100 core_affinity -log True
python3 process.py 100 koff_intercept -log True
python3 process.py 100 koff_slope -log True
python3 process.py 100 n_TF -log True
python3 process.py 100 switching_ratio -log True
python3 process.py 100 diffusion -log True
python3 process.py 100 DNA_concentration -log True
python3 process.py 100 local_volume -log True
python3 process.py 100 simulation
python3 process.py 100 simulation_mutated
python3 process.py 100 simulation_strong
cp n_flanks/n_* download/
cp core_affinity/cor* download/
cp koff_intercept/koff* download/
cp koff_slope/koff* download/
cp n_TF/n* download/
cp switching_ratio/switch* download/
cp diffusion/diff* download/
cp DNA_concentration/DNA* download/
cp local_volume/loc* download/
cp simulation/simulation_r* download/
cp simulation_mutated/simulation_m* download/
cp simulation_strong/simulation_s* download/
echo 'Job finished'
