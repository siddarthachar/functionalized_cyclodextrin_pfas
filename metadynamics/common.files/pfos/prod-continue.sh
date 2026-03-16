#!/bin/sh
# email on start, end, and abortion
#SBATCH --job-name=prod
##SBATCH --output=out.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sachar@rcc.uchicago.edu
#SBATCH --partition=gpu #andrewferguson-gpu
#SBATCH --account=pi-andrewferguson
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=20
#SBATCH --time=30:00:00

shopt -s expand_aliases
source ~/.bashrc
prepgmx
export PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/bin:$PATH"
export LD_LIBRARY_PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/lib:$LD_LIBRARY_PATH"

## Production run continue
OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
gmx mdrun -ntomp $OMP_NUM_THREADS -pin on -s sys-prod.tpr -deffnm sys-prod -nice 1 -plumed plumed.dat -cpi sys-prod.cpt
