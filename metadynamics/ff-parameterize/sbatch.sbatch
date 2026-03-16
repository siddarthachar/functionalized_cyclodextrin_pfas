#!/bin/bash

# email on start, end, and abortion
#SBATCH --job-name=NAME
##SBATCH --output=out.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sachar@uchicago.edu
#SBATCH --partition=caslake
#SBATCH --account=pi-andrewferguson
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
##SBATCH --mem=0 # use max memory on node
#SBATCH --time=30:00:00
##SBATCH --exclude=midway2-0661

# THIS EXAMPLE USES 1 GPU NODE - 1 MPI TASK - 4 THREADS PER TASK

# SET NUMBER OF MPI TASKS 
# SET NUMBER OF MD STEPS

#LOAD GROMACS MODULE 


#module load python
py # load python anaconda 2023
#module load python/anaconda-2020.02
module load amber/24
module load gaussian/16RevA.03
#module load openbabel/2.4

source activate gaussian

bash gen_top.sh
