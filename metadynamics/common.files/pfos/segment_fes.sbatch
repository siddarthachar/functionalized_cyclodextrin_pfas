#!/bin/bash

# email on start, end, and abortion
#SBATCH --job-name=rewight
##SBATCH --output=out.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sachar@rcc.uchicago.edu
#SBATCH --partition=andrewferguson
#SBATCH --account=pi-andrewferguson
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0 # use max memory on node
#SBATCH --time=0:30:00
##SBATCH --exclude=midway2-0661

bash segment_fes.sh