#!/bin/sh
# email on start, end, and abortion
#SBATCH --job-name=pfos-bcd
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

# # energy minimization
gmx grompp -f emin.mdp -c sys-neutral.gro -p sys.top -po sys-out.mdp -o sys.tpr -maxwarn 2 -r sys-neutral.gro
gmx mdrun -s sys.tpr -deffnm sys-emin

# equilibration
gmx grompp -f equil.mdp -c sys-emin.gro -p sys.top -po sys-equil-out.mdp -o sys-equil.tpr -maxwarn 2 -r sys-neutral.gro
gmx mdrun -s sys-equil.tpr -deffnm sys-equil

# NPT equilibration
gmx grompp -f equil-npt.mdp -c sys-equil.gro -p sys.top -po sys-equil-npt-out.mdp -o sys-equil-npt.tpr -maxwarn 2 -r sys-neutral.gro
gmx mdrun -s sys-equil-npt.tpr -deffnm sys-equil-npt

## Production run
gmx grompp -f prod.mdp -c  sys-equil-npt.gro -p sys.top -r sys-equil-npt.gro -po sys-prod-out.mdp -o sys-prod.tpr -maxwarn 2
export OMP_NUM_THREADS=16
gmx mdrun -ntomp $OMP_NUM_THREADS -pin on -s sys-prod.tpr -deffnm sys-prod -nice 1 -plumed plumed.dat
