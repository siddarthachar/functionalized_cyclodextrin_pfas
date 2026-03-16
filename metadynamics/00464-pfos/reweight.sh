shopt -s expand_aliases
source ~/.bashrc
prepgmx
export PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/bin:$PATH"
export LD_LIBRARY_PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/lib:$LD_LIBRARY_PATH"

rm -r reweight
mkdir reweight
cd reweight

cp ../plumed.dat reweight.dat
cp ../HILLS .
cp ../hills* .
cp ../COLVAR .
cp ../BIAS .
find ../ -type f -name "*.pdb" ! -name "*step*" -exec cp {} . \;


sed s/"COLVAR"/"COLVAR_reweight"/g reweight.dat > temp.dat
sed s/"PACE=500"/"PACE=500000000"/g temp.dat > temp1.dat
sed s/"HEIGHT=1.2"/"HEIGHT=0"/g temp1.dat > temp2.dat
sed s/"STRIDE=100"/"STRIDE=1"/g temp2.dat > temp3.dat
mv temp3.dat reweight.dat
rm -f temp.dat temp1.dat temp2.dat

gmx editconf -f ../sys-emin.gro -o sys-emin.pdb

fileName=sys-prod
#gmx mdrun -ntomp $OMP_NUM_THREADS -pin on -plumed reweight.dat -s ../$fileName.tpr -rerun ../$fileName.xtc
plumed driver --plumed reweight.dat --mf_xtc ../$fileName.xtc --trajectory-stride 0 --timestep 0.002 --pdb sys-emin.pdb
