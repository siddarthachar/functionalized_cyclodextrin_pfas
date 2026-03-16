export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

charge=0
multiplicity=1

mkdir 0_geoOpt
mkdir 1_resp
mkdir 2_top

cd 0_geoOpt

obabel ../struct.pdb -O opt.com

sed -i 1d opt.com
sed -i "5s/.*/$charge  $multiplicity/" opt.com
sed -i '1s@^@#p opt=(NoLinear,MaxStep=5) guess=indo b3lyp/6-31G(d)\n@' opt.com
sed -i '1s@^@%mem=15GB\n@' opt.com
sed -i '1s@^@%nprocshared=10\n@' opt.com

g16 < opt.com > opt.log

obabel -ig09 opt.log -O resp.com

echo "g16.gesp" >> resp.com 
echo -en "\n" >> resp.com 
echo "g16.gesp" >> resp.com 
echo -en "\n" >> resp.com

sed -i 1d resp.com
sed -i '1s@^@# b3lyp/6-31G(d) pop=MK iop(6/33=2,6/41=10,6/42=17,6/50=1) SCF=tight\n@' resp.com
sed -i '1s@^@%chk=molecule.chk\n@' resp.com
sed -i '1s@^@%mem=15GB\n@' resp.com
sed -i '1s@^@%nprocshared=10\n@' resp.com

cp resp.com ../1_resp

cd ../1_resp

# obabel -ig09 opt.log -O resp.com

# echo "g16.gesp" >> resp.com 
# echo -en "\n" >> resp.com 
# echo "g16.gesp" >> resp.com 
# echo -en "\n" >> resp.com

# sed -i '6s/.*/-1  1/' resp.com  # replace this
# sed -i "5s/.*/$charge  $multiplicity/" resp.com 
# sed -i 1d resp.com
# sed -i '1s@^@# b3lyp/6-31G(d) pop=MK iop(6/33=2,6/41=10,6/42=17,6/50=1) SCF=tight\n@' resp.com
# sed -i '1s@^@%chk=molecule.chk\n@' resp.com
# sed -i '1s@^@%mem=15GB\n@' resp.com
# sed -i '1s@^@%nprocshared=10\n@' resp.com

g16 < resp.com > resp.log
cp g16.gesp ../2_top

cd ../2_top
antechamber -fi gesp -i g16.gesp -fo mol2 -o MOL.mol2 -at gaff2 -nc $charge
parmchk2 -i MOL.mol2 -f mol2 -o MOL.frcmod -a Y
acpype -i MOL.mol2 -c user -o gmx -n $charge

python neutralize.py MOL.acpype/MOL_GMX.itp $charge

cd ../..
