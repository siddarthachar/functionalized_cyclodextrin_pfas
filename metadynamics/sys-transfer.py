import os 

# Small helper to copy a prepared example system into another directory. This
# is mainly a convenience script for duplicating a working input set.
syst = ['../00464-pfos/']

files = ['prod-equil.sh','emin.mdp','sys-neutral.gro','sys.top','equil.mdp','equil-npt.mdp','prod.mdp','posre_C.itp','sys-ffatoms.itp','pfos-anion-ff.itp', 'sds-ff.itp', 'bcd-ff.itp','tip3p.itp','ions.itp','analyte-posres.itp','reweight.sh', 'plumed.dat','*.pdb','probe-ff.itp']

for f in files:
    os.system(f'cp {syst[0]}/{f} 00464-pfos/')
