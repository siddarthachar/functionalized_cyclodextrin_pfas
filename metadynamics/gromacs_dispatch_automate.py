import os
import subprocess

# Legacy Python dispatcher that assembles one simulation directory at a time:
# 1. copy common MD input files,
# 2. pull in the guest force field from the parameterization workflow,
# 3. build/solvate/neutralize the host-guest complex,
# 4. generate PLUMED and restraint helper files.
# The shell versions in this repository are the more current entry points, but
# this script is still useful as a readable record of the workflow.
sim_main_path = '/project2/andrewferguson/sachar/test/rigid_CD_gaff/'
ff_path = '/project2/andrewferguson/sachar/test/ff-parameterize/'
systems = ['CCl2CCl3-prim-bcd-pfos','Cl-prim-bcd-pfos','CCl2CCl2CCl3-prim-bcd-pfos','CCl3-prim-bcd-pfos','Cl-both-bcd-pfos','CH3-prim-bcd-pfos','CH2CH3-prim-bcd-pfos','CH2CH2CH3-prim-bcd-pfos']
files2remove = ['bcd_30.pdb','segment_fes.sbatch','run_segment_fes.sh','log_analysis.py','prod-continue-gas.sh',
               'minima_conf_30234ps.pdb','segment_fes.sh','compare_dG_gas_ref.py','struct1.pdb','analyte-posres.itp','posre_C.itp']
os.chdir(sim_main_path)
for i, system in enumerate(systems):
    os.chdir(system)
    print(system,' | Starting ------------------------')

    # Copy the common simulation templates for the analyte family and define the
    # atom-number offset needed once the guest atoms are prepended to the host.
    if 'pfos' in system:
        os.system(f'copy_gmx_from.py {sim_main_path}/common.files/pfos/ .')
        analyte_pdb = 'pfos.pdb'
        increment = 30
    elif 'sds' in system:
        os.system(f'copy_gmx_from.py {sim_main_path}/common.files/sds/ .')
        analyte_pdb = 'sds.pdb'
        increment = 43

    # Remove templates and analysis helpers that are not needed for the active
    # system directory.
    [os.system(f'rm {x}') for x in files2remove]

    # Pull the guest force field generated in the separate RESP/GAFF workflow.
    if 'MOL.mol2' not in os.listdir(f'{ff_path}/{system}/2_top'):
        print(f'FF not generated for {system} - please go back and check')
        print('Skipping now ..............')
        os.chdir('../')
        continue
    else:
        print(f'FF found for {system}')
        os.system(f'cp {ff_path}/{system}/2_top/MOL.acpype/MOL_GMX.itp probe-ff.itp')
        os.system('replace_MOL_itp.sh probe-ff.itp') # replaces all occurrences of MOL with PRO
        
    # Copy the probe atom-type declarations into the system include file so the
    # combined topology knows about the guest-specific parameters.
    os.system('transfer_probe-ff2sys-ff.py probe-ff.itp sys-ffatoms.itp')

    # Build the initial host-guest coordinates inside struct/.
    os.makedirs('struct', exist_ok=True)
    os.chdir('struct')
    os.system(f'cp {ff_path}/{system}/2_top/MOL.acpype/MOL_GMX.gro probe.gro')

    os.system('replace_MOL_gro.sh probe.gro') # replaces all occurrences of MOL with PRO
    os.system('mv probe_pro.gro probe.gro') # rename the pro file to probe.gro
    os.system(f'cp ../{analyte_pdb} struct')
    
    os.system('gmx editconf -f probe.gro -o probe.pdb -box 7 7 7 -center 3.5 3.5 3.5') # uncomment at the end ##########
    os.system('gmx editconf -f pfos.pdb -o pfos.pdb -translate 0 0 0.5') # uncomment at the end ##########
    os.system(f'cat {analyte_pdb} probe.pdb > merged.pdb')
    os.system('fix_pdb.sh merged.pdb')
    os.system('gmx editconf -f merged.pdb -o merged.gro')
    os.system('gmx genconf -f merged.gro -o merged.gro -renumber')

    # Solvate the merged system and add ions so GROMACS can start from a
    # neutral periodic box.
    os.system('gmx solvate -cp merged.gro -cs spc216.gro -o sys-solvated.gro -p ../sys.top')
    os.system('gmx grompp -f ../emin.mdp -c sys-solvated.gro -p ../sys.top -po sys-neutral.mdp -o sys-neutral.tpr -maxwarn 6 -r sys-solvated.gro')
    os.system('gmx genion -s sys-neutral.tpr -p ../sys.top -o sys-neutral.gro -pname NA -pq 1 -np 1') # Can you automatically make it select SOL? ####### TODO
    os.system('cp sys-neutral.gro ../')
    os.chdir('../')

    # Reindex the probe-only PDB so its atom numbers match the combined system,
    # then generate the PLUMED plane definitions and host restraints from that
    # shifted numbering.
    os.system(f'python ../increment_pbd_file.py struct/probe.pdb {increment} probe_shifted.pdb')
    os.system(f'python ../get_plane_CD_plumed.py probe_shifted.pdb plumed.dat 30 42') # replace the 42 with the number that you want
    os.system(f'python ../generate_posres_C.py probe_shifted.pdb {increment} posre_C.itp')


    os.chdir('../')
