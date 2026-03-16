#!/bin/bash
# Reference setup script for alpha-cyclodextrin host-guest systems.
# It assembles each simulation directory from force-field outputs plus the
# shared MD templates, then submits the production/equilibration job.
shopt -s expand_aliases
source ~/.bashrc
prepgmx
export PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/bin:$PATH"
export LD_LIBRARY_PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/lib:$LD_LIBRARY_PATH"

#!/bin/bash

sim_main_path='/project2/andrewferguson/sachar/test/rigid_CD_gaff/'
ff_path='/project2/andrewferguson/sachar/test/ff-parameterize/'
systems=('F-1-acd-pfos' 'F-1-acd-sds' 'NH2-4-acd-pfos' 'NH2-4-acd-sds')
files2remove=('bcd_30.pdb' 'bcd_43.pdb' 'segment_fes.sbatch' 'run_segment_fes.sh' 'log_analysis.py' 'prod-continue-gas.sh' 'minima_conf_30234ps.pdb' 'segment_fes.sh' 'compare_dG_gas_ref.py' 'struct1.pdb' 'analyte-posres.itp' 'posre_C.itp')

cd "$sim_main_path" || exit 1

for system in "${systems[@]}"; do
    # Check if directory exists; if not, create it
    if [[ ! -d "$system" ]]; then
        echo "Directory $system does not exist. Creating it now."
        mkdir -p "$system"
    fi

    cd "$system" || exit 1
    echo "$system | Starting ------------------------"

    # The system name is used to infer both the analyte family and the probe
    # directory that contains the GAFF/RESP parameterization output.
    if [[ "$system" == *"-pfos" ]]; then
        analyte="pfos"
        probe="${system%-pfos}"  # Remove '-pfos' from system name
        cp "$sim_main_path/common.files/pfos/"* .
        analyte_pdb='pfos.pdb'
        increment=30
    elif [[ "$system" == *"-sds" ]]; then
        analyte="sds"
        probe="${system%-sds}"  # Remove '-sds' from system name
        cp "$sim_main_path/common.files/sds/"* .
        analyte_pdb='sds.pdb'
        increment=43
    else
        echo "Unknown analyte in system name: $system"
        cd .. || exit 1
        continue
    fi

    # Clean out templates that should not be carried into a fresh setup.
    for file in "${files2remove[@]}"; do
        rm "$file"
    done

    # Confirm that the guest parameterization workflow already finished.
    if [[ ! -f "$ff_path/$probe/2_top/MOL.mol2" ]]; then
        echo "FF not generated for $system - please go back and check"
        echo "Skipping now .............."
        cd .. || exit 1
        continue
    else
        echo "FF found for $system"
        cp "$ff_path/$probe/2_top/MOL.acpype/MOL_GMX.itp" probe-ff.itp
        replace_MOL_itp.sh probe-ff.itp # Replace all occurrences of MOL with PRO
    fi

    # Merge the probe atom-type declarations into the combined system include.
    python ~/bin/transfer_probe-ff2sys-ff.py probe-ff.itp sys-ffatoms.itp

    # Build the initial host-guest coordinate set.
    mkdir -p struct
    cd struct || exit 1
    cp "$ff_path/$probe/2_top/MOL.acpype/MOL_GMX.gro" probe.gro

    replace_MOL_gro.sh probe.gro # Replace all occurrences of MOL with PRO
    mv probe_pro.gro probe.gro # Rename the pro file to probe.gro
    cp "../$analyte_pdb" .

    # Prepare files for docking
    gmx editconf -f probe.gro -o probe.pdb -box 7 7 7 -center 3.5 3.5 3.5
    gmx editconf -f $analyte_pdb -o $analyte_pdb -translate 0 0 0.5
    cat "$analyte_pdb" probe.pdb > merged.pdb
    fix_pdb.sh merged.pdb
    gmx editconf -f merged.pdb -o merged.gro
    gmx genconf -f merged.gro -o merged.gro -renumber

    # Solvate the merged structure and neutralize the box.
    gmx solvate -cp merged.gro -cs spc216.gro -o sys-solvated.gro -p ../sys.top
    gmx grompp -f ../emin.mdp -c sys-solvated.gro -p ../sys.top -po sys-neutral.mdp -o sys-neutral.tpr -maxwarn 6 -r sys-solvated.gro
    # gmx genion -s sys-neutral.tpr -p ../sys.top -o sys-neutral.gro -pname NA -pq 1 -np 1
    printf "SOL\n" | gmx genion -s sys-neutral.tpr -p ../sys.top -o sys-neutral.gro -pname NA -pq 1 -np 1
    cp sys-neutral.gro ../
    cd .. || exit 1

    # Reindex the probe-only PDB to the final host-guest numbering, derive the
    # host backbone atoms, and use that numbering to generate PLUMED plus the
    # host position restraints.
    python ../increment_pbd_file.py struct/probe.pdb "$increment" probe_shifted.pdb

    # TODO: for gcd write a command here to get the backbone gcd C atoms
    python ../get_lowest_n_carbons.py probe_shifted.pdb 36 # 48 is for gCD, 36 for aCD Change for others accordingly

    python ../get_plane_gCD_plumed.py probe_shifted.pdb plumed_template.dat "$increment" 
    # update the existing plumed.dat 
    update_plumed.sh "$analyte"

    python ../generate_posres_gcd_acd_C.py probe_shifted.pdb posre_C.itp # TODO: need to change 42 to something else based on the base atom
    # the generate_posres_gcd_C.py code works for both acd and gcd
    
    # Make the new host restraint file available when POSRES is enabled.
    echo -e "\n#ifdef POSRES\n#include \"posre_C.itp\"\n#endif" >> probe-ff.itp
    
    # cp ../prod-equil.sh .
    # cp ../prod-continue.sh .
    # sbatch prod-equil.sh

    # Swap in the afGPU job scripts used on the target cluster and submit.
    cp ../prod-equil-afGPU.sh prod-equil.sh
    cp ../prod-continue-afGPU.sh prod-continue.sh
    sbatch prod-equil.sh

    cd .. || exit 1

done
