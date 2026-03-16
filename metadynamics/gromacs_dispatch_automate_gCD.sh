#!/bin/bash
# please change according to your systems - just use these scripts for reference. Do not run directly. 
# Reference setup script for gamma-cyclodextrin systems. It is intentionally
# verbose because the backbone selection and PLUMED plane generation are the
# parts readers usually need to adapt.
shopt -s expand_aliases
source ~/.bashrc
prepgmx
export PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/bin:$PATH" # path to gromacs
export LD_LIBRARY_PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/lib:$LD_LIBRARY_PATH"

#!/bin/bash

sim_main_path='.' # main path to simulations
ff_path='/project2/andrewferguson/sachar/test/ff-parameterize/' # contains gaff2 parameterization with RESP partial charges
systems=('') # list out all systems to iterate

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

    # Infer analyte type and guest parameterization directory from the name.
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

    # Remove helper files that should not be part of a clean setup.
    for file in "${files2remove[@]}"; do
        rm "$file"
    done

    # Require a completed guest parameterization before building the system.
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

    # Merge probe atom types into the shared system include file.
    python ~/bin/transfer_probe-ff2sys-ff.py probe-ff.itp sys-ffatoms.itp

    # Build the initial host-guest structure.
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

    # Solvate and neutralize before the MD stages begin.
    gmx solvate -cp merged.gro -cs spc216.gro -o sys-solvated.gro -p ../sys.top
    gmx grompp -f ../emin.mdp -c sys-solvated.gro -p ../sys.top -po sys-neutral.mdp -o sys-neutral.tpr -maxwarn 6 -r sys-solvated.gro
    # gmx genion -s sys-neutral.tpr -p ../sys.top -o sys-neutral.gro -pname NA -pq 1 -np 1
    printf "SOL\n" | gmx genion -s sys-neutral.tpr -p ../sys.top -o sys-neutral.gro -pname NA -pq 1 -np 1
    cp sys-neutral.gro ../
    cd .. || exit 1

    # Reindex the probe PDB to the final host-guest numbering and then derive
    # the gamma-CD backbone atoms used for both PLUMED planes and restraints.
    python ../increment_pbd_file.py struct/probe.pdb "$increment" probe_shifted.pdb

    # TODO: for gcd write a command here to get the backbone gcd C atoms
    python ../get_lowest_n_carbons.py probe_shifted.pdb 48 # 48 is for gCD. Change for others accordingly

    python ../get_plane_gCD_plumed.py probe_shifted.pdb plumed_template.dat "$increment" 
    # update the existing plumed.dat 
    update_plumed.sh "$analyte"

    python ../generate_posres_gcd_acd_C.py probe_shifted.pdb posre_C.itp # TODO: need to change 42 to something else based on the base atom
    
    # Append the host restraint include so POSRES activates it automatically.
    echo -e "\n#ifdef POSRES\n#include \"posre_C.itp\"\n#endif" >> probe-ff.itp
    
    cp ../prod-equil-afGPU.sh .
    cp ../prod-continue-afGPU.sh .
    sbatch prod-equil-afGPU.sh

    cd .. || exit 1

done
