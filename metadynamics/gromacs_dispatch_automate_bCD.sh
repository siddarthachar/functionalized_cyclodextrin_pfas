#!/bin/bash
# please change according to your systems - just use these scripts for reference. Do not run directly. 
# Reference setup script for beta-cyclodextrin systems. The intent is to show
# the order of operations needed to assemble and launch a parallel-bias MetaD
# calculation from force-field outputs plus common templates.
shopt -s expand_aliases
source ~/.bashrc
prepgmx
export PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/bin:$PATH" # path to gromacs
export LD_LIBRARY_PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/lib:$LD_LIBRARY_PATH"

#!/bin/bash

sim_main_path='.' # main path to simulations
ff_path='/project2/andrewferguson/sachar/test/ff-parameterize/' # contains gaff2 parameterization with RESP partial charges
systems=('bcd-pfos') # list out all systems to iterate

cd "$sim_main_path" || exit 1

for system in "${systems[@]}"; do
    # Check if directory exists; if not, create it
    if [[ ! -d "$system" ]]; then
        echo "Directory $system does not exist. Creating it now."
        mkdir -p "$system"
    fi

    cd "$system" || exit 1
    echo "$system | Starting ------------------------"

    # Infer analyte type and atom-number offset from the system name.
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

    # Remove helper files that should not be copied into a fresh system folder.
    for file in "${files2remove[@]}"; do
        rm "$file"
    done

    # The guest topology must already exist from the RESP/GAFF workflow.
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

    # Insert the guest atom-type definitions into the combined topology include.
    python ~/bin/transfer_probe-ff2sys-ff.py probe-ff.itp sys-ffatoms.itp

    # Build the initial docked/merged coordinates in struct/.
    mkdir -p struct
    cd struct || exit 1
    cp "$ff_path/$probe/2_top/MOL.acpype/MOL_GMX.gro" probe.gro

    replace_MOL_gro.sh probe.gro # Replace all occurrences of MOL with PRO
    mv probe_pro.gro probe.gro # Rename the pro file to probe.gro
    cp "../$analyte_pdb" .

    # Prepare files for docking
    gmx editconf -f probe.gro -o probe.pdb -box 7 7 7 -center 3.5 3.5 3.5
    gmx editconf -f $analyte_pdb -o $analyte_pdb -translate 0 0 1.0
    cat "$analyte_pdb" probe.pdb > merged.pdb
    fix_pdb.sh merged.pdb
    gmx editconf -f merged.pdb -o merged.gro
    gmx genconf -f merged.gro -o merged.gro -renumber

    # Solvate the box and add ions before minimization.
    gmx solvate -cp merged.gro -cs spc216.gro -o sys-solvated.gro -p ../sys.top
    gmx grompp -f ../emin.mdp -c sys-solvated.gro -p ../sys.top -po sys-neutral.mdp -o sys-neutral.tpr -maxwarn 6 -r sys-solvated.gro
    # gmx genion -s sys-neutral.tpr -p ../sys.top -o sys-neutral.gro -pname NA -pq 1 -np 1
    # printf "SOL\n" | gmx genion -s sys-neutral.tpr -p ../sys.top -o sys-neutral.gro -pname NA -pq 1 -np 1
    printf "SOL\n" | gmx genion -s sys-neutral.tpr -p ../sys.top -o sys-neutral.gro -pname NA -pq 1 -nname CL -nq -1 -neutral # more general
    cp sys-neutral.gro ../
    cd .. || exit 1

    # Reindex the host atoms to the final host-guest numbering, fill the PLUMED
    # plane definitions, and write the host restraint include file.
    python ../increment_pbd_file.py struct/probe.pdb "$increment" probe_shifted.pdb
    # for gcd write a command here to get the backbone gcd C atoms
    python ../get_plane_CD_plumed.py probe_shifted.pdb plumed_template.dat "$increment" 42 # TODO: need to change 42 to something else based on the base atom
    # update the existing plumed.dat 
    update_plumed.sh "$analyte"

    python ../generate_posres_C.py probe_shifted.pdb posre_C.itp 42 # TODO: need to change 42 to something else based on the base atom
    
    # Include the rigid-host restraints when the topology is built with POSRES.
    echo -e "\n#ifdef POSRES\n#include \"posre_C.itp\"\n#endif" >> probe-ff.itp
    # if gpu
    cp ../prod-equil.sh .
    cp ../prod-continue.sh .
    # sbatch prod-equil.sh
    
    # if afGPU
    # cp ../prod-equil-afGPU.sh .
    # cp ../prod-continue-afGPU.sh .
    # sbatch prod-equil-afGPU.sh

    cd .. || exit 1
done
