import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
import sys
import warnings
import os
warnings.filterwarnings("ignore")

# Usage:
#   python generate_posres_gcd_acd_C.py probe_shifted.pdb posre_C.itp
# This variant is used for alpha/gamma-CD systems, where the restrained carbon
# set is prepared beforehand and stored in backbone_C_CD.dat.
probe_file = sys.argv[1]
posre_file = sys.argv[2]
# max_C_index = int(sys.argv[3]) # 42 for bCD
C_numbers = os.popen('cat backbone_C_CD.dat').read().split() # file that contains backbone C atoms as C1 .. Cn
C_just_numbers = [int(x[1:]) for x in C_numbers] # stripping in the C to just have the int


# Load the probe-only structure whose numbering already matches the combined
# host-guest system after the analyte atom offset has been applied.
u = mda.Universe(probe_file)
positions = u.atoms.positions
guessed_elements = [guess_atom_element(atom.name) for atom in u.atoms]
symbols = guessed_elements

# Select only the carbon atoms explicitly listed in backbone_C_CD.dat.
carbon_indices = [
    i + 1 for i, atom in enumerate(u.atoms)
    if symbols[i] == 'C' and atom.name.startswith('C') and int(atom.name[1:]) in C_just_numbers
]

# Write a standard GROMACS position-restraint table.
with open(posre_file, 'w') as f:
    f.write(f"; generated for C atoms within backbone (part of beta-CD).\n\n")
    f.write("[ position_restraints ]\n")
    f.write(";  i funct       fcx        fcy        fcz\n")
    for index in carbon_indices:
        f.write(f"  {index:4d}    1       1000       1000       1000\n")

# Echo the selected atoms so the caller can inspect the backbone definition.
print(f"Carbon indices used for position restraints: {carbon_indices}")
