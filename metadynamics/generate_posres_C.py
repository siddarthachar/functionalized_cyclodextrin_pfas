import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
import sys
import warnings
warnings.filterwarnings("ignore")

# Usage:
#   python generate_posres_C.py probe_shifted.pdb posre_C.itp 42
# The script scans the shifted cyclodextrin PDB and writes a GROMACS
# position-restraint include file for the carbon atoms that define the
# comparatively rigid ring scaffold.
probe_file = sys.argv[1]
posre_file = sys.argv[2]
max_C_index = int(sys.argv[3]) # 42 for bCD

# Load the probe-only structure; the atom numbering in this file should already
# match the numbering used in the combined simulation system.
u = mda.Universe(probe_file)
positions = u.atoms.positions
guessed_elements = [guess_atom_element(atom.name) for atom in u.atoms]
symbols = guessed_elements

# Keep only ring carbons whose names follow the C<number> convention and fall
# below the user-supplied cutoff. This is how the bCD backbone is selected.
carbon_indices = [
    i + 1 for i, atom in enumerate(u.atoms)
    if symbols[i] == 'C' and atom.name.startswith('C') and int(atom.name[1:]) < max_C_index
]

# Write a standard GROMACS position-restraint table.
with open(posre_file, 'w') as f:
    f.write(f"; generated for C atoms less than {max_C_index} (part of beta-CD). Change {max_C_index} to another number\n\n")
    f.write("[ position_restraints ]\n")
    f.write(";  i funct       fcx        fcy        fcz\n")
    for index in carbon_indices:
        f.write(f"  {index:4d}    1       1000       1000       1000\n")

# Echo the selected atoms so the caller can sanity-check the restraint set.
print(f"Carbon indices used for position restraints: {carbon_indices}")
