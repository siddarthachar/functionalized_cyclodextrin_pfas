import sys
import numpy as np

# Usage:
#   python get_lowest_n_carbons.py probe_shifted.pdb 48
# The output file backbone_C_CD.dat is then consumed by the gamma/alpha-CD
# helper scripts to define the restrained backbone atoms and PLUMED planes.
pdb_filename = sys.argv[1]
backbone_C_total = int(sys.argv[2])

def get_lowest_n_carbons(pdb_file, n):
    # Store the original ATOM record alongside coordinates so the selected
    # names can be written back out exactly as they appeared in the PDB.
    carbon_atoms = []
    
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and line[76:78].strip() == "C":  # Check if it's a Carbon atom
                atom_name = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                carbon_atoms.append((atom_name, x, y, z, line))
    
    # Sort carbons from the lower rim upwards. For gamma/alpha-CD this acts as a
    # simple geometric way of identifying the backbone-side atoms.
    carbon_atoms.sort(key=lambda atom: atom[3])
    
    # The current script keeps the first n+1 entries; preserve this behavior
    # because downstream files may have been prepared with the same convention.
    lowest_carbons = carbon_atoms[:n+1]
    
    return [atom[4] for atom in lowest_carbons]  # Return original PDB lines

lowest_c_atoms = get_lowest_n_carbons(pdb_filename, backbone_C_total)

C_numbers = []
# Extract the atom names only (for example C1, C2, ...), because the next
# scripts just need the cyclodextrin carbon labels rather than full PDB lines.
for atom_line in lowest_c_atoms:
    C_numbers += [atom_line.split()[2]]

output_filename="backbone_C_CD.dat"
with open(output_filename, "w") as output_file:
    for atom_index in C_numbers:
        output_file.write(f"{atom_index}\n")
