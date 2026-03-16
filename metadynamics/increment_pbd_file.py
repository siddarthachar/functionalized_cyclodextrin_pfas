import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
import sys

# Usage:
#   python increment_pbd_file.py struct/probe.pdb 30 probe_shifted.pdb
# The shifted PDB is used only as a numbering/reference file for PLUMED and the
# restraint generators; the coordinates are unchanged.
input_file = sys.argv[1]
increment = int(sys.argv[2]) # use 30 if incremenet is 29
out_file = sys.argv[3]

def increment_atom_numbers(input_pdb, output_pdb, increment):
    """
    Reads a PDB file, increments atom numbers by a specified value, and writes a new PDB file.
    
    Parameters:
    - input_pdb (str): Path to the input PDB file.
    - output_pdb (str): Path for the output PDB file.
    - increment (int): Value to add to the atom numbers.
    """
    # MDAnalysis gives convenient access to atom metadata while preserving the
    # coordinates written by the docking/structure-preparation step.
    u = mda.Universe(input_pdb)

    with open(output_pdb, 'w') as f:
        # Copy any header records verbatim until atom records begin.
        with open(input_pdb, 'r') as original:
            for line in original:
                if line.startswith(('ATOM', 'HETATM', 'TER', 'END')):
                    break
                f.write(line)

        # Rewrite each atom record with the requested atom-id offset. The
        # residue information and coordinates are left as-is.
        for atom in u.atoms:
            new_atom_id = atom.id + increment
            line = f"ATOM  {new_atom_id:5d} {atom.name:<4} {atom.resname:<3} {atom.resid:4d}    " \
                   f"{atom.position[0]:8.3f}{atom.position[1]:8.3f}{atom.position[2]:8.3f}" \
                   f"{atom.occupancy:6.2f}{atom.bfactor:6.2f}          {guess_atom_element(atom.name):>2}\n"
            f.write(line)

        f.write("END\n")

# The historical workflow passes the analyte atom count and this script applies
# increment-1 internally so the first probe atom starts immediately after the
# analyte numbering.
increment_atom_numbers(input_file,
                       out_file, increment=increment-1)
