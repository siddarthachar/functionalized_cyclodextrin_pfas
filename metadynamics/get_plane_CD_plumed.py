import MDAnalysis as mda
from sklearn.cluster import KMeans
import numpy as np
import sys
import warnings
from dscribe.descriptors import SOAP
from ase import Atoms
from MDAnalysis.topology.guessers import guess_atom_element
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# This script fills the "# fill here" placeholder in a PLUMED template with
# three center-of-mass definitions corresponding to the bottom, middle, and top
# carbon layers of cyclodextrin.
removeNum = lambda s: ''.join([i for i in s if not i.isdigit()])
# Map atom names that MDAnalysis does not always infer correctly.
def get_element(atom_name):
    element_map = {
        "CL": "Cl",  # Add more mappings as needed
        "BR": "Br"
    }
    return element_map.get(atom_name, guess_atom_element(atom_name))

probe_file = sys.argv[1]
plumed_file = sys.argv[2]
index_shift = int(sys.argv[3])
max_carbon = int(sys.argv[4])

# Load the probe-only PDB and guess elements so SOAP can distinguish chemically
# similar atoms by local environment instead of coordinates alone.
u = mda.Universe(probe_file)
positions = u.atoms.positions
guessed_elements = [get_element(removeNum(atom.name)) for atom in u.atoms]

# Limit the layer analysis to the cyclodextrin backbone carbons.
carbon_indices = [
    i for i, atom in enumerate(u.atoms)
    if guessed_elements[i] == 'C' and atom.name.startswith('C') and int(atom.name[1:]) < max_carbon # change this to in the list loaded
]
carbon_positions = positions[carbon_indices]

# Convert to an ASE object because DScribe expects that interface.
ase_atoms = Atoms(symbols=guessed_elements, positions=positions)

# SOAP descriptors encode each carbon's local environment. Combining them with
# Cartesian positions makes it easier to separate the three rings robustly.
soap_descriptor = SOAP(
    species=list(set(guessed_elements)),  
    r_cut=4.0,                   
    n_max=6,                      
    l_max=6,                      
    average='off',                
    sparse=False
)

soap_vectors = soap_descriptor.create(ase_atoms, centers=carbon_positions)
scaler = StandardScaler()
scaled_positions = scaler.fit_transform(carbon_positions)

# Cluster in the joint descriptor/geometry space.
combined_features = np.hstack((soap_vectors, scaled_positions))

# Three clusters correspond to the lower, middle, and upper carbon layers.
kmeans = KMeans(n_clusters=3, random_state=42)
clusters = kmeans.fit_predict(combined_features)

# Use the mean z coordinate of each cluster to map arbitrary KMeans labels onto
# physically meaningful layer names.
z_coords = carbon_positions[:, 2]
cluster_z_means = [np.mean(z_coords[clusters == i]) for i in range(3)]
sorted_clusters = np.argsort(cluster_z_means)

layer_mapping = {
    sorted_clusters[0]: "Bottom",
    sorted_clusters[1]: "Middle",
    sorted_clusters[2]: "Top"
}

# Shift to the atom numbering used by the host-guest system in PLUMED/GROMACS.
layer_indices = {"Top": [], "Middle": [], "Bottom": []}
for idx, cluster in zip(carbon_indices, clusters):
    layer = layer_mapping[cluster]
    layer_indices[layer].append(idx + index_shift)

# PLUMED expects explicit atom lists for the three reference planes.
str_top = ','.join(map(str, sorted(layer_indices["Bottom"])))
str_mid = ','.join(map(str, sorted(layer_indices["Middle"])))
str_bot = ','.join(map(str, sorted(layer_indices["Top"])))

replacement_text = f"""planeBottom: COM ATOMS={str_bot}
planeMiddle: COM ATOMS={str_mid}
planeTop: COM ATOMS={str_top}
"""

# Replace the placeholder block in the template rather than rebuilding the
# entire file, so all other metadynamics settings remain untouched.
with open(plumed_file, "r") as file:
    lines = file.readlines()

with open(plumed_file, "w") as file:
    for line in lines:
        if "# fill here" in line:
            file.write(replacement_text)
        else:
            file.write(line)

print("Updated plumed.dat successfully!")
