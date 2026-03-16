import os
import numpy as np
import matplotlib.pyplot as plt
import nglview
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy.optimize import minimize

def get_nrotor(file):
    mol = Chem.MolFromPDBFile(file)
    return Descriptors.NumRotatableBonds(mol), Descriptors.MolWt(mol)

nr_pfos, m_pfos = get_nrotor('/project2/andrewferguson/sachar/test/dE_analysis/gas-phase/pfos/no-sol.pdb')
nr_sds, m_sds = get_nrotor('/project2/andrewferguson/sachar/test/dE_analysis/gas-phase/sds/no-sol.pdb')
nr_tcaa, m_tcaa = get_nrotor('/project2/andrewferguson/sachar/test/dE_analysis/gas-phase/tcaa/no-sol.pdb')
nr_bcd, m_bcd = get_nrotor('/project2/andrewferguson/sachar/test/dE_analysis/gas-phase/bcd/no-sol.pdb')
nr_gcd, m_gcd = get_nrotor('/project2/andrewferguson/sachar/test/dE_analysis/gas-phase/gcd/no-sol.pdb')

## comparing just gas phase energies 
main_path='/project2/andrewferguson/sachar/test/dE_analysis/gas-phase'
bcd = float(os.popen(f'grep Energy {main_path}/bcd/sys-emin.log').read().split('\n')[-2].split()[-1])
gcd = float(os.popen(f'grep Energy {main_path}/gcd/sys-emin.log').read().split('\n')[-2].split()[-1])
pfos = float(os.popen(f'grep Energy {main_path}/pfos/sys-emin.log').read().split('\n')[-2].split()[-1])
sds = float(os.popen(f'grep Energy {main_path}/sds/sys-emin.log').read().split('\n')[-2].split()[-1])
tcaa=float(os.popen(f'grep Energy {main_path}/tcaa/sys-emin.log').read().split('\n')[-2].split()[-1])

pfos_bcd = float(os.popen(f'grep Energy {main_path}/pfos-f-bcd/sys-emin.log').read().split('\n')[-2].split()[-1])
sds_bcd = float(os.popen(f'grep Energy {main_path}/sds-f-bcd/sys-emin.log').read().split('\n')[-2].split()[-1])
tcaa_bcd = float(os.popen(f'grep Energy {main_path}/tcaa-f-bcd/sys-emin.log').read().split('\n')[-2].split()[-1])
tcaa_gcd = float(os.popen(f'grep Energy {main_path}/tcaa-f-gcd/sys-emin.log').read().split('\n')[-2].split()[-1])

dH_pfos_bcd = pfos_bcd - bcd - pfos
dH_sds_bcd = sds_bcd - bcd - sds
dH_tcaa_bcd = tcaa_bcd - bcd - tcaa
dH_tcaa_gcd = tcaa_gcd - gcd - tcaa

# Sample data
dG_temp = np.array([-36.626, -27.967, -7.2495, -9.87590])
dE = np.array([dH_pfos_bcd, dH_sds_bcd, dH_tcaa_bcd, dH_tcaa_gcd])

# Properties for each ligand-receptor pair
m_values = [m_pfos, m_sds, m_tcaa, m_tcaa]
nr_values = [nr_pfos, nr_sds, nr_tcaa, nr_tcaa]
m_receptors = [m_bcd, m_bcd, m_bcd, m_gcd]
nr_receptors = [nr_bcd, nr_bcd, nr_bcd, nr_gcd]

# Define function to calculate dS
def calculate_dS(w1, w2):
    return np.array([
        w1 * np.log(m) + w1 * np.log(m_rec) + w2 * nr + w2 * nr_rec
        for m, nr, m_rec, nr_rec in zip(m_values, nr_values, m_receptors, nr_receptors)
    ])

# Define the objective function to minimize
def objective(weights):
    w1, w2 = weights
    dS = calculate_dS(w1, w2)
    dG_est = dE + dS
    # Calculate the sum of squared differences
    return np.sqrt(np.sum((dG_est - dG_temp) ** 2))

# Perform minimization
result = minimize(objective, x0=(1, 1), bounds=[(0, None), (0, None)])
w1_opt, w2_opt = result.x

# Calculate dG_est with optimized weights
dS_opt = calculate_dS(w1_opt, w2_opt)
dG_est_opt = dE + dS_opt
labels = ['PFOS-$\\beta CD$', 'SDS-$\\beta CD$', 'TCAA-$\\beta CD$', 'TCAA-$\\gamma CD$']
plt.bar(labels, dS_opt)
plt.ylabel('dS approx. (arb)')
plt.savefig('dS_compare.png',dpi=400)
plt.close()
dG_gas = [-90.67, -101.37, -0.42, -14.61]
# Plot results

plt.plot(dG_est_opt, dG_temp, 'bo', label='approx')
plt.plot(dG_gas, dG_temp, 'rx', label='dG gas-phase metaD')

for i, label in enumerate(labels):
    plt.text(dG_est_opt[i] + 0.5, dG_temp[i] + 2, label, fontsize=8, color='b')  # Adjust offset for clarity
    plt.text(dG_gas[i] + 0.5, dG_temp[i] + 2, label, fontsize=8, color='r')  # Adjust offset for clarity

# Linear fit line
m, b = np.polyfit(dG_est_opt, dG_temp, 1)
dG_range = np.linspace(min(dG_est_opt) - 40, max(dG_est_opt) + 40)
dG_proxy = m * dG_range + b
plt.plot(dG_range, dG_proxy, 'c--')
plt.legend(frameon=False)

plt.title('Gas phase')
plt.ylabel('$\\Delta G$')
plt.xlabel('$\\Delta G_{proxy}$')
plt.savefig('compare_dG_gas.png',dpi=400)

print("Optimized weights:")
print("w1 =", w1_opt)
print("w2 =", w2_opt)