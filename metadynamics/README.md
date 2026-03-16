# Metadynamics

This directory contains the setup templates, helper scripts, and example systems used to build GROMACS + PLUMED metadynamics simulations for cyclodextrin host-guest complexes.

The files are organized around three related tasks:

1. parameterize a guest molecule
2. assemble a solvated host-guest simulation directory
3. run and reweight the metadynamics trajectory

## Directory Layout

### Top-level helper scripts

- `gromacs_dispatch_automate.py`
  Legacy Python dispatcher that documents the older end-to-end setup workflow for beta-CD systems.
- `gromacs_dispatch_automate_aCD.sh`
  Shell workflow for alpha-CD systems. Copies templates, imports guest force fields, builds the solvated box, generates PLUMED and restraint files, and submits the production job.
- `gromacs_dispatch_automate_bCD.sh`
  Beta-CD setup script in the same style.
- `gromacs_dispatch_automate_gCD.sh`
  Gamma-CD setup script in the same style.
- `increment_pbd_file.py`
  Renumbers the probe-only PDB so its atom IDs line up with the combined host-guest system. Used before generating PLUMED atom lists and position restraints.
- `get_plane_CD_plumed.py`
  Builds three cyclodextrin plane definitions for beta-CD-like systems by clustering backbone carbons and replacing the `# fill here` placeholder in a PLUMED template.
- `get_plane_gCD_plumed.py`
  Same idea for alpha/gamma-CD systems, but it uses a precomputed backbone carbon list from `backbone_C_CD.dat`.
- `generate_posres_C.py`
  Generates a GROMACS `posre_C.itp` restraint file for beta-CD backbone carbons using a simple atom-number cutoff.
- `generate_posres_gcd_acd_C.py`
  Generates `posre_C.itp` for alpha/gamma-CD systems using `backbone_C_CD.dat`.
- `get_lowest_n_carbons.py`
  Selects the lowest-lying cyclodextrin carbons from a shifted PDB and writes their atom names to `backbone_C_CD.dat`.
- `sys-transfer.py`
  Convenience script for copying a known-good example system into another directory.

### Common templates

- `common.files/pfos/`
  Shared templates for PFOS-containing systems.
- `common.files/sds/`
  Shared templates for SDS-containing systems.

These folders contain the reusable MD inputs and analysis helpers that the dispatch scripts copy into newly created system directories:

- `emin.mdp`, `equil.mdp`, `equil-npt.mdp`, `prod.mdp`
  GROMACS parameter files for minimization, equilibration, NPT equilibration, and production metadynamics.
- `sys.top`, `sys-ffatoms.itp`, `tip3p.itp`, `ions.itp`
  Topology includes for the solvated system.
- `pfos.pdb` or `sds.pdb`
  Guest/analyte coordinates.
- `plumed_template.dat`
  PLUMED template with a placeholder that is filled by the helper scripts.
- `prod-equil.sh`, `prod-continue.sh`
  Job scripts used to start or continue production.
- `reweight.sh`, `run_reweight.sh`, `reweight.sbatch`
  Reweighting workflow after the biased simulation finishes.
- `segment_fes.sh`, `run_segment_fes.sh`, `segment_fes.sbatch`
  Free-energy-surface segmentation utilities.
- `dG_interval.py`, `dG-plot.py`, `log_analysis.py`, `compare_dG_gas_ref.py`
  Post-processing and analysis helpers.
- `analyte-posres.itp`, `probe-posre.itp` or `posre_analyte.itp`
  Restraints applied to selected molecule subsets during setup or equilibration.

### Example system directories

- `00464-pfos/`
  Example PFOS system with a generated probe force field (`probe-ff.itp`) and a ready-to-run `plumed.dat`.
- `bcd-pfos/`
  Example beta-CD + PFOS system using a host force field file `bcd-ff.itp`.
- `bcd-sds/`
  Example beta-CD + SDS system using `sds-ff.itp`.

Each example system directory contains the files required to run a prepared simulation:

- `sys-neutral.gro`
  Neutral solvated coordinates ready for GROMACS preprocessing.
- `sys.top`
  Main topology file.
- `sys-ffatoms.itp`
  Atom-type include file that must contain both host and guest atom types.
- `plumed.dat`
  PLUMED control file with the host plane definitions already inserted.
- `posre_C.itp`
  Position restraints for the selected cyclodextrin backbone carbons.
- `prod-equil.sh`
  Batch script used to start the job on the target cluster.
- `reweight.sh`
  Post-processing script for reweighting.

### Guest force-field generation

- `ff-parameterize/`
  RESP/GAFF-style guest parameterization workflow.

Important files there:

- `struct.pdb`
  Starting guest structure.
- `0_geoOpt/opt.com`, `0_geoOpt/resp.com`
  Quantum-chemistry input files used during geometry optimization and RESP preparation.
- `2_top/MOL.mol2`, `2_top/MOL.frcmod`
  Generated guest topology inputs.
- `2_top/MOL.acpype/`
  ACPYPE-converted GROMACS topology and coordinate files such as `MOL_GMX.itp` and `MOL_GMX.gro`.
- `gen_top.sh`
  Topology-generation driver script.
- `sbatch.sbatch`, `sbatch_skip0.sbatch`
  Cluster submission scripts.

## Typical Workflow

### 1. Parameterize the guest

Populate or update `ff-parameterize/struct.pdb`, then run the parameterization workflow in `ff-parameterize/` so that `2_top/MOL.acpype/MOL_GMX.itp` and `2_top/MOL.acpype/MOL_GMX.gro` exist.

Without those files, the dispatch scripts will skip the system.

### 2. Assemble the simulation directory

Pick the dispatcher that matches the host family:

- `gromacs_dispatch_automate_aCD.sh`
- `gromacs_dispatch_automate_bCD.sh`
- `gromacs_dispatch_automate_gCD.sh`

These scripts do the following:

1. copy the shared template files from `common.files/pfos/` or `common.files/sds/`
2. copy the guest force field from `ff-parameterize/.../MOL.acpype/`
3. merge host and analyte structures
4. solvate and neutralize the system with GROMACS
5. renumber the probe-only structure with `increment_pbd_file.py`
6. generate `plumed.dat` and `posre_C.itp`
7. submit `prod-equil.sh`

### 3. Run or continue production

The prepared system directory normally contains:

- `prod-equil.sh` to start the first production segment
- `prod-continue.sh` in the template sets to continue an existing run

### 4. Reweight and analyze

After the simulation completes, use the reweighting and FES tools copied from `common.files/...`:

- `reweight.sh`
- `run_reweight.sh`
- `segment_fes.sh`
- `dG_interval.py`
- `dG-plot.py`
- `log_analysis.py`

## Direct Use Of The Helper Scripts

Examples from the repository workflow:

```bash
python metadynamics/increment_pbd_file.py struct/probe.pdb 30 probe_shifted.pdb
python metadynamics/get_lowest_n_carbons.py probe_shifted.pdb 36
python metadynamics/get_plane_gCD_plumed.py probe_shifted.pdb plumed_template.dat 30
python metadynamics/generate_posres_gcd_acd_C.py probe_shifted.pdb posre_C.itp
```

For beta-CD-like systems:

```bash
python metadynamics/get_plane_CD_plumed.py probe_shifted.pdb plumed.dat 30 42
python metadynamics/generate_posres_C.py probe_shifted.pdb posre_C.itp 42
```

## Environment Requirements

This directory assumes access to:

- `gmx` from GROMACS
- `plumed`
- Python packages used by the helper scripts:
  `MDAnalysis`, `numpy`, `scikit-learn`, `dscribe`, and `ase`
- cluster tooling such as `sbatch`
- a local shell environment that provides helper commands referenced in the scripts, for example `replace_MOL_itp.sh`, `replace_MOL_gro.sh`, `fix_pdb.sh`, `update_plumed.sh`, and `transfer_probe-ff2sys-ff.py`

Several scripts contain hard-coded project paths from the original cluster environment. If you run them elsewhere, update those path variables first.

## Notes

- The top-level dispatchers are designed for a specific HPC environment and are not immediately portable without path cleanup.
- The example system folders are the quickest way to see a fully assembled input set.
- The smallest reliable entry points for understanding the workflow are `gromacs_dispatch_automate_aCD.sh`, `increment_pbd_file.py`, `get_plane_*_plumed.py`, and `generate_posres*.py`.
