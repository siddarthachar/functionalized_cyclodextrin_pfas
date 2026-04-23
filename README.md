# Selective PFAS Detection with Functionalized Cyclodextrin Probes Designed via Bayesian Optimization

[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](#citation)

![Table of contents graphic](toc_scale.jpg)

This repository contains the computational assets associated with a study on cyclodextrin-based molecular recognition of PFAS, with an emphasis on improving PFOS selectivity over structurally similar surfactants such as SDS.

## Citation

**Title:** Selective PFAS Detection with Functionalized Cyclodextrin Probes Designed via Bayesian Optimization

**DOI:** Pending. Replace the badge and this line with the final DOI link after publication, for example:

```md
[![DOI](https://img.shields.io/badge/DOI-10.XXXX%2FXXXXX-blue)](https://doi.org/10.XXXX/XXXXX)
```

## Paper Summary

Per- and polyfluoroalkyl substances (PFAS) are persistent environmental contaminants that demand highly selective molecular recognition strategies for field-deployable detection. `beta`-Cyclodextrin-based field-effect transistor (FET) sensors demonstrate high sensitivity to perfluorooctanesulfonic acid (PFOS), achieving sub-ppt detection limits, yet exhibit limited selectivity in the presence of structurally similar surfactants such as sodium dodecyl sulfate (SDS). Here, we screen a synthetically accessible library of 1,629 functionalized `alpha`-, `beta`-, and `gamma`-cyclodextrins to quantify competitive binding thermodynamics using docking and all-atom molecular dynamics simulations. We identify host architectures with sub-nanomolar PFOS affinity and high selectivity, and use regression analysis to connect binding behavior to structural and electronic descriptors. Together, these results establish quantitative structure-selectivity relationships for cyclodextrin-based PFOS recognition and provide design principles for next-generation PFAS sensing materials.

## Repository Overview

The repository is organized into four main parts:

### [`chem_space_data/`](chem_space_data/)

Chemical-space datasets and structure libraries for the screened cyclodextrin hosts.

- `chem_space.pkl` stores the enumerated library and associated thermodynamic and descriptor fields.
- `chem_space_export.csv` is a tabular export of selected screened candidates for quick inspection outside Python.
- `chem_space_pdb_files/` contains the candidate 3D structures as PDB files.
- `prim_cleaved_structs/` contains reference alpha-, beta-, and gamma-cyclodextrin scaffolds.
- `analyze_chem_space.py` provides a lightweight way to inspect the dataset contents.

This is the best starting point if you want to understand the screened design space or inspect individual candidates.

### [`metadynamics/`](metadynamics/)

Simulation setup files and helper scripts for the all-atom molecular dynamics and metadynamics workflows used to evaluate host-guest binding.

- example system directories such as `bcd-pfos/`, `bcd-sds/`, and `00464-pfos/`
- shared GROMACS and PLUMED templates in `common.files/`
- helper scripts for generating PLUMED plane definitions and cyclodextrin backbone restraints
- `ff-parameterize/` for guest force-field generation and topology preparation

This is the best starting point if you want to reproduce or inspect the simulation setup workflow.

### [`bayesianoptimization/`](bayesianoptimization/)

Gaussian-process and candidate-ranking utilities used to connect descriptors and MD-derived data to chemical-space prioritization.

- `morganKernel.py` defines additive GP models
- `training.py` contains the model-training loop
- `get_candidates_delta.ipynb` is the interactive analysis notebook
- `data/` stores the serialized training, MD, and candidate datasets used by the modeling workflow

This is the best starting point if you want to inspect the surrogate-modeling and candidate-selection components.

### [`lasso/`](lasso/)

Descriptor-analysis notebook for sparse, interpretable regression on the screened chemical space.

- `lasso.ipynb` loads `chem_space.pkl`
- computes RDKit molecular descriptors and charge-based features for substituent sets
- supports feature selection and regression-style analysis alongside the Bayesian-optimization workflow

This is the best starting point if you want a more interpretable descriptor-based model rather than the Gaussian-process workflow.

## How To Read The CSV

The file [`chem_space_data/chem_space_export.csv`](chem_space_data/chem_space_export.csv) is a compact export of selected probe designs. It currently contains 79 probe entries plus a header row.

The columns are:

- `probe ID`
  Repository-style identifier for the candidate, matching the numeric naming convention used elsewhere in the chemical-space data.
- `CD type`
  Cyclodextrin family for the probe, reported as `alpha-CD`, `beta-CD`, or `gamma-CD`.
- `primary`
  A string-encoded list of substituents on the primary face of the cyclodextrin. These entries are written as SMILES-like fragments inside a Python-style list.
- `secondary`
  A string-encoded list of substituents on the secondary face of the cyclodextrin, in the same format.
- `dG_md`
  MD-derived PFOS binding free energy stored as `[mean, uncertainty]`.
- `ddG_md`
  Relative selectivity-style free energy term, also stored as `[mean, uncertainty]`.
- `Kd_md`
  Dissociation constant for PFOS from the MD workflow, stored as `[mean, uncertainty]`.
- `Kd_SDS/Kd_PFOS`
  Selectivity ratio between SDS and PFOS binding, stored as `[mean, uncertainty]`. Larger values indicate stronger preference for PFOS over SDS.

Two formatting details are important:

- The `primary` and `secondary` columns are not plain text labels. They are serialized lists of substituent strings.
- The thermodynamic columns are not single numbers. Each cell is a two-element array, where the first value is the central estimate and the second value is the uncertainty.

Example interpretation of one row:

- `probe ID = 00001`
- `CD type = gamma-CD`
  This probe is built on a gamma-cyclodextrin scaffold.
- `primary = ["[Br]", ...]`
  Every primary-site substituent in that candidate is bromine.
- `secondary = ["[OH]", ...]`
  The secondary sites remain hydroxylated.
- `dG_md = [-31.52, 0.40]`
  Mean binding free energy of about `-31.5` with uncertainty `0.4` in the stored units.
- `Kd_SDS/Kd_PFOS = [22.5, 8.8]`
  SDS is predicted to bind more weakly than PFOS by roughly a factor of `22.5`, with the listed uncertainty.

## Suggested Entry Points

- Start with [`chem_space_data/README.md`](chem_space_data/README.md) to understand the screened library and stored structures.
- Read [`metadynamics/README.md`](metadynamics/README.md) for the simulation setup and analysis workflow.
- Read [`bayesianoptimization/README.md`](bayesianoptimization/README.md) for the GP modeling and prioritization workflow.
- Open [`lasso/lasso.ipynb`](lasso/lasso.ipynb) for descriptor extraction and sparse regression analysis.

## Notes

- The repository is focused on computational workflow components and intermediate data products rather than a polished software package.
- Several simulation scripts reflect the original HPC environment and may require path updates before reuse on another system.
