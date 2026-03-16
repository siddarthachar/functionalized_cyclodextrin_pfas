# Chem Space Data

This directory stores the enumerated cyclodextrin chemical space and the structure files associated with it. It is the data source that feeds the screening and Bayesian-optimization workflow.

## What Is Here

- `chem_space.pkl`
  Main serialized chemical-space object. This is the file used by the analysis script in this directory and by the Bayesian-optimization assets elsewhere in the repository.
- `analyze_chem_space.py`
  Lightweight inspection script for `chem_space.pkl`. It prints counts, data coverage, a few example entries, and top-ranked systems by selected fields.
- `chem_space_pdb_files/`
  PDB structures for the enumerated candidates. The filenames are numeric IDs such as `00000.pdb`, `00001.pdb`, and so on.
- `prim_cleaved_structs/`
  Reference cyclodextrin scaffolds with primary-side substituents removed:
  `acd_prim_removed.pdb`, `bcd_prim_removed.pdb`, and `gcd_prim_removed.pdb`.

## What The Data Represents

From the analysis code, each entry in `chem_space.pkl` is keyed by a system identifier and includes fields such as:

- `CD`
  Cyclodextrin family.
- `primary`
  Primary-side substitution description.
- `secondary`
  Secondary-side substitution description.
- `dG_md`
  MD-derived free energy.
- `ddG_md`
  Relative free energy.
- `Kb_md`
  Binding constant.
- `Kd_md`
  Dissociation constant.
- `Kd_SDS/Kd_PFOS`
  Selectivity-style ratio field.
- `delta_dG`
- `delta_ddG`

Not every field is necessarily populated for every entry; `analyze_chem_space.py` explicitly reports coverage.

## How To Use It

### Inspect the chemical space

From this directory:

```bash
cd chem_space_data
python analyze_chem_space.py
```

Or pass an explicit path:

```bash
python analyze_chem_space.py chem_space.pkl
```

The script reports:

- total number of systems
- how many systems have observed MD values versus generated-only values
- counts by cyclodextrin type
- field coverage across the dataset
- example entries
- top entries by `dG_md`, `ddG_md`, `Kb_md`, and `Kd_SDS/Kd_PFOS`

### Retrieve a candidate structure

If you know the numeric candidate ID, open the corresponding PDB file:

```bash
ls chem_space_pdb_files/00042.pdb
```

This is useful when you want to move from a ranked candidate in the data object to a concrete 3D structure for follow-up simulation or visualization.

### Reuse the reference scaffolds

The files in `prim_cleaved_structs/` are helpful when generating or modifying candidate structures by cyclodextrin family:

- `acd_prim_removed.pdb`
- `bcd_prim_removed.pdb`
- `gcd_prim_removed.pdb`

## Relationship To Other Directories

- `chem_space.pkl` here is the same kind of dataset referenced in [`bayesianoptimization/data`](../bayesianoptimization/data/).
- The structure library in `chem_space_pdb_files/` mirrors the same candidate-ID style used by `bayesianoptimization/data/chem_space_structures/`.
- This directory is data-focused; model fitting and candidate ranking happen in [`bayesianoptimization`](../bayesianoptimization/).

## Requirements

The inspection script only needs the Python standard library:

- `pathlib`
- `pickle`
- `sys`

That makes this the easiest place to sanity-check the chemical-space dataset before using the heavier Bayesian-optimization workflow.
