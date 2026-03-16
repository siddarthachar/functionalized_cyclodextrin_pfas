# Bayesian Optimization

This directory contains the Gaussian-process and candidate-selection pieces used to rank cyclodextrin variants from precomputed molecular descriptors and MD-derived binding data.

## What Is Here

### Core code

- `morganKernel.py`
  Defines additive Gaussian-process models built with GPyTorch and BoTorch. The models split the feature vector into a "primary" block and a "secondary" block, then apply separate Matern kernels and add them together.
- `training.py`
  Contains `train_model(...)`, a small training loop for exact GP models using Adam and `ExactMarginalLogLikelihood`.
- `get_candidates_delta.ipynb`
  Notebook for interactive analysis and candidate selection. This is the main place to combine the stored data objects with the GP models and inspect predicted rankings.
- `README`
  Older placeholder note. The more complete directory documentation is in this file.

### Data subdirectory

`data/` holds the serialized inputs used by the notebook and training code.

- `data/training_data.pkl`
  Training-ready data object used to fit the GP model.
- `data/md_data.pkl`
  MD-derived measurements used to build or validate the training set.
- `data/chem_space.pkl`
  Full screened chemical space with descriptors and metadata.
- `data/candidates.pkl`
  Candidate subset selected from the larger chemical space.
- `data/dG_data_pfoa.pkl`
  Free-energy dataset for PFOA-related runs.
- `data/dG_data_pfoa_2.pkl`
  Alternate or updated PFOA free-energy dataset.
- `data/chem_space_structures/`
  PDB files for the screened chemical-space structures. Filenames are numeric IDs such as `00000.pdb`.
- `data/prim_cleaved_structs/`
  Reference host structures with primary-side groups removed:
  `acd_prim_removed.pdb`, `bcd_prim_removed.pdb`, and `gcd_prim_removed.pdb`.
- `data/labels/`
  Label assets used by the data-generation workflow.
- `data/README.md`
  Short note indicating that `chem_space.pkl` is generated from a fingerprint-generation script.

## How The Code Fits Together

The intended flow is:

1. Prepare the descriptor and MD data objects in `data/`.
2. Load them in `get_candidates_delta.ipynb`.
3. Instantiate one of the additive GP models from `morganKernel.py`.
4. Train the model with `train_model(...)` from `training.py`.
5. Score the full chemical space in `data/chem_space.pkl`.
6. Save or inspect shortlisted candidates, for example in `data/candidates.pkl`.

The model classes assume the input feature vector is concatenated in a fixed order:

- the first `primary_dim` columns correspond to the primary substituent features
- the next `secondary_dim` columns correspond to the secondary substituent features

If you change the descriptor layout, you must update those dimensions when building the model.

## How To Use It

### 1. Install the Python dependencies

At minimum, the code in this directory expects:

- `torch`
- `gpytorch`
- `botorch`
- `jupyter` for the notebook workflow

You will also need whatever packages were used to create the pickle files in `data/`.

### 2. Open the notebook workflow

From the repository root:

```bash
cd bayesianoptimization
jupyter notebook get_candidates_delta.ipynb
```

Use the notebook when you want the full interactive workflow: loading the pickles, fitting the GP, and exploring candidate rankings.

### 3. Reuse the GP model in a script

Minimal training example:

```python
import torch
import gpytorch
from gpytorch.mlls import ExactMarginalLogLikelihood

from morganKernel import AdditiveGPModel_botorch
from training import train_model

train_x = torch.randn(20, 16)
train_y = torch.randn(20)
likelihood = gpytorch.likelihoods.GaussianLikelihood()

model = AdditiveGPModel_botorch(
    train_x=train_x,
    train_y=train_y,
    likelihood=likelihood,
    primary_dim=8,
    secondary_dim=8,
)

mll = ExactMarginalLogLikelihood(likelihood, model)
model = train_model(model, mll, train_x, train_y, max_iter=500, lr=0.1)
```

### 4. Work with the stored structures

If you need the 3D structures associated with a chemical-space ID, look in:

```bash
bayesianoptimization/data/chem_space_structures/<ID>.pdb
```

Example:

```bash
ls bayesianoptimization/data/chem_space_structures/00042.pdb
```

## Notes And Assumptions

- This directory is centered on model definition and candidate ranking, not raw descriptor generation.
- The notebook is likely the main analysis entry point; the Python modules are lightweight helpers around that workflow.
- The pickles are binary data products and are not self-describing from the command line. If you modify their schema, update this README and any loading code in the notebook.
