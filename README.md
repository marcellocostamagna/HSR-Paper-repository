# HSR-Paper-repository [![DOI](https://zenodo.org/badge/756043970.svg)](https://doi.org/10.5281/zenodo.14631654) 

This repository collects the scripts used to generate results shown in the HSR paper (TODO:Link).

## Getting Started

### Environment Setup

Before running any experiments, make sure to create the appropriate conda environment by running the following command:

```bash
conda env create -f environment.yml
conda activate HSR_results
```

You can test the correct build by running:

```bash
hsr -v
```

### Run experiments

This repository includes several folders containing data and scripts to reproduce the results presented in the HSR paper[TODO:link to paper].

To run any experiment move to the respective folder:

```bash
cd folder_name
```

and follow the instructions provided there.

**REQUIREMENT** : To run the DUD-E experiments you need to download the default [DUD-E database](https://dude.docking.org/) and move it to the DUD-E folder.


### Additional files

- usr.py: in-house implementation of the Ultrafast Shape Recognition (USR) method
- csr.py: in-house implementation of the Chiral Shape Recognition (CSR) method
- usr_optiso.py: in-house implementation of the USR:OptIso method
- perturbations.py: scripts collecting diffferent functions to perform trasformations such as
    - Translation of molecule coordinates
    - Rotations of molecule coordinates
    - Reflection of molecule coordinates (used to generate enantiomers)
    - Scaling
    - Permutation
    - etc...

