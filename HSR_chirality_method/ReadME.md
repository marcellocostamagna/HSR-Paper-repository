## Description

This folder contains scripts to generate the images used fo Figure 3 (TODO: check consistency) in the HSR Paper (TODO: link), which show HSR method to distinguish chiral molecules.

## Instructions

Run the command:

```bash
pymol enantiomers.py
```

To generate the different enantiomers, with or without chirality detection, modify the beginnig script `enantiomers.py`:

```python
# MOLECULE SELECTION AND CHIRALITY
# Select one or the other to visualize an enantiomer at the time
# file_name = 'helicene_M.sdf'
file_name = 'helicene_P.sdf'

# Select chirality to be True or False
chirality = True
# chirality = False
```

In this case, the script will generate the enantiomer P, showing the Principal Components obtained through the chirality detection method.