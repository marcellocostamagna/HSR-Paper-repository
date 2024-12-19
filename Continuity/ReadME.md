## Description

This folder contains scripts for experiments designed to show HSR greater robustness to discontinuity issues compared to USR. The similarity comparison between two molecule in two different conformations are conducted with both *HSR* and *USR* similarity meausures.

## Instructions

*HSR* 

Run the command:

```bash
python HSR.py
```

*USR* 

Run the command:

```bash
python USR.py
```

*HSR visualization*

To obtain the visualization of the reference points for the two sets of conformers run the command:

```bash
pymol HSR_conf_reference_points.py
```
To change the conformers modify the following section of the script `HSR_conf_reference_points.py`:

```python
# Choose what confermers to process: A or B
conformer = 'A' # 'B'
# Process both molecules
if conformer == 'A':
    # Conformer A
    file1 = f'{cwd}/conformers/conformersA_H.sdf'
    file2 = f'{cwd}/conformers/conformersA_Met.sdf'
elif conformer == 'B':
    # Conformer B
    file1 = f'{cwd}/conformers/conformersB_H.sdf'
    file2 = f'{cwd}/conformers/conformersB_Met.sdf'
```



