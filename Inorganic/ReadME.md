## Description

This folder contains scripts for experiments designed to show HSR ability to compare inorganic and organometallic molecules in two main contexts. 
1- Three pairs of organometallic isomers are compared with the methods: *HSR*, *USR* and *USRCAT*.
2- Three structural organometallic analogues are compared using the same methods.

Additionally, a separate script is provided to illustrate USRCAT’s dependency on molecular connectivity.

## Instructions

### Organometallic Isomers

*HSR*

Run the command:

```bash
python inorganic_molecules_HSR.py
```

*USR and USRCAT*

Run the command:

```bash
python inorganic_molecules_USR.py
```
### Structural analogues

The user must indipendently retreive the sdf files of the heaammine complexes present in the CSD entries: ADIYES, BIWHEY, and CAFWEM, and name them `Co_3.sdf`, `Ir_3.sdf`, and `Co_2.sdf`, respectively.

Once the files have been placed in the currect directory run the command:

```bash
python structural_analogues.py
```

### USRCAT's connectivity dependency*

Run the command:

```bash
python USRCAT_connectivity_dependency.py
```