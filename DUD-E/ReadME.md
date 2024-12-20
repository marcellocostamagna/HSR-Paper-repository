## Description

This folder contains scripts for conducting enrichment factor experiments of the Database of Useful Decoys: Enhanced (DUD-E) using three different packages: *HSR*, *RDKIT*, and *ODDT*.

## Instructions

**REQUIREMENT**: Ensure you have downloaded the [DUD-E database](https://dude.docking.org/) and placed it in this directory, it is named `all.tar.gz`.

Unzip the file by running the command:

```bash
tar -xvzf all.tar.gz && find . -name "*.sdf.gz" -exec gunzip {} \;
```

A folder name `all` should  have been created inside the directory. If you choose to rename this folder, you must also update the folder name in each script accordingly.

*HSR* 

Run the command:

```bash
python dude_hsr_enrichments.py
```

*RDKIT* 

Run the command:

```bash
python dude_rdkit_enrichments.py
```
*ODDT* 

Run the command:

```bash
python dude_oddt_enrichments.py
```
