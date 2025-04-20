## Description

This folder contains scripts to test the speed performance of HSR for generating fingerprints and performing similarity comparisons.

## Instructions

### Step 1: Collect ZINC20 SDF URLs

Run the following command:

```bash
pymol save_zinc_sdf_links.py
```

This script collects all the URLs of `sdf.gz` files available in the ZINC20 database. The process may take a few minutes but significantly speeds up the subsequent steps.

### Step 2: Run a Single Benchmark Test

Once the `zinc_sdf_links.txt` file has been generated, run:

```bash
pymol hsr_speed_test.py
```
This script selects a random subset of molecules (the number is defined in the script) from the ZINC20 database using the previously generated list of URLs. It returns information about the time taken to generate the HSR fingerprints and perform similarity comparisons.

### Step 3: Run Multiple Benchmark Repetitions

To generate average performance statistics over multiple runs of `hsr_speed_test.py`, use:

```bash
pymol run_hsr_speed_tests.py
```

By default this command will run 3 repetitions. To specify a different number of repetitions—for example, 10—use the -n option:

```bash
pymol run_hsr_speed_tests.py -n 10
```