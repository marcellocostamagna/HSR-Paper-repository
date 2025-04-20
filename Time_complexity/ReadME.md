## Description

This folder contains a script to evaluate the increase in runtime complexity when generating HSR fingerprints with additional atomic features.

The experiment measures the average time required to generate fingerprints in 3D, 4D, 5D, and 6D feature spaces, using the same molecule. The goal is to quantify the performance cost introduced by including more atomic descriptors.

To reproduce the results averaged over 1000 runs, simply execute:

```bash
python time_complexity.py
```

The output will report the average fingerprint generation time for each dimensionality, along with the relative percentage increase compared to the 3D baseline.