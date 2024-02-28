## Description

This folder contains script of an experiment designed to show the dependency of the HSR similarity score to the positon of a feature change (specifically the change of an hydrogen to a deuterium) and the script tp generate the image showing the result (Figure 8 in th HSR paper (TODO:consictency and link)).


## Instructions

TO obtain the similarity scores of the comparison betwen the original molecules and all the ones showing a deuterium in different positions run the command:

```bash
python position_dependency.py
```

To get the relative image run the command:


```bash
pymol heat_map.py
```
