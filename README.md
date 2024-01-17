# Installation
```
pip install git+https://github.com/jonasbarth/DIAMOnD
```
The repository contains code for the DIAMOnD, DiaBLE, Diffusion algorithms.

# DIAMOnD
```
diamond("<your network file>", "<your seed genes file>", <number of genes to predict>)
```
# DiaBLE
```
from diamond import diable

diable("<your network file>", "<your seed genes file>", <number of genes to predict>)
```

# Diffusion
This runs the cytoscape diffusion algorithm. Default time value is `time=0.01`. It will return a pandas dataframe with columns `gene, heat`, sorted in descending order of the `heat` column. 
```
from diamond import diffusion

diffusion("<your network file>", "<your seed genes file>", <number of genes to predict>, time=<time>)
```
