# Installation
```
pip install git+https://github.com/jonasbarth/DIAMOnD
```

# DIAMOnD

DIAMOnD.py runs the DIAMOnD algorithm as described in
 
 A DIseAse MOdule Detection (DIAMOnD) Algorithm derived from a
 systematic analysis of connectivity patterns of disease proteins in
 the Human Interactome. PlOS Comp Bio (in press), 2015.

by Susan Dina Ghiassian, Joerg Menche & Albert-Laszlo Barabasi

The DIAMOnD website can be found at: 

Instruction to use the source code:
1. Download the code.
2. Make sure you are in the main directory where the code is.
3. Run the following.</br>
 <em><pre>python3 DIAMOnD.py  network_file seed_file  n  alpha(optional)  outfile_name(optional)</pre></em>

# -------------------

Directory Example

contains two input files:
1. seed_genes.txt (list of genes associated with a phenotype of interest) 
2. PPI.txt (Protein-protein interaction network. note that gene IDs should be consistent in the two input files)</br>
The following command will generate the first 100 DIAMOnD nodes and save them in a file) </br>
<em><pre>python3  DIAMOnD.py  Example/PPI.txt  Example/seed_genes.txt  100</pre></em>

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
