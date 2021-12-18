# NAMD-fep

This is a detailed workflow, based on NAMD tutorial "A Tutorial on Alchemical Free Energy Perturbation Calculations in NAMD:" from http://www.ks.uiuc.edu/Training/Tutorials/; and the paper describing Feprepare web server https://pubs.acs.org/doi/abs/10.1021/acs.jcim.1c00215;
Feprepare web server https://feprepare.vi-seem.eu/

If you have access to google colab, or any other open souced cloud platforms, with a NAMD installed, it is then possible for you to run the whole FEP process there instead of your loptop, with usually a faster performance in simulation speed.

The most important file of this tutorial is the jupyter notebook, inside there are 8 configuration files, 4 for complex leg, 4 for solvent leg, to finish all the testing simulation. For HPC or cloud NAMD simulation, use the 3 configuration files listed above instead. (you can generate 4th one for HPC/Cloud by yourself)

# Just follow the step by step tutorial inside the jupyter notebook

# Usage

It is assumed you already compiled namd on your labtop (Or just use a binary version, i.e., a pre-cmopiled version) from https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD
```
git clone https://github.com/quantaosun/NAMD-FEP
```

```
cd NAMD-FEP
```
It is assumed you have installed jupyter notebook
```
jupyter notebook NAMD-FEP_local.ipynb
```
The jupyter notebook will guide you through a very short FEP calculation and analysis, for a practical use, please change the simulation steps to much larger number that, and submit the job to a HPC cluster or clould platforms. A fair fep simulation would cost a day or two days, depending the size of your protein, and mostly on how good is the computing resources you could get access to.
