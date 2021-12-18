# NAMD-fep

# The primary goal of this repository is to calculate the difference of binding free energy of a pair of small molecules against the same protein target, i.e., the delt delt G of binding, which is of significant importance in hit-to-lead drug discovery.

This is a detailed workflow, based on NAMD tutorial "A Tutorial on Alchemical Free Energy Perturbation Calculations in NAMD:" from http://www.ks.uiuc.edu/Training/Tutorials/; and the paper describing Feprepare web server https://pubs.acs.org/doi/abs/10.1021/acs.jcim.1c00215;
Feprepare web server https://feprepare.vi-seem.eu/

If you have access to google colab, or any other open souced cloud platforms, with a NAMD installed, it is then possible for you to run the whole FEP process there instead of your loptop, with usually a faster performance in simulation speed.

The most important file of this tutorial is the jupyter notebook, inside there are 8 configuration files, 4 for complex leg, 4 for solvent leg, to finish all the testing simulation. For HPC or cloud NAMD simulation, use the 3 configuration files listed above instead. (you can generate 4th one for HPC/Cloud by yourself)


![image](https://user-images.githubusercontent.com/75652473/146633202-94569a82-c2cf-457a-95c0-754dfee4d7ae.png)


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
In a test run, examples of analysis would be like :
![Screenshot from 2021-12-18 15-27-46](https://user-images.githubusercontent.com/75652473/146633327-6e5e4e86-d76f-4758-aff3-78c31e51532d.png)
![Screenshot from 2021-12-18 15-28-06](https://user-images.githubusercontent.com/75652473/146633332-b4f62f43-5a55-493d-a0ee-0bf797862681.png)

so, the delt delt G would be ΔΔG FEP = ΔG complex- ΔG solvent = -7.1 kcal/mol (This is just a test due to very short simulation time)

