# NAMD-FEP

The most important file of this tutorial is the jupyter notebook, inside there are 8 configuration files, 4 for complex leg, 4 for solvent leg, to finish all the testing simulation.

The primary goal of this repository is to calculate the difference of binding free energy of a pair of small molecules against the same protein target, i.e., the delt delt G of binding, which is of significant importance in hit-to-lead drug discovery.

 Pros and Cons, FEP is more accurate than docking, but it cost more time, and could only handle ligands with similar scaffold that can be aligned.


This is a detailed workflow, based on NAMD tutorial "A Tutorial on Alchemical Free Energy Perturbation Calculations in NAMD:" from http://www.ks.uiuc.edu/Training/Tutorials/; and the paper describing Feprepare web server  J. Chem. Inf. Model. 2021, 61, 9, 4131–4138 (https://pubs.acs.org/doi/10.1021/acs.jcim.1c00215);
Feprepare web server https://feprepare.vi-seem.eu/


Free energy pertubation,basicly, it involves one protein target, with a hybridized ligand (merged from a pair of similar ligands of interests), then we calculate the energy difference when gradually turning off interaction of first ligand while turning on the interaction of the second. The reason why we do it slowly is a request of sampling strategy, you don't have to understand 100% before you could do it, just like you don't have to understand 100% of the mechansim of a chemical reaction before you could actually finish the reaction. But it is always good if you can. 

![9999999999999999999999999999999](https://user-images.githubusercontent.com/75652473/146633817-a19cd8fc-3355-44c1-a50d-98c1e22caaaf.png)
The first image above is a hybrized ligand in water system, while the next is hybridized ligand and protein in water, so what we do is we separately simulate these to systems, the ddG then will be processed with ddG = dG complex (second image) - dG solvent (first image). For a better explanation, refer to http://www.alchemistry.org/wiki/Example:_Relative_Binding_Affinity (The link do not talk about NAMD, but the fundamental theory are all the same)

The next image is an image of PDB 1MQ5 with a hybridized ligand, it will serve as the input as the FEP calculation. The ligand topologies are generated from ligpargen web server http://zarbi.chem.yyale.edu/ligpargen/, and the hybridization of the ligands and the overall input generation of this protein-ligand complex is done with help of Feprepare web server https://feprepare.vi-seem.eu/.

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

Based on above two images, the delt delt G would be ΔΔG FEP = ΔG complex- ΔG solvent = -7.1 kcal/mol (This is just a test due to very short simulation time)

If you have access to google colab, or any other open souced cloud platforms, with a NAMD installed, it is then possible for you to run the whole FEP process there instead of your loptop, with usually a faster performance in simulation speed.
