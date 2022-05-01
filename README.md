# NAMD-FEP

# If you are familar with NAMD and have good HPC resources, please consider using Charmm GUI RBF workflow instead of this one.

The primary goal of this repository is to calculate the difference of binding free energy of a pair of small molecules against the same protein target, i.e., the ∆∆G of binding, which is of significant importance in hit-to-lead drug discovery.

There are many vague tutorials with either organic compound only system, or too much theoretical discussion, but without a clear one of actually showing people how to actually do an FEP with regard to protein-ligand complex, so I decide to make one, what you need to do is download the "ipynb" format file, and launch that file with your Jupyter Notebook, so you can run an actual FEP workflow, with stability and repeatability.

The work here is, based on NAMD tutorial "A Tutorial on Alchemical Free Energy Perturbation Calculations in NAMD:" from http://www.ks.uiuc.edu/Training/Tutorials/; and the paper describing Feprepare web server  J. Chem. Inf. Model. 2021, 61, 9, 4131–4138 (https://pubs.acs.org/doi/10.1021/acs.jcim.1c00215);
Feprepare webserver https://feprepare.vi-seem.eu/

Free energy perturbation, basically, involves one protein target, with a hybridized ligand (merged from a pair of similar ligands of interests), then we calculate the energy difference when gradually turning off the interaction of the first ligand while turning on the interaction of the second. The reason why we do it slowly is a request of sampling strategy, you don't have two understand 100% before you could do it, just like you don't have to understand 100% of the mechanism of a chemical reaction before you could actually finish the reaction. But it is always good if you can.

![9999999999999999999999999999999](https://user-images.githubusercontent.com/75652473/146633817-a19cd8fc-3355-44c1-a50d-98c1e22caaaf.png)
The first image above is a hybridized ligand in the water system, while the next is hybridized ligand and protein in water, so what we do is we separately simulate these to systems, the ∆∆G then will be processed with ∆∆G = ∆G complex (second image) - ∆G solvent (first image). For a better explanation, refer to http://www.alchemistry.org/wiki/Example:_Relative_Binding_Affinity (The link do not talk about NAMD, but the fundamental theory are all the same)

The next image is an image of PDB 1MQ5 with a hybridized ligand, it will serve as one of the two inputs in the whole FEP calculation. The ligand topologies are generated from LigParGen webserver http://zarbi.chem.yyale.edu/ligpargen/, and the hybridization of the ligands and the overall input generation of this protein-ligand complex is done with help of Feprepare webserver https://feprepare.vi-seem.eu/.

![image](https://user-images.githubusercontent.com/75652473/146633202-94569a82-c2cf-457a-95c0-754dfee4d7ae.png)

# Just follow the step by step tutorial inside the Jupyter Notebook

# Usage

It is assumed you already compiled NAMD on your laptop (Or just use a binary version, i.e., a pre-compiled version) from https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD
```
git clone https://github.com/quantaosun/NAMD-FEP
```

```
cd NAMD-FEP
```
It is assumed you have installed Jupyter Notebook
```
Jupyter Notebook NAMD-FEP_local.ipynb
```
In a test run, examples of analysis would be like :

Analysis and result of the complex leg

![Screenshot from 2021-12-18 15-28-06](https://user-images.githubusercontent.com/75652473/146633332-b4f62f43-5a55-493d-a0ee-0bf797862681.png)

Analysis and result of the solvent leg

![Screenshot from 2021-12-18 15-27-46](https://user-images.githubusercontent.com/75652473/146633327-6e5e4e86-d76f-4758-aff3-78c31e51532d.png)

# Based on the above two images, the ∆∆G would be ΔΔG FEP = ΔG complex- ΔG solvent = -7.1 kcal/mol (This is just a test due to very short simulation time)

If you have access to google Colab, or any other open souced cloud platforms, with a NAMD installed, it is then possible for you to run the whole FEP process there instead of your laptop, with usually a faster performance in simulation speed.

 Pros and Cons, FEP is more accurate than docking, but it cost more time, and could only handle ligands with a similar scaffold that can be aligned.
