# NAMD-FEP

<img width="1489" height="590" alt="image" src="https://github.com/user-attachments/assets/a8d4088a-076c-4ca0-92fa-ca8235001ee7" />


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

# Local workflow (no FepPrepare dependency)

The supported direction for new projects is a hybrid local workflow:

1. Prepare and inspect the protein and aligned ligands locally with VMD/`psfgen`.
   `tools/psfgen_ligand.tcl` is a starting template, not a parameter generator.
2. Supply validated ligand topology and parameters (for example, locally generated
   CGenFF-compatible files). NAMD and VMD do not assign reliable small-molecule
   parameters automatically.
3. Generate independent windows from an existing NAMD FEP template:

   ```bash
   python3 tools/rbfe_workflow.py validate ligand.pdb ligand.str
   python3 tools/rbfe_workflow.py generate \
     --template tools/independent_window.namd \
     --output-dir windows/forward --windows 16
   ```

   The template must contain `@ALCH_LAMBDA@` and `@ALCH_LAMBDA2@` markers
   where the window-specific `alchLambda` and `alchLambda2` values belong.
4. Run one window at a time, or schedule a small independent batch on a V100.
   Four CPU workers and one CUDA device are the default:

   ```bash
   python3 tools/rbfe_workflow.py command \
     --namd /opt/NAMD_3.0/namd3 --config windows/forward/window_000.namd
   ```

   The command is printed for review; the utility does not launch jobs or
   require replica exchange.

This removes FepPrepare and CHARMM-GUI from the normal orchestration path while
retaining CHARMM-GUI as a fallback for membranes, unusual residues, or systems
that cannot be assembled locally. FepPrepare source code is not assumed to be
reusable; licensing must be verified before incorporating any code.

Before migrating a production system, validate both solvent and complex legs:
check topology completeness, forward/reverse agreement, replicate convergence,
and consistent energies against the existing T4L example. The local workflow
does not replace chemical validation of ligand parameters.

## Local ligand preparation

`tools/ligand_prepare.py` provides alignment, hybrid-PDB construction, and
NAMD FEP B-factor generation without FepPrepare:

```bash
python3 tools/ligand_prepare.py align \
  --reference ligand_a.pdb --mobile ligand_b.pdb --output ligand_b_aligned.pdb

python3 tools/ligand_prepare.py hybrid \
  --ligand-a ligand_a.pdb --ligand-b ligand_b_aligned.pdb \
  --output ligand_hybrid.pdb --fep-output ligand_hybrid.fep
```

The hybrid writer uses the B-factor column referenced by `alchCol B`:
common atoms are `0.00`, atoms disappearing from ligand A are `-1.00`, and
atoms appearing from ligand B are `1.00`. It preserves A records and appends
B-only records. The generated PDB still requires a compatible dual-topology
PSF and CHARMM/CGenFF parameter files; these scripts do not invent atom types,
charges, bonds, angles, or dihedrals. Atom matching is based on unique atom
name plus element, so inspect the hybrid before production.

# Following the tutorial without Maestro or the Feprepare web server

The PDF tutorial (`NAMD-FEP-tutorial.pdf`, 1MQ5 XLC -> 1MQ6 XLD) was written
around two proprietary/remote resources: Schrodinger Maestro for protein
preparation and ligand alignment, and the Feprepare web server
(https://feprepare.vi-seem.eu/) for building the dual-topology hybrid and the
NAMD inputs. Both are optional. Every step below maps a tutorial step to a
free, local, open-source replacement; only NAMD itself is a license download.

**Step 1 - protein preparation (replaces Maestro PrepWizard).** Download the
biological assembly from the PDB (for the tutorial: 1MQ5 and 1MQ6), keep the
relevant chains, and remove waters and co-factors except the reference ligand.
Instead of Maestro, use one of:

- CHARMM-GUI PDB Reader (free web) to fill missing residues, protonate, and
  generate a clean PDB/PSF, or
- `pdbfixer` (OpenMM, open source) plus `reduce` or VMD `psfgen` for
  hydrogens, or
- VMD `psfgen` directly on a curated PDB if no loops need rebuilding.

Whatever the tool, reproduce the Maestro outputs manually: check the
protonation states (the tutorial neutralized both ligands and recorded the
protonation-state penalties S1/S2 so that
`ddG = ddG_FEP + (S2 - S1)`), delete protein hydrogens if your ligand
pipeline requires it, and strip `CONECT` records from the bottom of the PDB
with any text editor instead of Sublime. Save `protein.pdb`, `XLC.pdb`
(reference), and `XLD.pdb` (mutation).

**Step 1b - ligand alignment (replaces Maestro "flexible ligand alignment" /
"superimpose structures").** Use the bundled alignment tool (Kabsch/Horn
superposition over atoms with matching names; only standard Python required):

```bash
python3 tools/ligand_prepare.py align \
  --reference XLC.pdb --mobile XLD.pdb --output XLD_aligned.pdb
```

RDKit (`rdMolAlign`) or Open Babel are free alternatives if you prefer a
chemistry-aware (substructure-based) alignment; what matters is that the
shared scaffold is tightly superimposed and that the reference ligand's
coordinates are untouched so the complex stays consistent with the protein.
Note the tool matches atoms by unique name + element: give the common core
identical atom names in both PDB files before aligning.

**Step 2 - ligand parameters.** LigParGen (free web server) remains the
simplest option used by the tutorial: upload `XLC.pdb` and `XLD_aligned.pdb`
with the correct total charge (0 in the tutorial) and download the RTF/TOP and
PRM files. Fully local alternatives: CGenFF via the free SilcsBio
`cgenff` binary or `charmm-gui.org` ligand modeler, or `parmed`/`antechamber`
plus a conversion to CHARMM format. These scripts never invent atom types,
charges, bonds, or dihedrals - validated parameters are your responsibility.

**Step 3 - hybrid ligand and FEP flags (replaces the Feprepare web server).**
Build the dual-topology PDB plus the `.fep` flag file locally:

```bash
python3 tools/ligand_prepare.py hybrid \
  --ligand-a XLC.pdb --ligand-b XLD_aligned.pdb \
  --output ligand_hybrid.pdb --fep-output ligand_hybrid.fep
```

This writes the B-factor column that `alchCol B` reads: `0.00` for common
atoms, `-1.00` for atoms that disappear (reference-only), `1.00` for atoms
that appear (mutation-only). You must still create the dual-topology PSF that
pairs these coordinates with the two ligand topologies: adapt
`tools/psfgen_ligand.tcl` (run with `vmd -dispdev text -e`) so that the
topology files for both ligands are loaded and the hybrid segment is written,
then solvate/ionize the solvent leg and the protein-complex leg with VMD
(`solvate`, `autoionize`) or CHARMM-GUI. If your system defeats local assembly
(membranes, unusual residues), CHARMM-GUI's RBFE workflow is the fallback, as
noted above.

**Steps 4-6 - NAMD inputs and simulation.** Nothing here required Feprepare:
the repository already contains the four modified configuration files used by
the tutorial (`*_nvt_equil_test.namd`, `*_npt_equil_test.namd`,
`*_md_forward_test.namd`, `*_md_backward_test.namd` for the complex and
solvent legs) together with `toppar_modified.zip` and `parameter_patch2.txt`.
Copy them into your `complex/` and `solvent/` directories, update the
`structure`/`coordinates`/`parameters` paths, the periodic cell vectors (from
the solvated systems), and the step counts. The notebook
`NAMD-FEP_local.ipynb` then performs the identical short test run
(`namd2 nvt_equil_test.namd > nvt_test.log`, then NPT, forward, backward) for
both legs; scale the step counts up for production. If you prefer
Feprepare-style independent lambda windows, generate them from the template
instead of running the single sweep file:

```bash
python3 tools/rbfe_workflow.py generate \
  --template tools/independent_window.namd \
  --output-dir windows/forward --windows 16
python3 tools/rbfe_workflow.py command \
  --namd /path/to/namd3 --config windows/forward/window_000.namd
```

For cluster submission, any scheduler works (the tutorial's PBS example, or
Slurm); FileZilla is only a file-transfer convenience - `scp`/`rsync` do the
same job.

**Step 7 - analysis.** VMD's ParseFEP plugin (`Extensions > Analysis >
Analyze FEP Simulation`) is free and works without Maestro/Feprepare: select
the forward and backward `.fepout` files per leg, set the simulation
temperature, and take `ddG_FEP = dG_complex - dG_solv`. See
`T4L_RBFE/output_files/` for reference outputs and the expected BAR table
format. Free non-VMD options include `alchemlyb` (parsing + BAR/MBAR in
Python) or the bundled `Decomp_barchart.ipynb` notebook for per-residue
decomposition plots. Finally apply the protonation-state correction:
`ddG = ddG_FEP + (S2 - S1)`.
