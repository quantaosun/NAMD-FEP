#!/bin/bash

# Lambda values from replica exchange setup
LAMBDA_VALUES=(0.0000 0.045 0.090 0.14546 0.22425 0.30303 0.38182 0.46061 0.5394 0.61819 0.6970 0.77576 0.85455 0.910 0.955 1.0000)

echo "Generating equilibration files..."

# Forward equilibration (λ=0.0000)
cat > equ_solv_forward.namd << 'EOF'
# EQUILIBRATION FOR FORWARD FEP (λ=0.0000)
structure               ligand.psf
coordinates             ligand.pdb
source                  ligand.str

set temp                300
set outputname          ligand_eq_forward

temperature             $temp
outputName              $outputname
outputEnergies          100
outputPressure          100

paraTypeCharmm          on;
parameters              ../../toppar/par_all36m_prot.prm
parameters              ../../toppar/par_all36_na.prm
parameters              ../../toppar/par_all36_carb.prm
parameters              ../../toppar/par_all36_lipid.prm
parameters              ../../toppar/par_all36_cgenff.prm
parameters              ../../toppar/par_interface.prm
parameters              ../../toppar/toppar_all36_moreions.str
parameters              ../../toppar/toppar_all36_nano_lig.str
parameters              ../../toppar/toppar_all36_nano_lig_patch.str
parameters              ../../toppar/toppar_all36_synthetic_polymer.str
parameters              ../../toppar/toppar_all36_synthetic_polymer_patch.str
parameters              ../../toppar/toppar_all36_polymer_solvent.str
parameters              ../../toppar/toppar_water_ions.str
parameters              ../../toppar/toppar_dum_noble_gases.str
parameters              ../../toppar/toppar_ions_won.str
parameters              ../../toppar/cam.str
parameters              ../../toppar/toppar_all36_prot_arg0.str
parameters              ../../toppar/toppar_all36_prot_c36m_d_aminoacids.str
parameters              ../../toppar/toppar_all36_prot_fluoro_alkanes.str
parameters              ../../toppar/toppar_all36_prot_heme.str
parameters              ../../toppar/toppar_all36_prot_na_combined.str
parameters              ../../toppar/toppar_all36_prot_retinol.str
parameters              ../../toppar/toppar_all36_prot_model.str
parameters              ../../toppar/toppar_all36_prot_modify_res.str
parameters              ../../toppar/toppar_all36_na_nad_ppi.str
parameters              ../../toppar/toppar_all36_na_rna_modified.str
parameters              ../../toppar/toppar_all36_lipid_sphingo.str
parameters              ../../toppar/toppar_all36_lipid_archaeal.str
parameters              ../../toppar/toppar_all36_lipid_bacterial.str
parameters              ../../toppar/toppar_all36_lipid_cardiolipin.str
parameters              ../../toppar/toppar_all36_lipid_cholesterol.str
parameters              ../../toppar/toppar_all36_lipid_dag.str
parameters              ../../toppar/toppar_all36_lipid_inositol.str
parameters              ../../toppar/toppar_all36_lipid_lnp.str
parameters              ../../toppar/toppar_all36_lipid_lps.str
parameters              ../../toppar/toppar_all36_lipid_mycobacterial.str
parameters              ../../toppar/toppar_all36_lipid_miscellaneous.str
parameters              ../../toppar/toppar_all36_lipid_model.str
parameters              ../../toppar/toppar_all36_lipid_prot.str
parameters              ../../toppar/toppar_all36_lipid_tag.str
parameters              ../../toppar/toppar_all36_lipid_yeast.str
parameters              ../../toppar/toppar_all36_lipid_hmmm.str
parameters              ../../toppar/toppar_all36_lipid_detergent.str
parameters              ../../toppar/toppar_all36_lipid_ether.str
parameters              ../../toppar/toppar_all36_lipid_oxidized.str
parameters              ../../toppar/toppar_all36_carb_glycolipid.str
parameters              ../../toppar/toppar_all36_carb_glycopeptide.str
parameters              ../../toppar/toppar_all36_carb_imlab.str
parameters              ../../toppar/toppar_all36_label_spin.str
parameters              ../../toppar/toppar_all36_label_fluorophore.str
parameters              ../morph.prm


exclude                 scaled1-4
1-4scaling              1.0
switching               on
cutoff                  12.0
switchdist              10.0
pairlistdist            16.0
stepspercycle           20

timestep                2.0
rigidBonds              water
nonbondedFreq           1
fullElectFrequency      1

reassignFreq            500
reassignTemp            $temp

cellBasisVector1     $a   0.0   0.0
cellBasisVector2    0.0    $b   0.0
cellBasisVector3    0.0   0.0    $c
cellOrigin          0.0   0.0 $zcen
wrapAll                 on

PME                     yes
PMEGridSpacing          1.0

langevin                on
langevinDamping         1.0
langevinTemp            $temp
langevinPiston          on
langevinPistonTarget    1.01325
langevinPistonPeriod    100.0
langevinPistonDecay     100.0
langevinPistonTemp      $temp

alch                    on
singleTopology          on
alchType                FEP
alchFile                ligand.pdb
alchCol                 B
alchLambda              0.0000
alchLambda2             0.0000

minimize                200
run                     10000
EOF

# Backward equilibration (λ=1.0000)
cat > equ_solv_backward.namd << 'EOF'
# EQUILIBRATION FOR BACKWARD FEP (λ=1.0000)
structure               ligand.psf
coordinates             ligand.pdb
source                  ligand.str

set temp                300
set outputname          ligand_eq_backward

temperature             $temp
outputName              $outputname
outputEnergies          100
outputPressure          100

paraTypeCharmm          on;
parameters              ../../toppar/par_all36m_prot.prm
parameters              ../../toppar/par_all36_na.prm
parameters              ../../toppar/par_all36_carb.prm
parameters              ../../toppar/par_all36_lipid.prm
parameters              ../../toppar/par_all36_cgenff.prm
parameters              ../../toppar/par_interface.prm
parameters              ../../toppar/toppar_all36_moreions.str
parameters              ../../toppar/toppar_all36_nano_lig.str
parameters              ../../toppar/toppar_all36_nano_lig_patch.str
parameters              ../../toppar/toppar_all36_synthetic_polymer.str
parameters              ../../toppar/toppar_all36_synthetic_polymer_patch.str
parameters              ../../toppar/toppar_all36_polymer_solvent.str
parameters              ../../toppar/toppar_water_ions.str
parameters              ../../toppar/toppar_dum_noble_gases.str
parameters              ../../toppar/toppar_ions_won.str
parameters              ../../toppar/cam.str
parameters              ../../toppar/toppar_all36_prot_arg0.str
parameters              ../../toppar/toppar_all36_prot_c36m_d_aminoacids.str
parameters              ../../toppar/toppar_all36_prot_fluoro_alkanes.str
parameters              ../../toppar/toppar_all36_prot_heme.str
parameters              ../../toppar/toppar_all36_prot_na_combined.str
parameters              ../../toppar/toppar_all36_prot_retinol.str
parameters              ../../toppar/toppar_all36_prot_model.str
parameters              ../../toppar/toppar_all36_prot_modify_res.str
parameters              ../../toppar/toppar_all36_na_nad_ppi.str
parameters              ../../toppar/toppar_all36_na_rna_modified.str
parameters              ../../toppar/toppar_all36_lipid_sphingo.str
parameters              ../../toppar/toppar_all36_lipid_archaeal.str
parameters              ../../toppar/toppar_all36_lipid_bacterial.str
parameters              ../../toppar/toppar_all36_lipid_cardiolipin.str
parameters              ../../toppar/toppar_all36_lipid_cholesterol.str
parameters              ../../toppar/toppar_all36_lipid_dag.str
parameters              ../../toppar/toppar_all36_lipid_inositol.str
parameters              ../../toppar/toppar_all36_lipid_lnp.str
parameters              ../../toppar/toppar_all36_lipid_lps.str
parameters              ../../toppar/toppar_all36_lipid_mycobacterial.str
parameters              ../../toppar/toppar_all36_lipid_miscellaneous.str
parameters              ../../toppar/toppar_all36_lipid_model.str
parameters              ../../toppar/toppar_all36_lipid_prot.str
parameters              ../../toppar/toppar_all36_lipid_tag.str
parameters              ../../toppar/toppar_all36_lipid_yeast.str
parameters              ../../toppar/toppar_all36_lipid_hmmm.str
parameters              ../../toppar/toppar_all36_lipid_detergent.str
parameters              ../../toppar/toppar_all36_lipid_ether.str
parameters              ../../toppar/toppar_all36_lipid_oxidized.str
parameters              ../../toppar/toppar_all36_carb_glycolipid.str
parameters              ../../toppar/toppar_all36_carb_glycopeptide.str
parameters              ../../toppar/toppar_all36_carb_imlab.str
parameters              ../../toppar/toppar_all36_label_spin.str
parameters              ../../toppar/toppar_all36_label_fluorophore.str
parameters              ../morph.prm


exclude                 scaled1-4
1-4scaling              1.0
switching               on
cutoff                  12.0
switchdist              10.0
pairlistdist            16.0
stepspercycle           20

timestep                2.0
rigidBonds              water
nonbondedFreq           1
fullElectFrequency      1

reassignFreq            500
reassignTemp            $temp

cellBasisVector1     $a   0.0   0.0
cellBasisVector2    0.0    $b   0.0
cellBasisVector3    0.0   0.0    $c
cellOrigin          0.0   0.0 $zcen
wrapAll                 on

PME                     yes
PMEGridSpacing          1.0

langevin                on
langevinDamping         1.0
langevinTemp            $temp
langevinPiston          on
langevinPistonTarget    1.01325
langevinPistonPeriod    100.0
langevinPistonDecay     100.0
langevinPistonTemp      $temp

alch                    on
singleTopology          on
alchType                FEP
alchFile                ligand.pdb
alchCol                 B
alchLambda              1.0000
alchLambda2             1.0000

minimize                200
run                     10000
EOF

echo "Generating base FEP configuration files..."

# Base forward FEP config
cat > fep_base_forward.conf << 'EOF'
structure               ligand.psf
coordinates             ligand.pdb
source                  ligand.str

set temp                300

paraTypeCharmm          on;
parameters              ../../toppar/par_all36m_prot.prm
parameters              ../../toppar/par_all36_na.prm
parameters              ../../toppar/par_all36_carb.prm
parameters              ../../toppar/par_all36_lipid.prm
parameters              ../../toppar/par_all36_cgenff.prm
parameters              ../../toppar/par_interface.prm
parameters              ../../toppar/toppar_all36_moreions.str
parameters              ../../toppar/toppar_all36_nano_lig.str
parameters              ../../toppar/toppar_all36_nano_lig_patch.str
parameters              ../../toppar/toppar_all36_synthetic_polymer.str
parameters              ../../toppar/toppar_all36_synthetic_polymer_patch.str
parameters              ../../toppar/toppar_all36_polymer_solvent.str
parameters              ../../toppar/toppar_water_ions.str
parameters              ../../toppar/toppar_dum_noble_gases.str
parameters              ../../toppar/toppar_ions_won.str
parameters              ../../toppar/cam.str
parameters              ../../toppar/toppar_all36_prot_arg0.str
parameters              ../../toppar/toppar_all36_prot_c36m_d_aminoacids.str
parameters              ../../toppar/toppar_all36_prot_fluoro_alkanes.str
parameters              ../../toppar/toppar_all36_prot_heme.str
parameters              ../../toppar/toppar_all36_prot_na_combined.str
parameters              ../../toppar/toppar_all36_prot_retinol.str
parameters              ../../toppar/toppar_all36_prot_model.str
parameters              ../../toppar/toppar_all36_prot_modify_res.str
parameters              ../../toppar/toppar_all36_na_nad_ppi.str
parameters              ../../toppar/toppar_all36_na_rna_modified.str
parameters              ../../toppar/toppar_all36_lipid_sphingo.str
parameters              ../../toppar/toppar_all36_lipid_archaeal.str
parameters              ../../toppar/toppar_all36_lipid_bacterial.str
parameters              ../../toppar/toppar_all36_lipid_cardiolipin.str
parameters              ../../toppar/toppar_all36_lipid_cholesterol.str
parameters              ../../toppar/toppar_all36_lipid_dag.str
parameters              ../../toppar/toppar_all36_lipid_inositol.str
parameters              ../../toppar/toppar_all36_lipid_lnp.str
parameters              ../../toppar/toppar_all36_lipid_lps.str
parameters              ../../toppar/toppar_all36_lipid_mycobacterial.str
parameters              ../../toppar/toppar_all36_lipid_miscellaneous.str
parameters              ../../toppar/toppar_all36_lipid_model.str
parameters              ../../toppar/toppar_all36_lipid_prot.str
parameters              ../../toppar/toppar_all36_lipid_tag.str
parameters              ../../toppar/toppar_all36_lipid_yeast.str
parameters              ../../toppar/toppar_all36_lipid_hmmm.str
parameters              ../../toppar/toppar_all36_lipid_detergent.str
parameters              ../../toppar/toppar_all36_lipid_ether.str
parameters              ../../toppar/toppar_all36_lipid_oxidized.str
parameters              ../../toppar/toppar_all36_carb_glycolipid.str
parameters              ../../toppar/toppar_all36_carb_glycopeptide.str
parameters              ../../toppar/toppar_all36_carb_imlab.str
parameters              ../../toppar/toppar_all36_label_spin.str
parameters              ../../toppar/toppar_all36_label_fluorophore.str
parameters              ../morph.prm


cellBasisVector1     $a   0.0   0.0
cellBasisVector2    0.0    $b   0.0
cellBasisVector3    0.0   0.0    $c
cellOrigin          0.0   0.0 $zcen
wrapAll                 on

langevin                on
langevinDamping         1.0
langevinTemp            $temp
langevinPiston          on
langevinPistonTarget    1.01325
langevinPistonPeriod    100.0
langevinPistonDecay     100.0
langevinPistonTemp      $temp

PME                     yes
PMEGridSpacing          1.0

exclude                 scaled1-4
1-4scaling              1.0
switching               on
cutoff                  12.0
switchdist              10.0
pairlistdist            14.0
stepspercycle           20

timestep                2.0
rigidBonds              all

alch                    on
singleTopology          on
alchType                FEP
alchFile                ligand.pdb
alchCol                 B
alchOutFreq             1000

outputEnergies          1000
outputPressure          1000
dcdfreq                 1000
xstFreq                 1000

binCoordinates          ligand_eq_forward.coor
binVelocities           ligand_eq_forward.vel
extendedSystem          ligand_eq_forward.xsc
EOF

# Base backward FEP config
cat > fep_base_backward.conf << 'EOF'
structure               ligand.psf
coordinates             ligand.pdb
source                  ligand.str

set temp                300

paraTypeCharmm          on;
parameters              ../../toppar/par_all36m_prot.prm
parameters              ../../toppar/par_all36_na.prm
parameters              ../../toppar/par_all36_carb.prm
parameters              ../../toppar/par_all36_lipid.prm
parameters              ../../toppar/par_all36_cgenff.prm
parameters              ../../toppar/par_interface.prm
parameters              ../../toppar/toppar_all36_moreions.str
parameters              ../../toppar/toppar_all36_nano_lig.str
parameters              ../../toppar/toppar_all36_nano_lig_patch.str
parameters              ../../toppar/toppar_all36_synthetic_polymer.str
parameters              ../../toppar/toppar_all36_synthetic_polymer_patch.str
parameters              ../../toppar/toppar_all36_polymer_solvent.str
parameters              ../../toppar/toppar_water_ions.str
parameters              ../../toppar/toppar_dum_noble_gases.str
parameters              ../../toppar/toppar_ions_won.str
parameters              ../../toppar/cam.str
parameters              ../../toppar/toppar_all36_prot_arg0.str
parameters              ../../toppar/toppar_all36_prot_c36m_d_aminoacids.str
parameters              ../../toppar/toppar_all36_prot_fluoro_alkanes.str
parameters              ../../toppar/toppar_all36_prot_heme.str
parameters              ../../toppar/toppar_all36_prot_na_combined.str
parameters              ../../toppar/toppar_all36_prot_retinol.str
parameters              ../../toppar/toppar_all36_prot_model.str
parameters              ../../toppar/toppar_all36_prot_modify_res.str
parameters              ../../toppar/toppar_all36_na_nad_ppi.str
parameters              ../../toppar/toppar_all36_na_rna_modified.str
parameters              ../../toppar/toppar_all36_lipid_sphingo.str
parameters              ../../toppar/toppar_all36_lipid_archaeal.str
parameters              ../../toppar/toppar_all36_lipid_bacterial.str
parameters              ../../toppar/toppar_all36_lipid_cardiolipin.str
parameters              ../../toppar/toppar_all36_lipid_cholesterol.str
parameters              ../../toppar/toppar_all36_lipid_dag.str
parameters              ../../toppar/toppar_all36_lipid_inositol.str
parameters              ../../toppar/toppar_all36_lipid_lnp.str
parameters              ../../toppar/toppar_all36_lipid_lps.str
parameters              ../../toppar/toppar_all36_lipid_mycobacterial.str
parameters              ../../toppar/toppar_all36_lipid_miscellaneous.str
parameters              ../../toppar/toppar_all36_lipid_model.str
parameters              ../../toppar/toppar_all36_lipid_prot.str
parameters              ../../toppar/toppar_all36_lipid_tag.str
parameters              ../../toppar/toppar_all36_lipid_yeast.str
parameters              ../../toppar/toppar_all36_lipid_hmmm.str
parameters              ../../toppar/toppar_all36_lipid_detergent.str
parameters              ../../toppar/toppar_all36_lipid_ether.str
parameters              ../../toppar/toppar_all36_lipid_oxidized.str
parameters              ../../toppar/toppar_all36_carb_glycolipid.str
parameters              ../../toppar/toppar_all36_carb_glycopeptide.str
parameters              ../../toppar/toppar_all36_carb_imlab.str
parameters              ../../toppar/toppar_all36_label_spin.str
parameters              ../../toppar/toppar_all36_label_fluorophore.str
parameters              ../morph.prm


cellBasisVector1     $a   0.0   0.0
cellBasisVector2    0.0    $b   0.0
cellBasisVector3    0.0   0.0    $c
cellOrigin          0.0   0.0 $zcen
wrapAll                 on

langevin                on
langevinDamping         1.0
langevinTemp            $temp
langevinPiston          on
langevinPistonTarget    1.01325
langevinPistonPeriod    100.0
langevinPistonDecay     100.0
langevinPistonTemp      $temp

PME                     yes
PMEGridSpacing          1.0

exclude                 scaled1-4
1-4scaling              1.0
switching               on
cutoff                  12.0
switchdist              10.0
pairlistdist            14.0
stepspercycle           20

timestep                2.0
rigidBonds              all

alch                    on
singleTopology          on
alchType                FEP
alchFile                ligand.pdb
alchCol                 B
alchOutFreq             1000

outputEnergies          1000
outputPressure          1000
dcdfreq                 1000
xstFreq                 1000

binCoordinates          ligand_eq_backward.coor
binVelocities           ligand_eq_backward.vel
extendedSystem          ligand_eq_backward.xsc
EOF

echo "Generating FEP window files..."

# Forward FEP windows (λ_i → λ_{i+1})
for i in {0..14}; do
    lambda1=${LAMBDA_VALUES[$i]}
    lambda2=${LAMBDA_VALUES[$((i+1))]}
    cat > fep_forward_${i}.conf << EOF
source fep_base_forward.conf
outputName              fep_forward_${i}
alchLambda              ${lambda1}
alchLambda2             ${lambda2}
run                     100000
EOF
done

# Backward FEP windows (λ_i → λ_{i-1})
for i in {15..1}; do
    lambda1=${LAMBDA_VALUES[$i]}
    lambda2=${LAMBDA_VALUES[$((i-1))]}
    cat > fep_backward_${i}.conf << EOF
source fep_base_backward.conf
outputName              fep_backward_${i}
alchLambda              ${lambda1}
alchLambda2             ${lambda2}
run                     100000
EOF
done

echo "Generated:"
echo "- 2 equilibration files (forward/backward)"
echo "- 2 base FEP configs"
echo "- 15 forward FEP windows"
echo "- 15 backward FEP windows"
