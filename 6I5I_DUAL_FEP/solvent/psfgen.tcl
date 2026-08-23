package require psfgen
topology /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP/hybrid/hybrid.rtf
segment LIG { first none; last none; pdb /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP/hybrid/hybrid.pdb }
coordpdb /home/aistudio/work/NAMD-FEP/6I5I_DUAL_FEP/hybrid/hybrid.pdb LIG
guesscoord
writepsf ligand.psf
writepdb ligand.pdb
exit
