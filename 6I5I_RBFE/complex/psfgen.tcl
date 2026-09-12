package require psfgen
topology /home/aistudio/work/NAMD-FEP/toppar/top_all36_prot.rtf
topology /home/aistudio/work/NAMD-FEP/6I5I_RBFE/complex/ligand.rtf
pdbalias atom ILE CD1 CD
pdbalias atom SER HG HG1
pdbalias residue HIS HSD

segment P1 { first NTER; last CTER; pdb /home/aistudio/work/NAMD-FEP/6I5I_RBFE/complex/protein_seg1.pdb }
coordpdb /home/aistudio/work/NAMD-FEP/6I5I_RBFE/complex/protein_seg1.pdb P1
segment P2 { first NTER; last CTER; pdb /home/aistudio/work/NAMD-FEP/6I5I_RBFE/complex/protein_seg2.pdb }
coordpdb /home/aistudio/work/NAMD-FEP/6I5I_RBFE/complex/protein_seg2.pdb P2
segment LIG { first none; last none; pdb /home/aistudio/work/NAMD-FEP/6I5I_RBFE/complex/ligand.pdb }
coordpdb /home/aistudio/work/NAMD-FEP/6I5I_RBFE/complex/ligand.pdb LIG

guesscoord
writepsf complex.psf
writepdb complex.pdb
exit
