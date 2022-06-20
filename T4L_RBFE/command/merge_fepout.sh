#!/bin/bash  
#Script to merge output files of .fepout
#including complex_forward.fepout, complex_backward.fepout, solvent_forward.fepout, and solvent_backward.fepout, total 4 files.

snum=1
enum=16

while [[ $snum -le $enum ]];  
do
cat md_forward_$snum.fepout >> complex_forward.fepout #here should be changed for different .fepout
((snum++))
done
