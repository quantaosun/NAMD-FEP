#!/bin/bash  
#Script to merge all forward/backward output files of .fepout

snum=1
enum=16

while [[ $snum -le $enum ]];  
do
cat md_forward_$snum.fepout >> complex_forward.fepout
((snum++))
done
