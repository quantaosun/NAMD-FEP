#!/bin/bash  
#Script to start running backward  

read -p "Enter forward starting number: " snum  
read -p "Enter forward ending number: " enum  

while [[ $snum -le $enum ]];  
do
/home/aistudio/NAMD_3.0alpha8_Linux-x86_64-multicore-CUDA/namd3 +p24 +setcpuaffinity --CUDASOAintegrate on +devices 0 md_backward_$snum.namd > md_backward_$snum.log  
echo "md_backward_$snum.namd is running."  
((snum++))  
done  

