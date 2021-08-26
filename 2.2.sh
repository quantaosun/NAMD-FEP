module load intel/18.0.1.163 openmpi/3.1.2-intel namd/2.13-mpi

export Honey_NAMD=/apps/namd/2.13-mpi/arch/Linux-x86_64-icc   # this need to be modified

###################################################################
mpirun --oversubscribe -np 4 $Fuck_NAMD/namd2 +ppn 7 nvt_equil_modi.namd > nvt_equil_modi.out &&

#NPT###############################################################

mpirun --oversubscribe -np 4 $Fuck_NAMD/namd2 +ppn 7 npt_equil_modi.namd > npt_equil_modi.out &&

# Production#######################################################

mpirun --oversubscribe -np 4 $Fuck_NAMD/namd2 +ppn 7 md_forward_modi.namd > md_forward_modi.out &&

cd ../../solvent/solvent-input-files &&

mpirun --oversubscribe -np 4 $Fuck_NAMD/namd2 +ppn 7 nvt_equil_modi.namd > nvt_equil_modi.out &&

#NPT###############################################################

mpirun --oversubscribe -np 4 $Fuck_NAMD/namd2 +ppn 7 npt_equil_modi.namd > npt_equil_modi.out &&

# Production#######################################################

mpirun --oversubscribe -np 4 $Fuck_NAMD/namd2 +ppn 7 md_forward_modi.namd > md_forward_modi.out 
