#!/bin/bash
#PBS -N ratecal 
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l walltime=96:00:00
#PBS -q longn
#PBS -joe
#PBS -V


cd $PBS_O_WORKDIR
export I_MPI_FABRICS shm:dapl
export I_MPI_MPD_TMPDIR /scratch/$USER

echo Working directory is $PBS_O_WORKDIR

./rateliq.x >&  out  
