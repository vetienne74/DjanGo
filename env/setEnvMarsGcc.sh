
# Load needed modules

# MPI config
#export DJANGO_MPI_INVOKER='mpirun --oversubscribe'
export DJANGO_MPI_INVOKER='mpirun'

# OpenMP config
export DJANGO_NTHREADS=2
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$DJANGO_NTHREADS

# C++ compiler
export DJANGO_CPP=mpic++
export DJANGO_CPP_OPENACC_FLAG=
#export DJANGO_CPP_FLAGCOMP='-g -O3 -mavx2 -fopenmp'
export DJANGO_CPP_FLAGCOMP='-w -g -O3 -fopenmp -ffast-math -fpermissive -mavx2'
export DJANGO_CPP_LIB=

# Fortran compiler
export DJANGO_FC='mpif90'

# CUDA compiler
export DJANGO_CUDA=
export DJANGO_CUDA_FLAGCOMP=
export DJANGO_CUDA_LIB=

# HIP compiler
export DJANGO_HIP=
export DJANGO_HIP_FLAGCOMP=
export DJANGO_HIP_LIB=

# display Hpcscan settings
sh ./displayDjanGoEnv.sh

