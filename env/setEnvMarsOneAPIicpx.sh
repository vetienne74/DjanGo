
# Load needed modules
. /opt/intel/oneapi/setvars.sh --force

# MPI config
export DJANGO_MPI_INVOKER='mpirun'

# OpenMP config
export DJANGO_NTHREADS=4
export KMP_AFFINITY=granularity=fine,compact
export OMP_NUM_THREADS=$DJANGO_NTHREADS

# C++ compiler
export DJANGO_CPP='mpiicpc -cxx=icpx'
export DJANGO_CPP_OPENACC_FLAG=
export DJANGO_CPP_FLAGCOMP='-w -g -O3 -fopenmp -xHost -ffast-math'
export DJANGO_CPP_LIB=

# Fortran compiler
export DJANGO_FC='mpiifort'

# CUDA compiler
export DJANGO_CUDA=
export DJANGO_CUDA_FLAGCOMP=
export DJANGO_CUDA_LIB=

# HIP compiler
export DJANGO_HIP=
export DJANGO_HIP_FLAGCOMP=
export DJANGO_HIP_LIB=

# display settings
sh ./displayDjanGoEnv.sh

