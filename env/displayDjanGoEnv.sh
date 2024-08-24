
echo DjanGo is set on machine `hostname` with:
echo - C++ Compiler
$DJANGO_CPP --version

echo - FORTRAN Compiler
$DJANGO_FC --version

if [ -f $DJANGO_CUDA ]
then
    echo - No CUDA compiler set
else
    echo - CUDA Compiler
    $DJANGO_CUDA --version
fi
if [ -f $DJANGO_HIP ]
then
    echo - No HIP compiler set
else
    echo - HIP Compiler
    $DJANGO_HIP --version
fi
echo - MPI library
$DJANGO_MPI_INVOKER --version
echo - Number of OpenMP threads $DJANGO_NTHREADS

echo You are ready to go!

