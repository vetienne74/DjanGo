
# TODO: fix huge performance issues with multiples threads
export OMP_NUM_THREADS=1
export GOMP_CPU_AFFINITY=0,1

sh clean.sh

# pre-processing
sh prepare.sh

# run DjanGo
# loop on all *.xml files present in the directory
for test_xml in $( ls django.config.*.xml); do

    echo '* RUN DJANGO...'
    mpirun -n 1 ${DJANGO_DIR}/bin/django -xml $test_xml
    
    #sh display.sh $test_xml
    cp pr.time.rec.django.out.bin ${test_xml}.pr.out.bin
    #cp src.wavelet.django.out ${test_xml}.src

done

# post-processing
#sh analytic_1d.sh


