
# TODO: fix huge performance issues with multiples threads
export OMP_NUM_THREADS=1
export GOMP_CPU_AFFINITY=0,1

# run DjanGo
# loop on all *.xml files present in the directory
for test_xml in $( ls django.config.run*.xml); do

    #sh clean.sh
    
    echo '* RUN DJANGO...'
    mpirun -n 1 ${DJANGO_DIR}/bin/django -xml $test_xml
    
    #sh display.sh $test_xml
    cp pr.time.rec.django.out.bin ${test_xml}.pr.out.bin
    rm $test_xml

    # clean DjanGo output files
    rm -f *django.out*


done

# post-processing
#sh analytic_1d.sh


