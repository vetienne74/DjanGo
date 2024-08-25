
nt=801
ntrace=11

#file_in=pr.time.rec.django.out.bin
file_in=$1
transp < ${file_in} > ${file_in}.tr n1=$ntrace
xwigb < ${file_in}.tr n1=$nt title=$1 &

exit
ximage < pr.time.snapshot.django.out.bin n1=2521 &
