
nt=801
ntrace=11

#file_in=pr.time.rec.django.out.bin
file_in=$1
transp < ${file_in} > ${file_in}.tr n1=$ntrace
xwigb < ${file_in}.tr n1=$nt title=$1 &

exit

transp < pr.freq.grid.cornflex.out > pr.freq.grid.cornflex.out.tr n1=2
ntmp=$((n1 * 2))
ximage <  pr.freq.grid.cornflex.out.tr n1=$ntmp perc=99 &
