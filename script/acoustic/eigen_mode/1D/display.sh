
nt=801
ntrace=101
file_in=pr.time.rec.django.out.bin
transp < ${file_in} > ${file_in}.tr n1=$ntrace perc
ximage < ${file_in}.tr n1=$nt bclip=0.96 wclip=-0.96 title=$1 &

exit

file_in=vz.time.rec.django.out.bin
transp < ${file_in} > ${file_in}.tr n1=$ntrace perc
ximage < ${file_in}.tr n1=$nt perc=99 title=$1 &

exit

transp < pr.freq.grid.cornflex.out > pr.freq.grid.cornflex.out.tr n1=2
ntmp=$((n1 * 2))
ximage <  pr.freq.grid.cornflex.out.tr n1=$ntmp perc=99 &
