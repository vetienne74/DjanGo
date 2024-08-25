
ntrace=6
nt=1501

data_file=pr.time.rec.django.out.bin
transp < ${data_file} >  ${data_file}.tr n1=${ntrace}
xwigb < ${data_file}.tr n1=${nt} perc=${perc} title='Vx' &

exit

n1=225
perc=99

data_file=pr.time.snapshot.django.out.bin
ximage < ${data_file} n1=${n1} perc=${perc} title=$1 &

exit

data_file=vx.time.rec.cornflex.out
transp < ${data_file} >  ${data_file}.tr n1=${ntrace}
ximage < ${data_file}.tr n1=${nt} perc=${perc} title='Vx' &

data_file=vz.time.rec.cornflex.out
transp < ${data_file} >  ${data_file}.tr n1=${ntrace}
ximage < ${data_file}.tr n1=${nt} perc=${perc} title='Vz' &

exit

transp < pr.freq.grid.cornflex.out > pr.freq.grid.cornflex.out.tr n1=2
ximage <  pr.freq.grid.cornflex.out.tr n1=$n1 perc=99 &
