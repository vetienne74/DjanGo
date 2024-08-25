
n1=499
perc=99

data_file=pr.time.snapshot.django.out.bin
ximage < ${data_file} n1=${n1} perc=${perc} title='Pressure' &

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
