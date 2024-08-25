
n1=101
perc=99

data_file=$1
ximage < ${data_file} n1=${n1} bclip=1.3 wclip=-1.3 title=${data_file} &

exit

data_file=vx.time.rec.django.out.bin
ximage < ${data_file} n1=${n1} perc=${perc} title=${data_file} &

data_file=vz.time.rec.django.out.bin
ximage < ${data_file} n1=${n1} perc=${perc} title=${data_file} &

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
