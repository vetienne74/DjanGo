

nt=132301
ntrace=1

xwigb < pr.time.rec.django.out n1=$nt title=$1 &

exit

ximage < pr.time.snapshot.django.out n1=1002 &

exit

transp < pr.time.rec.django.out > pr.time.rec.django.out.tr n1=$ntrace perc
xwigb < pr.time.rec.django.out.tr n1=$nt title=$1 &

exit

transp < pr.freq.grid.django.out > pr.freq.grid.django.out.tr n1=2
ntmp=$((n1 * 2))
ximage <  pr.freq.grid.django.out.tr n1=$ntmp perc=99 &
