
nt=1501
nrec=6
perc=99

data_file=$1
transp < ${data_file} > ${data_file}.tr n1=${nrec}
xwigb < ${data_file}.tr n1=${nt} perc=${perc} title=$1 &

