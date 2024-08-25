ntrace=6
nt=1501

data_file=$1
transp < ${data_file} >  ${data_file}.tr n1=${ntrace}
xwigb < ${data_file}.tr n1=${nt} perc=${perc} title=${data_file} &


