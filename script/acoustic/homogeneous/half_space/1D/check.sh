
failed="####### FAILED #######"
passed="------- PASSED -------"

if [ -z "$rms_max" ] 
then
    rms_max=0.1
fi

nrec=11

echo $(pwd)
grep Equation tmp
grep Numerical tmp
    
# check pr

file_1=pr.time.rec.django.out.bin
file_2=pr.time.rec.analytic_homo.out.bin

#${CWPROOT}/bin/transp < $file_1 > ${file_1}.tr n1=$nrec
${DJANGO_DIR}/bin/transpose_2d $file_1 $nrec
mv transpose_2d.out ${file_1}.tr

${DJANGO_DIR}/bin/check_rms >> tmp_check <<EOF
${file_1}.tr
$file_2
1.0
0
$rms_max
EOF

tmp_check_str=$(grep "RMS BELOW LIMIT" tmp_check)
if [ -z "$tmp_check_str" ] 
then
    echo $(pwd)
    echo $1
    echo $file_1 validation against analytic_homo $failed
    cat tmp_check
else
    echo $file_1 validation against analytic_homo $passed
    echo $(grep -B1 rms: tmp_check)
fi

rm tmp_check
rm ${file_1}.tr
