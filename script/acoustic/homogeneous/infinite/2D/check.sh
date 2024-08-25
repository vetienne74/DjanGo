
failed="####### FAILED #######"
passed="------- PASSED -------"

if [ -z "$rms_max" ] 
then
    rms_max=0.1
fi

nrec=6

echo $(pwd)
grep Equation tmp
grep Numerical tmp

# check vx traces
#----------------

file_1=vx.time.rec.django.out.bin
file_s=./specfem2d/vx.time.rec.specfem2d.out
# convert to double if needed
echo DOUBLE ${double_precision}
if [ "$double_precision" -eq 0 ]
then
    file_2=${file_s}
else
    ${DJANGO_DIR}/bin/single2double ${file_s}	
    file_2=./single2double.out.bin
fi
    
if [ -f "$file_1" ]
then
    ${DJANGO_DIR}/bin/transpose_2d $file_1 $nrec
    mv transpose_2d.out ${file_1}.tr
    
    ${DJANGO_DIR}/bin/check_rms >> tmp_check <<EOF
${file_1}.tr
$file_2
-1.0
0
$rms_max
EOF

    tmp_check_str=$(grep "RMS BELOW LIMIT" tmp_check)
    if [ -z "$tmp_check_str" ] 
    then
	echo $(pwd)
	echo $1
	echo $file_1 validation against specfem2d $failed
	cat tmp_check
    else
	echo $file_1 validation against specfem2d $passed
	echo $(grep -B1 rms: tmp_check)
    fi
    
    rm tmp_check
    rm ${file_1}.tr
    
fi

# check vz traces
#----------------

file_1=vz.time.rec.django.out.bin
file_s=./specfem2d/vz.time.rec.specfem2d.out
# convert to double if needed    
if [ "$double_precision" -eq 0 ]
then
    file_2=${file_s}
else
    ${DJANGO_DIR}/bin/single2double ${file_s}	
    file_2=./single2double.out.bin
fi
    
if [ -f "$file_1" ]
then
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
	echo $file_1 validation against specfem2d $failed
	cat tmp_check
    else
	echo $file_1 validation against specfem2d $passed
	echo $(grep -B1 rms: tmp_check)
    fi
    
    rm tmp_check
    rm ${file_1}.tr

fi

# check pr traces
#----------------

file_1=pr.time.rec.django.out.bin
file_s=./specfem2d/pr.time.rec.specfem2d.out
# convert to double if needed
if [ "$double_precision" -eq 0 ]
then
    file_2=${file_s}
else
    ${DJANGO_DIR}/bin/single2double ${file_s}	
    file_2=./single2double.out.bin
fi
    
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
    echo $file_1 validation against specfem2d $failed
    cat tmp_check
else
    echo $file_1 validation against specfem2d $passed
    echo $(grep -B1 rms: tmp_check)
fi

rm tmp_check
rm ${file_1}.tr

# check pr snapshot
#------------------

file_1=pr.time.snapshot.django.out.bin
file_s=REF.fdm.staggered.O4.acoustic.O1.pml.xml.snapshot.pr
# convert to double if needed
echo DOUBLE ${double_precision}
if [ "$double_precision" = "0" ]
then
    file_2=${file_s}
else
    ${DJANGO_DIR}/bin/single2double ${file_s}	
    file_2=./single2double.out.bin
fi

${DJANGO_DIR}/bin/check_rms >> tmp_check <<EOF
$file_1
$file_2
1.0
0
$rms_max
EOF

tmp_check_str=$(grep "RMS BELOW LIMIT" tmp_check)
if [ -z "$tmp_check_str" ] 
then
    echo $file_1 validation against reference snapshot $failed
    cat tmp_check
else
    echo $file_1 validation against reference snapshot $passed
    echo $(grep -B1 rms: tmp_check)
fi

rm tmp_check
