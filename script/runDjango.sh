
export report_file=$1
export test_name=$2
export test_dir=$3
export test_opt=$4

failed="FAILED"
passed="PASSED"
warning="WARNING"

echo ""
echo "Running $test_name" $test_opt "..." | tee -a ${report_file} 

if [ -d $test_dir ]
then
    cd $test_dir
else
    echo $failed - $test_dir does not exist | tee -a ${report_file} 
    exit
fi

sh prepare.sh > tmp

# loop on all *.xml files present in the directory
for test_xml in $( ls django.config*.xml); do

    echo "* Running DjanGo with $test_xml ..." | tee -a ${report_file}
    mpirun -n 1 ${DJANGO_DIR}/bin/django -xml $test_xml $test_opt > tmp

    # check job completed
    test_case=Job_completed
    tmp_str=$(grep "DJANGO TERMINATED SUCCESSFULLY" tmp)
    if [ -z "$tmp_str" ] 
    then
	# JOB FAILED
	echo $test_case $failed >> ${report_file}
	cat tmp
    else
	# JOB COMPLETED
	echo $test_case $passed >> ${report_file}
	#echo "FILE: " $(pwd)/$test_xml
	#echo $(grep "TOTAL TIME" tmp)
	#echo $(grep "Speed (GFlop/s)" tmp)
	#echo $(grep "Speed (Gcell/s)" tmp)
	tmp_str2=$(grep "TOTAL TIME (sec)" tmp)
	if [ -z "$tmp_str2" ] 
	then
	    echo RUN TIME EXCEED 1 MINUTE $warning >> ${report_file}
	fi
    
	# check there is no nan or inf in output
	test_case=Nan_or_inf
	tmp_str=$(grep 'nan\|inf' tmp)
	if [ -z "$tmp_str" ] 
	then
	    echo $test_case $passed >> ${report_file} >> ${report_file}
	else
	    echo $test_xml $test_case $failed >> ${report_file}
	    echo $tmp_str >> ${report_file}
	fi
    
	# check memory has been freed
	test_case=Check_memory
	tmp_str=$(grep "MEMORY CLEAN" tmp)
	if [ -z "$tmp_str" ] 
	then
	    echo $test_xml $test_case $failed >> ${report_file}
	    grep "BYTES STILL" tmp >> ${report_file}
	else
	    echo $test_case $passed >> ${report_file}
	fi
	
	if [ "$test_opt" != "-dryrun" ]
	then
	    # check ouput binary file are Ok (no NaN, nor Inf)
	    sh ${DJANGO_DIR}/script/checkFiles.sh ${report_file}

	    # check numerical solution	
	    if [ -f "check.sh" ]
	    then
		sh check.sh $test_xml >> ${report_file}
	    fi

	    # clean DjanGo output files
	    rm -f *django.out.bin
	fi
    fi

done

# clean test dir
sh ${DJANGO_DIR}/script/clean.sh

rm tmp
