
report_file=$1

for file_name in $(ls *.out.bin); do

    ${DJANGO_DIR}/bin/check_file ${file_name} > tmp_check_file
    
    tmp_str=$(grep "FILE OK" tmp_check_file)
    if [ -z "$tmp_str" ] 
    then
	echo Check_file ${file_name} "FAILED" >> ${report_file}
	cat ./tmp_check_file
    else
	echo Check_file ${file_name} "PASSED" >> ${report_file}
    fi
    rm tmp_check_file
done


