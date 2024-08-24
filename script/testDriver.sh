##################################################################
#
#                    DJANGO VALIDATION TESTS
#
##################################################################

export rms_max=0.10
export test_1d=1
export test_2d=0
export test_3d=0
export appli_1d=0

export report_file=$1

#-----------------------------------------------------------------
# test DjanGo version
#-----------------------------------------------------------------

test_name=Django_version
echo | tee $report_file
echo Running $test_name ... | tee $report_file

mpirun -n 1 $DJANGO_DIR/bin/django -version > tmp

tmp_str=$(grep "D j a n G o" tmp)
if [ -z "$tmp_str" ]
then
    echo "> FAILED" | tee -a $report_file
else
    echo "> PASSED" | tee -a $report_file
fi
#cat tmp
rm tmp

#-----------------------------------------------------------------
# test eigen mode 1D dryrun
#-----------------------------------------------------------------

test_name=eigen_mode/1D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=-dryrun
if [ "$test_1d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#test_name=$test_dir
#test_opt=
#if [ "$test_1d" -eq 0 ]
#then
#    echo SKIP $test_name
#else
#    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
#fi

#-----------------------------------------------------------------
# End of validation tests
#-----------------------------------------------------------------

# END OF FILE