##################################################################
#
#                    DJANGO VALIDATION TESTS
#
##################################################################

export rms_max=0.11
export test_1d=0
export test_2d=0
export test_3d=0
export appli_1d=0

export report_file=$1

#-----------------------------------------------------------------
# test DjanGo version
#-----------------------------------------------------------------

test_name=Test_Django_version
echo
echo $test_name ...

mpirun -n 1 $DJANGO_DIR/bin/django -version > tmp

tmp_str=$(grep "D j a n G o" tmp)
if [ -z "$tmp_str" ]
then
    cat tmp
    echo $test_name FAILED >> $report_file
else
    echo $test_name PASSED >> $report_file
    cat tmp
fi

#-----------------------------------------------------------------
# test 1D
#-----------------------------------------------------------------

test_dir=$DJANGO_DIR/validation/modelling/eigen_mode/1D_case
test_name=$test_dir
test_opt=-dryrun
if [ "$test_1d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

test_name=$test_dir
test_opt=
if [ "$test_1d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# End of validation tests
#-----------------------------------------------------------------

rm tmp

# END OF FILE
