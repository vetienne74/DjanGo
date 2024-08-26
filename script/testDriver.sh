##################################################################
#
#                    DJANGO VALIDATION TESTS
#
##################################################################

export rms_max=0.10
export test_1d=1
export test_2d=1
export test_3d=0
export appli_1d=0

export report_file=$1

#-----------------------------------------------------------------
# test DjanGo executable without xml input file
#-----------------------------------------------------------------

test_name=Django_executable
echo | tee $report_file
echo Running $test_name ... | tee $report_file

mpirun -n 1 $DJANGO_DIR/bin/django -version > tmp

tmp_str=$(grep "D j a n G o" tmp)
if [ -z "$tmp_str" ]
then
    echo "FAILED" >> $report_file
else
    echo "PASSED" >> $report_file
fi
#cat tmp
rm tmp

#-----------------------------------------------------------------
# test acoustic eigen mode 1D [dryrun]
#-----------------------------------------------------------------

test_name=acoustic/eigen_mode/1D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=-dryrun
if [ "$test_1d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# test acoustic eigen mode 1D
#-----------------------------------------------------------------

test_name=acoustic/eigen_mode/1D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=
if [ "$test_1d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# test acoustic eigen mode 2D
#-----------------------------------------------------------------

test_name=acoustic/eigen_mode/2D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=
if [ "$test_2d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# test acoustic homogeneous infinite 1D
#-----------------------------------------------------------------

test_name=acoustic/homogeneous/infinite/1D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=
if [ "$test_1d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# test acoustic homogeneous infinite 2D
#-----------------------------------------------------------------

test_name=acoustic/homogeneous/infinite/2D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=
if [ "$test_2d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# test acoustic homogeneous half-space 1D
#-----------------------------------------------------------------

test_name=acoustic/homogeneous/half_space/1D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=
if [ "$test_1d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# test acoustic homogeneous half-space 2D
#-----------------------------------------------------------------

test_name=acoustic/homogeneous/half_space/2D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=
if [ "$test_2d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# test acoustic layer infinite 1D
#-----------------------------------------------------------------

test_name=acoustic/layer/infinite/1D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=
if [ "$test_1d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# test acoustic layer infinite 2D
#-----------------------------------------------------------------

test_name=acoustic/layer/infinite/2D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=
if [ "$test_2d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# test elastic homogeneous infinite 2D
#-----------------------------------------------------------------

test_name=elastic/homogeneous/infinite/2D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=
if [ "$test_2d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# test elastic homogeneous half-space 2D
#-----------------------------------------------------------------

test_name=elastic/homogeneous/half_space/2D
test_dir=$DJANGO_DIR/script/$test_name
test_opt=
if [ "$test_2d" -eq 0 ]
then
    echo SKIP $test_name
else
    sh ./runDjango.sh $report_file $test_name $test_dir $test_opt
fi

#-----------------------------------------------------------------
# End of validation tests
#-----------------------------------------------------------------

# END OF FILE
