
if [ -f $DJANGO_CPP ]
then
    echo It seems DjanGo environment has not been set.
    echo Source one of the env files in ../env or create a new one for your machine.
    exit
fi

#double_precision=1
double_precision=0
tester=`whoami`
machine=`hostname`

if [ "${double_precision}" -eq 0 ]
then
    reportFile=runValidationTests.${machine}.${tester}.out
else
    reportFile=runValidationTests.${machine}.${tester}.double.out
fi
reportFileWithPath=${DJANGO_DIR}/script/$reportFile
rm -f $reportFile

echo "========================================================================" 
echo "                         START VALIDATION TESTS                         " 
echo "========================================================================" 

start_time=$(date)
sh testDriver.sh $reportFileWithPath 
end_time=$(date)

echo "" 
echo "# Started  : " $start_time 
echo "# Ended    : " $end_time 
echo "On machine : " $machine 
echo "Done by    : " $tester 

echo "" 
echo '* SUMMARY' 
echo '# PASSED  :' $(grep PASSED $reportFileWithPath | wc -l) 
echo '# FAILED  :' $(grep FAILED $reportFileWithPath | wc -l) 
echo '# ERROR   :' $(grep 'E R R O R' $reportFileWithPath | wc -l) 
echo '# WARNING :' $(grep 'W A R N I N G' $reportFileWithPath | wc -l) 
echo "==> Results of tests are saved in ./"$reportFile 

echo "" 
echo "========================================================================" 
echo "                         END VALIDATION TESTS                           " 
echo "========================================================================" 

# END OF FILE

