
# set up environment
cd ../../env/
. ./setEnvMarsGcc.sh > ../misc/fileForReadme/setEnvDjanGo.txt
cat ../misc/fileForReadme/setEnvDjanGo.txt

# build DjanGo
cd ${DJANGO_DIR}/build
make clean
make -j $DJANGO_NTHREADS | tee ${DJANGO_DIR}/misc/fileForReadme/make.txt

# build Tools
cd ${DJANGO_DIR}/misc/tool/build
make clean
make -j $DJANGO_NTHREADS | tee ${DJANGO_DIR}/misc/fileForReadme/makeTools.txt

# command line parameters
cd ${DJANGO_DIR}/misc/fileForReadme
${DJANGO_DIR}/bin/django -h | tee commandLineParam.txt

# version
${DJANGO_DIR}/bin/django -v | tee version.txt

# run validation tests
cd ${DJANGO_DIR}/validation/
sh runValidationTests.sh | tee ${DJANGO_DIR}/misc/fileForReadme/runValidationTests.txt

