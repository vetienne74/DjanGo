
# clean bin dir
cd ../bin
rm -f *
cd -

# clean DjanGo objects
cd ../build
make clean
cd -

# clean tools objects

cd ../misc/build
make clean
cd -

# end of script
