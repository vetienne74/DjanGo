
sh clean.sh

# build xml config files
echo === build xml config files ===
python3 build_config_xml.py

# acquisition with cpml
${DJANGO_DIR}/bin/build_acqui_file <<EOF
1
1
2000.
0.
11
4500.
1000.
EOF

mv acquisition.config.out.ascii acqui_cpml.config.out.ascii

# acquisition without cpml
${DJANGO_DIR}/bin/build_acqui_file <<EOF
1
1
12000.
0.
11
14500.
1000.
EOF

mv acquisition.config.out.ascii acqui_no_cpml.config.out.ascii

# model with PML
# layer 1
${DJANGO_DIR}/bin/build_model_homo <<EOF
1
1000 
4000.
EOF
mv build_model_homo.out.bin vp1.out.bin

# layer 2
${DJANGO_DIR}/bin/build_model_homo <<EOF
1
1001
3000.
EOF
mv build_model_homo.out.bin vp2.out.bin
cat vp1.out.bin vp2.out.bin > vp.cpml.out.bin
rm vp1.out.bin vp2.out.bin

# model without PML
# layer 1
${DJANGO_DIR}/bin/build_model_homo <<EOF
1
2000 
4000.
EOF
mv build_model_homo.out.bin vp1.out.bin

# layer 2
${DJANGO_DIR}/bin/build_model_homo <<EOF
1
2001
3000.
EOF
mv build_model_homo.out.bin vp2.out.bin
cat vp1.out.bin vp2.out.bin > vp.no_cpml.out.bin
rm vp1.out.bin vp2.out.bin

# model without PML for FEM -> revert from FDM
# layer 1
${DJANGO_DIR}/bin/build_model_homo <<EOF
1
2000 
3000.
EOF
mv build_model_homo.out.bin vp1.out.bin

# layer 2
${DJANGO_DIR}/bin/build_model_homo <<EOF
1
2001
4000.
EOF
mv build_model_homo.out.bin vp2.out.bin
cat vp1.out.bin vp2.out.bin > vp.no_cpml2.out.bin
rm vp1.out.bin vp2.out.bin

# analytic solution
sh analytic_layer.sh
