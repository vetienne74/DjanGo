sh clean.sh

# build xml config files
echo === build xml config files ===
python3 build_config_xml.py

# to be used in the model 8km x 8km, without CPML
${DJANGO_DIR}/bin/build_acqui_file <<EOF
2
1 1
4500. 3500.
0. 0.
1  6
3500. 3500. 
0. 200.
EOF
mv acquisition.config.out.ascii acquisition2.config.out.ascii

${DJANGO_DIR}/bin/build_acqui_file <<EOF
2
1 1
1500. 500.
0. 0.
1  6
500. 500. 
0. 200.
EOF

${DJANGO_DIR}/bin/build_model_homo <<EOF
2
801 801
4000.
EOF

mv build_model_homo.out.bin vp.out.bin

${DJANGO_DIR}/bin/build_model_homo <<EOF
2
401 201
0.0
EOF

mv build_model_homo.out.bin vs.out.bin

${DJANGO_DIR}/bin/build_model_homo <<EOF
2
401 201
1.0
EOF

mv build_model_homo.out.bin rho.out.bin

# build rickerp.ascii for specfem2d
${DJANGO_DIR}/bin/wavelet <<EOF
2
10.0
0.001
1501
EOF
# src type (1=gauss,2=deriv gauss,3=ricker,4=deriv ricker, 5=monofreq)
# frequency
# dt
# nt
mv wavelet.out.bin   rickerp.out.bin
mv wavelet.out.ascii rickerp.out.ascii

# build ricker.bin for DjanGo 1st
${DJANGO_DIR}/bin/wavelet <<EOF
3
10.0
0.001
1501
EOF
# src type (1=gauss,2=deriv gauss,3=ricker,4=deriv ricker, 5=monofreq)
# frequency
# dt
# nt
mv wavelet.out.bin   ricker.out.bin
mv wavelet.out.ascii ricker.out.ascii

# build rickerd.bin for DjanGo 2nd
${DJANGO_DIR}/bin/wavelet <<EOF
4
10.0
0.001
1501
EOF
# src type (1=gauss,2=deriv gauss,3=ricker,4=deriv ricker, 5=monofreq)
# frequency
# dt
# nt
mv wavelet.out.bin   rickerd.out.bin
mv wavelet.out.ascii rickerd.out.ascii
