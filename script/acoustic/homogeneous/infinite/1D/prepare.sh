
sh clean.sh

# build xml config files
echo === build xml config files ===
python3 build_config_xml.py

# acquisition with cpml
${DJANGO_DIR}/bin/build_acqui_file <<EOF
1
1
10000.
0.
11
0.
2000.
EOF

mv acquisition.config.out.ascii acqui_cpml.config.out.ascii

# acquisition without cpml
${DJANGO_DIR}/bin/build_acqui_file <<EOF
1
1
20000.
0.
11
10000.
2000.
EOF

mv acquisition.config.out.ascii acqui_no_cpml.config.out.ascii

# analytic solution
sh analytic_homo.sh

# src wavelet 
${DJANGO_DIR}/bin/wavelet <<EOF
3
10.0
0.0001
40001
EOF
# src type (1=gauss,2=deriv gauss,3=ricker,4=deriv ricker, 5=monofreq)
# frequency
# dt
# nt
mv wavelet.out.bin ricker_dt_0_0001s.out.bin

${DJANGO_DIR}/bin/wavelet <<EOF
3
10.0
0.00015
25000
EOF
# src type (1=gauss,2=deriv gauss,3=ricker,4=deriv ricker, 5=monofreq)
# frequency
# dt
# nt
mv wavelet.out.bin ricker_dt_0_00015s.out.bin

${DJANGO_DIR}/bin/wavelet <<EOF
4
10.0
0.0001
40001
EOF
# src type (1=gauss,2=deriv gauss,3=ricker,4=deriv ricker, 5=monofreq)
# frequency
# dt
# nt
mv wavelet.out.bin rickerd_dt_0_0001s.out.bin

${DJANGO_DIR}/bin/wavelet <<EOF
4
10.0
0.00015
25000
EOF
# src type (1=gauss,2=deriv gauss,3=ricker,4=deriv ricker, 5=monofreq)
# frequency
# dt
# nt
mv wavelet.out.bin rickerd_dt_0_00015s.out.bin
