
# build xml config files
echo === build xml config files ===
python3 build_config_xml.py

# build acquisition
${DJANGO_DIR}/bin/build_acqui_file <<EOF
2
1 1
0.5 0.5
0. 0.
101  101
0. 0. 
0.01 0.01
EOF

# compute analytical solution
${DJANGO_DIR}/bin/eigen_sol <<EOF
acquisition.config.out.txt
0.0
3.0
0.3535
EOF


