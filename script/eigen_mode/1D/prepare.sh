
# build xml config files
echo === build xml config files ===
python3 build_config_xml.py

# build acquisition
${DJANGO_DIR}/bin/build_acqui_file <<EOF
1
1
0.5
0.
101
0.
0.01
EOF

#mv build_model_homo.out.bin vp.bin

# compute analytical solution
${DJANGO_DIR}/bin/eigen_sol <<EOF
acquisition.config.out.txt
0.0
4.0
0.005
EOF

exit

nt=801
file_in=pr.time.rec.eigen_sol.out.bin
ximage < ${file_in} n1=$nt perc=99 title=${file_in} &


