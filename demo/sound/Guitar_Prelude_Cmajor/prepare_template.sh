
# acquisition
${DJANGO_DIR}/bin/build_acqui_file <<EOF
1
1
0.50
0.
1
0.55
0.
EOF

# model
${DJANGO_DIR}/bin/build_model_homo <<EOF
1
651
_VP_
EOF

mv build_model_homo.out.bin vp.out.bin
