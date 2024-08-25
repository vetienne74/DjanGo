
nt=801

${DJANGO_DIR}/bin/analytic_layer <<EOF
acqui_cpml.config.out.ascii
4000.0
3000.0
10000.0
10.0
0.005
$nt
EOF

exit

file=pr.time.rec.analytic_layer.out
xwigb < $file n1=$nt title=$file &
