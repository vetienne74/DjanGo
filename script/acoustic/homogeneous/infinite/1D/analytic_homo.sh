
nt=801

${DJANGO_DIR}/bin/analytic_homo <<EOF
acqui_no_cpml.config.out.ascii
4000.0
10.0
0.005
$nt
0
EOF

exit

file=pr.time.rec.analytic_homo.out
xwigb < $file n1=$nt title=$file &
