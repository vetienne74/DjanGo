
nt=801

${DJANGO_DIR}/bin/analytic_homo <<EOF
acqui_cpml.config.out.txt
4000.0
10.0
0.005
$nt
1
EOF

exit

file=pr.time.rec.analytic_homo.out
xwigb < $file n1=$nt title=$file &
