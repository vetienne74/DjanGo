
fig='MOD-HOMO-AC-INF-2D-'
n1=1501
file=django.config.fdm.o02.ac.1st.xml.energy
pswigb < $file > $file.ps n1=${n1} title=${file}
convert $file.ps $file.jpg
mv $file.jpg ${fig}${file}.jpg
display ${fig}${file}.jpg &

file=django.config.fdm.o02.ac.2nd.xml.energy
pswigb < $file > $file.ps n1=${n1} title=${file}
convert $file.ps $file.jpg
mv $file.jpg ${fig}${file}.jpg
display ${fig}${file}.jpg &

rm -f *.ps
