import subprocess

prepare_template = 'prepare_tmp.sh'
xml_template     = 'django.config.template.xml'

t0=0
t60=20.0
loss1=18.42/t60
loss2=0.1

# new param
loss1=2.e-5
loss2=4.e-8
#loss2=0.0

# slow
#delta=0.75
# fast
delta=0.15

# E7#9 notes:
# E4 E5 G#5 D6 G6 E6
string=1
for freq in [329.63, 659.25, 830.61, 1174.66, 1567.98, 1318.51]:    
   
    xml_new = xml_template.replace('template', 'run'+str(string))
    f1 = open(xml_template, 'r')       
    f2 = open(xml_new, 'w')
    for line in f1:
        new_line=line.replace('_t0_', str(t0))
        new_line=new_line.replace('_VP_', str(freq))
        new_line=new_line.replace('_LOSS1_', str(loss1))
        new_line=new_line.replace('_LOSS2_', str(loss2))
        f2.write(new_line)
    f1.close()
    f2.close()
    
    t0=t0+delta
    string=string+1   
    subprocess.call(['sh','./run2.sh'])
