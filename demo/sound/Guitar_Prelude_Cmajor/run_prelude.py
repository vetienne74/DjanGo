import subprocess

prepare_template = 'prepare_template.sh'
xml_template     = 'django.config.template.xml'
string=1

# string length (m)
length=0.65

# loss

# guitar
t60=6.0

# organ
#t60=25.0

# string 1 = high E
# string 6 = low E

# one oactave higher 
# E3 A3 D4 G4 B5 E5
# for freq in [659.25, 493.88, 391.99, 293.66, 220.0, 164.81]:

# standard tuning
# E2 A2 D3 G3 B4 E4
#for freq in [329.63, 246.94, 196.00, 146.83, 110.0, 82.41]:    

for freq in [659.25, 493.88, 391.99, 293.66, 220.0, 164.81]:
#for freq in [164.81]:
    f1 = open(prepare_template, 'r')       
    f2 = open('prepare.sh', 'w')
    vp=length*freq
    print('*** frequency ', freq)
    print('*** vp        ', vp)
    for line in f1:
        new_line=line.replace('_VP_', str(vp))                             
        f2.write(new_line)
    f1.close()
    f2.close()

    xml_new = xml_template.replace('template', 'run'+str(string))
    f1 = open(xml_template, 'r')       
    f2 = open(xml_new, 'w')
    for line in f1:
        new_line=line.replace('_num_', str(string))
        new_line=new_line.replace('_LOSS_', str(t60))   
        f2.write(new_line)
    f1.close()
    f2.close()
    
    string=string+1
    #subprocess.call(['cat','./prepare.sh'])
    subprocess.call(['sh','./prepare.sh'])
    subprocess.call(['sh','./run2.sh'])
