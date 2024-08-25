
xml_template = 'django.template.xml'

#------------------------------------------------------------------------------
# FDM - build xml config files from template
#------------------------------------------------------------------------------

bound=['pml']
scheme=['fdm']
SCHt='staggered'
equation=['acoustic','acoustic_lossy']
equation_order=['1','2']
space_order=['2']
NZ="401"
NX="201"
DZ="5.0"
DX="10.0"
ACQUI="acquisition.config.out.ascii"
ACQUI2="acquisition2.config.out.ascii"
src_type="gaussian"

for boundary in bound:
    for SCHm in scheme:
        for ACCs in space_order:
            for EQt in equation:
                for EQo in equation_order:                          

                    # skip config
                    if EQo == '2' and boundary == 'sponge':
                        continue
                    if EQt == 'elastic' and EQo == '2':                                      
                        continue
                    
                    # create new file name
                    new_name = 'config.py'+'.'+SCHm+'.'+SCHt+'.'+'O'+ACCs+'.'+EQt+'.'+'O'+EQo+'.'+boundary
                    xml_new = xml_template.replace('template', new_name)

                    # open files
                    f2 = open(xml_new, 'w')
                    f1 = open(xml_template, 'r')

                    # create new file
                    for line in f1:
                        new_line = line.replace('_SCHm_', SCHm)
                        new_line = new_line.replace('_SCHt_', SCHt)
                        new_line = new_line.replace('_ACCs_', ACCs)
                        new_line = new_line.replace('_EQt_', EQt)   
                        new_line = new_line.replace('_EQo_', EQo)

                        # ricker
                        if EQo == '1':
                            new_line = new_line.replace('_RICKER_', 'ricker')
                        if EQo == '2':
                            new_line = new_line.replace('_RICKER_', 'rickerd')

                        new_line = new_line.replace('_NZ_', NZ)
                        new_line = new_line.replace('_NX_', NX)
                        new_line = new_line.replace('_DZ_', DZ)
                        new_line = new_line.replace('_DX_', DX)
                        new_line = new_line.replace('_ACQUI_', ACQUI)
                        new_line = new_line.replace('_SRC_TYPE_', src_type)

                        # boundary
                        if (boundary == 'pml'):
                            new_line=new_line.replace('_BOUND1_', 'boundary edge="zbeg" type="cpml" width="10" coef="1.e-12"')
                            new_line=new_line.replace('_BOUND2_', 'boundary edge="zend" type="cpml" width="10" coef="1.e-12"')
                            new_line=new_line.replace('_BOUND3_', 'boundary edge="xbeg" type="cpml" width="10" coef="1.e-12"')
                            new_line=new_line.replace('_BOUND4_', 'boundary edge="xend" type="cpml" width="10" coef="1.e-12"') 
                        elif (boundary == 'sponge'):
                            new_line=new_line.replace('_BOUND1_', 'boundary edge="zbeg" type="sponge" width="45" coef="0.0053"')
                            new_line=new_line.replace('_BOUND2_', 'boundary edge="zend" type="sponge" width="45" coef="0.0053"')
                            new_line=new_line.replace('_BOUND3_', 'boundary edge="xbeg" type="sponge" width="45" coef="0.0053"')
                            new_line=new_line.replace('_BOUND4_', 'boundary edge="xend" type="sponge" width="45" coef="0.0053"')
                        elif (boundary == 'freesurf'):
                            new_line=new_line.replace('_BOUND1_', 'boundary edge="zbeg" type="freesurf"')
                            new_line=new_line.replace('_BOUND2_', 'boundary edge="zend" type="freesurf"')
                            new_line=new_line.replace('_BOUND3_', 'boundary edge="xbeg" type="freesurf"') 
                            new_line=new_line.replace('_BOUND4_', 'boundary edge="xend" type="freesurf"') 
                        
                        f2.write(new_line)

                    # close files
                    f1.close()
                    f2.close()


# fine grid for O2 is used to get accurate solution
space_order=['2']
NZ="801"
NX="801"
DZ="2.5"
DX="2.5"
test=False

if test == True:
    for boundary in bound:
        for SCHm in scheme:
            for ACCs in space_order:
                for EQt in equation:
                    for EQo in equation_order:                          

                        # skip config
                        if EQo == '2' and boundary == 'sponge':
                            continue
                        if EQt == 'elastic' and EQo == '2':                                      
                            continue

                        # create new file name
                        new_name = 'config.py'+'.'+SCHm+'.'+SCHt+'.'+'O'+ACCs+'.'+EQt+'.'+'O'+EQo+'.'+boundary
                        xml_new = xml_template.replace('template', new_name)                                  

                        # open files
                        f2 = open(xml_new, 'w')
                        f1 = open(xml_template, 'r')

                        # create new file
                        for line in f1:
                            new_line = line.replace('_SCHm_', SCHm)
                            new_line = new_line.replace('_SCHt_', SCHt)
                            new_line = new_line.replace('_ACCs_', ACCs)
                            new_line = new_line.replace('_EQt_', EQt)   
                            new_line = new_line.replace('_EQo_', EQo)

                            # ricker
                            if EQo == '1':
                                new_line = new_line.replace('_RICKER_', 'ricker')
                            if EQo == '2':
                                new_line = new_line.replace('_RICKER_', 'rickerd')

                            new_line = new_line.replace('_NZ_', NZ)
                            new_line = new_line.replace('_NX_', NX)
                            new_line = new_line.replace('_DZ_', DZ)
                            new_line = new_line.replace('_DX_', DX)
                            new_line = new_line.replace('_ACQUI_', ACQUI)
                            new_line = new_line.replace('_SRC_TYPE_', src_type)

                            # boundary
                            if (boundary == 'pml'):
                                new_line=new_line.replace('_BOUND1_', 'boundary edge="zbeg" type="cpml" width="10" coef="1.e-12"')
                                new_line=new_line.replace('_BOUND2_', 'boundary edge="zend" type="cpml" width="10" coef="1.e-12"')
                                new_line=new_line.replace('_BOUND3_', 'boundary edge="xbeg" type="cpml" width="10" coef="1.e-12"')
                                new_line=new_line.replace('_BOUND4_', 'boundary edge="xend" type="cpml" width="10" coef="1.e-12"') 
                            elif (boundary == 'sponge'):
                                new_line=new_line.replace('_BOUND1_', 'boundary edge="zbeg" type="sponge" width="45" coef="0.0053"')
                                new_line=new_line.replace('_BOUND2_', 'boundary edge="zend" type="sponge" width="45" coef="0.0053"')
                                new_line=new_line.replace('_BOUND3_', 'boundary edge="xbeg" type="sponge" width="45" coef="0.0053"')
                                new_line=new_line.replace('_BOUND4_', 'boundary edge="xend" type="sponge" width="45" coef="0.0053"')  
                            elif (boundary == 'freesurf'):
                                new_line=new_line.replace('_BOUND1_', 'boundary edge="zbeg" type="freesurf"')
                                new_line=new_line.replace('_BOUND2_', 'boundary edge="zend" type="freesurf"')
                                new_line=new_line.replace('_BOUND3_', 'boundary edge="xbeg" type="freesurf"') 
                                new_line=new_line.replace('_BOUND4_', 'boundary edge="xend" type="freesurf"')

                            f2.write(new_line)

                        # close files
                        f1.close()
                        f2.close()
                    
#------------------------------------------------------------------------------
# FEM - build xml config files from template
#------------------------------------------------------------------------------

scheme=['fem']
fem_type=['discontinuous','continuous']
equation=['acoustic']
equation_order=['1']
space_order=['5']
NZ="21"
NX="21"
NZ2="81"
NX2="81"
DZ="100.0"
DX="100.0"
bound=['sponge']

for boundary in bound:
    for SCHt in fem_type:
        for SCHm in scheme:
            for ACCs in space_order:
                for EQt in equation:
                    for EQo in equation_order:                          

                        # create new file name
                        new_name = 'config.py'+'.'+SCHm+'.'+SCHt+'.'+'O'+ACCs+'.'+EQt+'.'+'O'+EQo+'.'+boundary
                        xml_new = xml_template.replace('template', new_name)

                        # skip config
                        if EQo == '2' and boundary == 'sponge':
                            continue
                        if EQt == 'elastic' and EQo == '2':                                      
                            continue

                        # open files
                        f2 = open(xml_new, 'w')
                        f1 = open(xml_template, 'r')

                        # create new file
                        for line in f1:
                            new_line = line.replace('_SCHm_', SCHm)
                            new_line = new_line.replace('_SCHt_', SCHt)
                            new_line = new_line.replace('_ACCs_', ACCs)
                            new_line = new_line.replace('_EQt_', EQt)   
                            new_line = new_line.replace('_EQo_', EQo)

                            # ricker
                            if EQo == '1':
                                new_line = new_line.replace('_RICKER_', 'ricker')
                            if EQo == '2':
                                new_line = new_line.replace('_RICKER_', 'rickerd')

                            new_line = new_line.replace('_SRC_TYPE_', src_type)
                            new_line = new_line.replace('_DZ_', DZ)
                            new_line = new_line.replace('_DX_', DX)

                            # boundary                            
                            if (boundary == 'sponge'):
                                new_line = new_line.replace('_NZ_', NZ)
                                new_line = new_line.replace('_NX_', NX)
                                new_line = new_line.replace('_ACQUI_', ACQUI)
                                new_line=new_line.replace('_BOUND1_', 'boundary edge="zbeg" type="sponge" width="6" coef="0.94"')
                                new_line=new_line.replace('_BOUND2_', 'boundary edge="zend" type="sponge" width="6" coef="0.94"')
                                new_line=new_line.replace('_BOUND3_', 'boundary edge="xbeg" type="sponge" width="6" coef="0.94"')
                                new_line=new_line.replace('_BOUND4_', 'boundary edge="xend" type="sponge" width="6" coef="0.94"') 
                            elif (boundary == 'freesurf'):
                                new_line = new_line.replace('_NZ_', NZ2)
                                new_line = new_line.replace('_NX_', NX2)
                                new_line = new_line.replace('_ACQUI_', ACQUI2)
                                new_line=new_line.replace('_BOUND1_', 'boundary edge="zbeg" type="freesurf"')
                                new_line=new_line.replace('_BOUND2_', 'boundary edge="zend" type="freesurf"')
                                new_line=new_line.replace('_BOUND3_', 'boundary edge="xbeg" type="freesurf"') 
                                new_line=new_line.replace('_BOUND4_', 'boundary edge="xend" type="freesurf"')

                            f2.write(new_line)

                        # close files
                        f1.close()
                        f2.close()              

        
