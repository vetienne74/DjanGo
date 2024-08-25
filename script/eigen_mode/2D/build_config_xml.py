
xml_template = 'django.template.xml'

#------------------------------------------------------------------------------
# FDM - build xml config files from template
#------------------------------------------------------------------------------

scheme=['fdm']
SCHt='staggered'
equation=['acoustic','acoustic_lossy']
equation_order=['1','2']
space_order=['2']
NZ="51"
NX="51"
DZ="0.02"
DX="0.02"
DZRAND="0.0"
DXRAND="0.0"

for SCHm in scheme:
    for ACCs in space_order:
        for EQt in equation:
            for EQo in equation_order:                          

                # create new file name
                new_name = 'config.py'+'.'+SCHm+'.'+SCHt+'.'+'O'+ACCs+'.'+EQt+'.'+'O'+EQo
                xml_new = xml_template.replace('template', new_name)                

                # skip some invalid config                                    
                if EQt == 'string' and EQo == '1':
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

                    new_line = new_line.replace('_NZ_', NZ)
                    new_line = new_line.replace('_NX_', NX)
                    new_line = new_line.replace('_DZ_', DZ)
                    new_line = new_line.replace('_DX_', DX)
                    new_line = new_line.replace('_DZRAND_', DZRAND)
                    new_line = new_line.replace('_DXRAND_', DXRAND)

                    # ricker
                    if EQo == '1':
                        new_line = new_line.replace('_RICKER_', 'ricker')
                    if EQo == '2':
                        new_line = new_line.replace('_RICKER_', 'rickerd')
                    
                    f2.write(new_line)

                # close files
                f1.close()
                f2.close()
             
#------------------------------------------------------------------------------
# FEM - build xml config files from template
#------------------------------------------------------------------------------

NZ="21"
NX="21"
DZ="0.05"
DX="0.05"

scheme=['fem']
fem_type=['discontinuous','continuous']
distrib=['uniform']
equation=['acoustic','acoustic_lossy']
equation_order=['1','2']
space_order=['1']
node_type=['GLL']
DZRAND="0.0"
DXRAND="0.0"

for SCHm in scheme:
    for ACCs in space_order:
        for EQt in equation:
            for EQo in equation_order:                          
                for DIST in distrib:
                    for NODE in node_type:
                        for SCHt in fem_type:

                            # skip some invalid config                              
                            if SCHt == 'continuous' and DIST == 'random':
                                continue                                                       
                            if SCHt == 'discontinuous' and EQo == '2':
                                continue
                            if EQt == 'acoustic_lossy' and EQo == '2':
                                continue
                            if EQt == 'acoustic_lossy' and  SCHt == 'continuous':
                                continue
                            
                            # create new file name                           
                            new_name = 'config.py'+'.'+SCHm+'.'+SCHt+'.'+'O'+ACCs+'.'+EQt+'.'+DIST+'.'+NODE+'.'+'O'+EQo+\
                            '.z'+DZ+'.x'+DX+'.zran'+DZRAND+'.xran'+DXRAND
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
                                new_line = new_line.replace('_DIST_', DIST)
                                new_line = new_line.replace('_NODE_', NODE)

                                new_line = new_line.replace('_NZ_', NZ)
                                new_line = new_line.replace('_NX_', NX)
                                new_line = new_line.replace('_DZ_', DZ)
                                new_line = new_line.replace('_DX_', DX) 
                                new_line = new_line.replace('_DZRAND_', DZRAND)
                                new_line = new_line.replace('_DXRAND_', DXRAND)
                    
                                # ricker
                                if EQo == '1':
                                    new_line = new_line.replace('_RICKER_', 'ricker')
                                if EQo == '2':
                                    new_line = new_line.replace('_RICKER_', 'rickerd')

                                # quadrature
                                if NODE == 'equidistant':
                                    new_line = new_line.replace('_QUAD_', 'GL')
                                if NODE == 'GLL':
                                    new_line = new_line.replace('_QUAD_', 'GLL')
                    
                                f2.write(new_line)

                            # close files
                            f1.close()
                            f2.close()

exit()
                            
#------------------------------------------------------------------------------
# FEM - build xml config files from template
#------------------------------------------------------------------------------

NZ="21"
NX="11"
DZ="0.05"
DX="0.1"

scheme=['fem']
fem_type=['discontinuous','continuous']
distrib=['uniform']
equation=['acoustic']
equation_order=['1','2']
space_order=['4']
node_type=['GLL']
DZRAND2=['0.0','0.2']
DXRAND="0.0"

for SCHm in scheme:
    for ACCs in space_order:
        for EQt in equation:
            for EQo in equation_order:                          
                for DIST in distrib:
                    for NODE in node_type:
                        for SCHt in fem_type:
                            for DZRAND in DZRAND2:
                                # create new file name                           
                                new_name = 'config.py'+'.'+SCHm+'.'+SCHt+'.'+'O'+ACCs+'.'+EQt+'.'+DIST+'.'+NODE+'.'+'O'+EQo+\
                                '.z'+DZ+'.x'+DX+'.zran'+DZRAND+'.xran'+DXRAND
                                xml_new = xml_template.replace('template', new_name) 

                                # skip some invalid config                              
                                if SCHt == 'continuous' and DIST == 'random':
                                    continue
                                # skip some invalid config                              
                                if SCHt == 'discontinuous' and EQo == '2':
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
                                    new_line = new_line.replace('_DIST_', DIST)
                                    new_line = new_line.replace('_NODE_', NODE)

                                    new_line = new_line.replace('_NZ_', NZ)
                                    new_line = new_line.replace('_NX_', NX)
                                    new_line = new_line.replace('_DZ_', DZ)
                                    new_line = new_line.replace('_DX_', DX) 
                                    new_line = new_line.replace('_DZRAND_', DZRAND)
                                    new_line = new_line.replace('_DXRAND_', DXRAND)

                                    # ricker
                                    if EQo == '1':
                                        new_line = new_line.replace('_RICKER_', 'ricker')
                                    if EQo == '2':
                                        new_line = new_line.replace('_RICKER_', 'rickerd')

                                    # quadrature
                                    if NODE == 'equidistant':
                                        new_line = new_line.replace('_QUAD_', 'GL')
                                    if NODE == 'GLL':
                                        new_line = new_line.replace('_QUAD_', 'GLL')

                                    f2.write(new_line)

                                # close files
                                f1.close()
                                f2.close()
                            
