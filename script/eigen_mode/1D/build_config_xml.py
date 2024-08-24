
xml_template = 'django.template.xml'

#------------------------------------------------------------------------------
# FDM - build xml config files from template
#------------------------------------------------------------------------------

scheme=['fdm']
SCHt='staggered'
equation=['acoustic','acoustic_lossy']
equation_order=['1','2']
space_order=['2']
dzrandom="0.0"

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
                    new_line = new_line.replace('_DZRANDOM_', dzrandom)     

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

scheme=['fem']
fem_type=['discontinuous','continuous']
distrib=['uniform']
equation=['acoustic','acoustic_lossy']
equation_order=['1','2']
space_order=['1']
node_type=['GLL']
dzrandom=['0.0']

for SCHm in scheme:
    for ACCs in space_order:
        for EQt in equation:
            for EQo in equation_order:                          
                for DIST in distrib:
                    for NODE in node_type:
                        for SCHt in fem_type:
                            for DZR in dzrandom:
                                # create new file name                           
                                new_name = 'config.py'+'.'+SCHm+'.'+SCHt+'.'+'O'+ACCs+'.'+EQt+'.'+DIST+'.'+NODE+'.'+'O'+EQo+'.rand'+DZR
                                xml_new = xml_template.replace('template', new_name) 

                                # skip some invalid config                              
                                if SCHt == 'continuous' and NODE == 'equidistant':
                                    continue
                                if SCHt == 'continuous' and EQt == 'acoustic_lossy':
                                    continue
                                if EQo == '2' and SCHt == 'discontinuous':
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
                                    new_line = new_line.replace('_DZRANDOM_', DZR)

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

