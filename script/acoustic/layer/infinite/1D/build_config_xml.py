
xml_template = 'django.template.xml'

#------------------------------------------------------------------------------
# build files for FDM
#------------------------------------------------------------------------------

method='fdm'
type_method='staggered'
h='10'
dt='0.0005'
for boundary in ['pml']:
    for weq in ['1st','2nd']:
        for order in ['12']:
            for src in ['src_func']:

                # create new file name
                new_name = 'config.py'+'.'+method+'.O'+order+'.'+weq+'.h'+h+'.'+boundary+'.'+src
                xml_new = xml_template.replace('template', new_name)
                
                f1 = open(xml_template, 'r')               
                f2 = open(xml_new, 'w')
                
                for line in f1:
                    new_line=line.replace('_ORDER_', order)
                    new_line=new_line.replace('_H_', h)
                    new_line=new_line.replace('_METHOD_', method)
                    new_line=new_line.replace('_TYPE_', type_method)   
                    new_line=new_line.replace('_DT_', dt)
                    
                    # wave eq. order
                    if (weq == '1st'):
                        new_line=new_line.replace('_WEQ_', '1')
                    elif (weq == '2nd'):
                        new_line=new_line.replace('_WEQ_', '2')      
                    
                    # source functioon
                    if (weq == '1st'):
                        if (src == 'src_func'):
                            new_line=new_line.replace('_SRC_', 'source type="point" sigma="5.0" subtype="explosive" function="ricker" frequency="10.0"')
                        elif (src == 'src_file1'):
                            new_line=new_line.replace('_SRC_', 'source type="point" sigma="5.0" subtype="explosive" function="file" name="ricker_dt_0_0001s.out.bin" dt="0.0001"')
                        elif (src == 'src_file2'):
                            new_line=new_line.replace('_SRC_', 'source type="point" sigma="5.0" subtype="explosive" function="file" name="ricker_dt_0_00015s.out.bin" dt="0.00015"')
                    elif (weq == '2nd'):
                        if (src == 'src_func'):
                            new_line=new_line.replace('_SRC_', 'source type="point" sigma="5.0" subtype="explosive" function="rickerd" frequency="10.0"')
                        elif (src == 'src_file1'):
                            new_line=new_line.replace('_SRC_', 'source type="point" sigma="5.0" subtype="explosive" function="file" name="rickerd_dt_0_0001s.out.bin" dt="0.0001"')
                        elif (src == 'src_file2'):
                            new_line=new_line.replace('_SRC_', 'source type="point" sigma="5.0" subtype="explosive" function="file" name="rickerd_dt_0_00015s.out.bin" dt="0.00015"')
                            
                    # boundary
                    if (boundary == 'pml'):
                        new_line=new_line.replace('_BOUND1_', 'boundary edge="zbeg" type="cpml" width="10" coef="1.e-12"')
                        new_line=new_line.replace('_BOUND2_', 'boundary edge="zend" type="cpml" width="10" coef="1.e-12"')
                        new_line=new_line.replace('_ACQUI_', 'acqui_cpml.config.out.ascii')
                        new_line=new_line.replace('_NZ_', '2001')
                        new_line=new_line.replace('_VP_FILE_', 'vp.cpml.out.bin')
                    else:
                        new_line=new_line.replace('_ACQUI_', 'acqui_no_cpml.config.out.ascii')
                        new_line=new_line.replace('_NZ_', '4001')
                        new_line=new_line.replace('_VP_FILE_', 'vp.no_cpml.out.bin')
                    f2.write(new_line)
f1.close()
f2.close()

#------------------------------------------------------------------------------
# build files for FEM
#------------------------------------------------------------------------------

method='fem'
type_method=['discontinuous']
h='10'
dt='0.0005'
for type_m in type_method:
    for boundary in ['sponge']:
        for weq in ['1st']:
            for order in ['4']:
                for src in ['src_func']:
                    f1 = open(xml_template, 'r')
                    new_name = 'config.py'+'.'+method+'.'+type_m+'.O'+order+'.'+weq+'.h'+h+'.'+boundary+'.'+src
                    xml_new = xml_template.replace('template', new_name)
                    f2 = open(xml_new, 'w')
                    for line in f1:
                        new_line=line.replace('_ORDER_', order)
                        new_line=new_line.replace('_H_', h)
                        new_line=new_line.replace('_METHOD_', method)
                        new_line=new_line.replace('_TYPE_', type_m)
                        new_line=new_line.replace('_DT_', dt)

                        # wave eq. order
                        if (weq == '1st'):
                            new_line=new_line.replace('_WEQ_', '1')
                        elif (weq == '2nd'):
                            new_line=new_line.replace('_WEQ_', '2')      

                        # source function
                        if (weq == '1st'):
                            if (src == 'src_func'):
                                new_line=new_line.replace('_SRC_', 'source type="gaussian" sigma="1.0" subtype="explosive" function="ricker" frequency="10.0"')
                            elif (src == 'src_file1'):
                                new_line=new_line.replace('_SRC_', 'source type="gaussian" sigma="1.0" subtype="explosive" function="file" name="ricker_dt_0_0001s.out.bin" dt="0.0001"')
                            elif (src == 'src_file2'):
                                new_line=new_line.replace('_SRC_', 'source type="gaussian" sigma="1.0" subtype="explosive" function="file" name="ricker_dt_0_00015s.out.bin" dt="0.00015"')
                        elif (weq == '2nd'):
                            if (src == 'src_func'):
                                new_line=new_line.replace('_SRC_', 'source type="gaussian" sigma="1.0" subtype="explosive" function="rickerd" frequency="10.0"')
                            elif (src == 'src_file1'):
                                new_line=new_line.replace('_SRC_', 'source type="gaussian" sigma="1.0" subtype="explosive" function="file" name="rickerd_dt_0_0001s.out.bin" dt="0.0001"')
                            elif (src == 'src_file2'):
                                new_line=new_line.replace('_SRC_', 'source type="gaussian" sigma="1.0" subtype="explosive" function="file" name="rickerd_dt_0_00015s.out.bin" dt="0.00015"')

                        # boundary
                        if (boundary == 'sponge'):
                            new_line=new_line.replace('_BOUND1_', 'boundary edge="zbeg" type="sponge" width="10" coef="0.92"')
                            new_line=new_line.replace('_BOUND2_', 'boundary edge="zend" type="sponge" width="10" coef="0.92"')                            
                            new_line=new_line.replace('_ACQUI_', 'acqui_cpml.config.out.ascii')
                            new_line=new_line.replace('_NZ_', '2001')
                            new_line=new_line.replace('_VP_FILE_', 'vp.cpml.out.bin')
                        else:
                            new_line=new_line.replace('_BOUND1_', 'boundary edge="zbeg" type="freesurf"')
                            new_line=new_line.replace('_BOUND2_', 'boundary edge="zend" type="freesurf"')
                            new_line=new_line.replace('_ACQUI_', 'acqui_no_cpml.config.out.ascii')
                            new_line=new_line.replace('_NZ_', '4001')
                            new_line=new_line.replace('_VP_FILE_', 'vp.no_cpml.out.bin')
                        f2.write(new_line)
f1.close()
f2.close()

