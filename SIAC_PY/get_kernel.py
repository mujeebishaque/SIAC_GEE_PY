import ee, math

def makeBRDFKernels(image):

    BR = 1.0
    HB = 2.0
    d2r = math.pi / 180.0
    zthresh = 0.00001

    sza = image.select('sza').float().multiply(ee.Number(d2r));
    vza = image.select('vza').float().multiply(ee.Number(d2r));
    vaa = image.select('vaa').float().multiply(ee.Number(d2r));
    saa = image.select('saa').float().multiply(ee.Number(d2r));
    raa = vaa.subtract(saa)
    raa_plus = raa.add(ee.Number(math.pi))
    w = vza.lt(0)
    vza = ee.Image(vza.where(w,vza.multiply(ee.Number(-1)))).rename('vza')
    raa = ee.Image(raa.where(w,raa_plus))
    w = sza.lt(0)
    sza = ee.Image(sza.where(w,sza.multiply(ee.Number(-1)))).rename('sza')
    raa = ee.Image(raa.where(w,raa_plus))
    raa = ee.Image(0).expression('raa % (2*pi)',{'raa': raa,'pi':math.pi}).rename('raa')     

    cos_vza = vza.cos().rename('cos_vza')
    sin_vza = vza.sin().rename('sin_vza')
    cos_sza = sza.cos().rename('cos_sza')
    sin_sza = sza.sin().rename('sin_sza')
    cos_raa = raa.cos().rename('cos_raa')
    sin_raa = raa.sin().rename('sin_raa')
    tanti   = sza.tan().rename('tanti')
    tantv   = vza.tan().rename('tantv')

    tan1    = ee.Image(ee.Image(0).expression('BR*tan1',{'tan1': tanti,'BR':BR}))
    angp1 = tan1.atan()
    sin1 = angp1.sin()
    cos1 = angp1.cos()
    w = cos1.lte(zthresh)
    cos1 = cos1.where(w,zthresh) 

    tan2 = ee.Image(ee.Image(0).expression('BR*tan1',{'tan1': tantv,'BR':BR}))
    angp2 = tan2.atan()
    sin2 = angp2.sin()
    cos2 = angp2.cos()
    
    w = cos2.lte(zthresh)
    cos2 = cos2.where(w,zthresh)

    cdict = {'cos1': cos_vza,'sin1': sin_vza,'cos2': cos_sza,'sin2': sin_sza,'cos3':cos_raa}
    cosphaang = ee.Image(0).expression('cos1*cos2 + sin1*sin2*cos3',cdict)

    w = cosphaang.lte(-1)
    cosphaang = ee.Image(cosphaang.where(w,-1))
    w = cosphaang.gte(1)
    cosphaang = ee.Image(cosphaang.where(w,1)).rename('cos_phaang')
    phaang = cosphaang.acos().rename('phaang')
    sinphaang = phaang.sin().rename('sin_phaang')
    

    cdict = {'cosphaang': cosphaang,'sinphaang': sinphaang,'pi': math.pi, 'phaang':phaang,
            'cos1': cos_vza, 'cos2': cos_sza}
    ross = ee.Image(0).expression('((pi/2. - phaang)*cosphaang+sinphaang)/(cos1+cos2) - pi/4.',cdict).rename('ross')

    cdict = {'tan1': tan1,'tan2': tan2,'cos3':cos_raa}
    temp = ee.Image(0).expression('tan1*tan1 + tan2*tan2 - 2*tan1*tan2*cos3',cdict)
    w = temp.lte(0)
    temp = temp.where(w,0)
    distance = temp.sqrt()

    cdict = {'cos1': cos1,'sin1': sin1,'cos2': cos2,'sin2': sin2,'cos3':cos_raa}
    temp = ee.Image(0).expression('1./cos1 + 1./cos2',cdict)

    cdict = {'tan1': tan1,'tan2': tan2,'cos3':cos_raa,'HB':HB,'distance':distance,'sin3':sin_raa,'temp':temp}
    cost = ee.Image(0).expression('HB * sqrt(distance * distance + tan1 * tan1 * tan2 * tan2 * sin3 * sin3) / temp',cdict)
    w = cost.lte(-1)
    cost = cost.where(w,-1)
    w = cost.gte(1)
    cost = cost.where(w,1)
    tvar = cost.acos()
    sint = tvar.sin()
    
    cdict = {'tvar': tvar,'sint': sint,'cost':cost,'pi':math.pi, 'temp':temp}
    overlap = ee.Image(0).expression('(1/pi) * (tvar - sint * cost) * temp',cdict)
    w = overlap.lte(0)
    overlap = overlap.where(w,0).rename('overlap')

    cdict = {'cos1': cos1,'sin1': sin1,'cos2': cos2,'sin2': sin2,'cos3':cos_raa}
    cosphaang2 = ee.Image(0).expression('cos1*cos2 + sin1*sin2*cos3',cdict)
    cdict = {'overlap': overlap,'cosphaang2': cosphaang2,'cos1':cos1,'cos2':cos2, 'temp':temp}
    li = ee.Image(0).expression('overlap - temp + 0.5 * (1. + cosphaang2) / cos1 / cos2',cdict).rename('li')
    isotropic = ee.Image.constant(1.).rename('isotropic')


    simu_names = ['simu_boa_b1', 'simu_boa_b2', 'simu_boa_b3', 'simu_boa_b4', 'simu_boa_b5', 'simu_boa_b6']
    simu_boas = []
    for i in range(5, -1, -1):
        band = i + 1
        simu = (image.select(mcd43_names[i])
                .add(image.select(mcd43_names[i+6]).multiply(ross))
                .add(image.select(mcd43_names[i+12]).multiply(li))
                .rename(simu_names[i]))
        simu_boas.append(simu)
    return image.addBands([isotropic, ross, li, cos_sza, cos_vza, cos_raa]).addBands(simu_boas)
    

mcd43_names = ['BRDF_Albedo_Parameters_Band1_iso', 
                    'BRDF_Albedo_Parameters_Band2_iso', 
                    'BRDF_Albedo_Parameters_Band3_iso', 
                    'BRDF_Albedo_Parameters_Band4_iso', 
                    'BRDF_Albedo_Parameters_Band6_iso', 
                    'BRDF_Albedo_Parameters_Band7_iso',
                    'BRDF_Albedo_Parameters_Band1_vol', 
                    'BRDF_Albedo_Parameters_Band2_vol', 
                    'BRDF_Albedo_Parameters_Band3_vol', 
                    'BRDF_Albedo_Parameters_Band4_vol',
                    'BRDF_Albedo_Parameters_Band6_vol', 
                    'BRDF_Albedo_Parameters_Band7_vol', 
                    'BRDF_Albedo_Parameters_Band1_geo', 
                    'BRDF_Albedo_Parameters_Band2_geo', 
                    'BRDF_Albedo_Parameters_Band3_geo', 
                    'BRDF_Albedo_Parameters_Band4_geo', 
                    'BRDF_Albedo_Parameters_Band6_geo', 
                    'BRDF_Albedo_Parameters_Band7_geo']