import ee, math
from SIAC_PY import s2_cloud as s2_cloud
from SIAC_PY import Sen2Cloud as Sen2Cloud
from SIAC_PY import Inverse_6S_AOT_brdf as inverse_6S_AOT
from SIAC_PY import Inverse_S2_TCWV as inverse_S2_TCWV
from SIAC_PY import get_kernel as kernel
from SIAC_PY import TWO_NN as tnn
from SIAC_PY import interp_MCD43A1 as interp_mcd43a1
from SIAC_PY import S2B_AC as s2b_ac

# s2_cloud = require('users/marcyinfeng/SIAC/:/S2_cloud')
# Sen2Cloud = require('users/marcyinfeng/SIAC/:/Sen2Cloud')
# inverse_6S_AOT = require('users/marcyinfeng/SIAC/:/Inverse_6S_AOT_brdf')
# inverse_S2_TCWV = require('users/marcyinfeng/SIAC/:/Inverse_S2_TCWV')
# kernel = require('users/marcyinfeng/SIAC/:/get_kernel')
# interp_mcd43a1 = require('users/marcyinfeng/SIAC/:/interp_MCD43A1')
# tnn = require('users/marcyinfeng/SIAC/:/Two_NN')

# s2b_ac = require('users/marcyinfeng/SIAC/:/S2B_AC')


mcd43_names  = kernel.mcd43_names
get_cloud    = s2_cloud.get_cloud
common_bands = s2_cloud.common_bands
s2_bands     = s2_cloud.s2_bands
l8_bands     = s2_cloud.l8_bands


def get_boa(image):
    image     = ee.Image(image)
    geom  = image.geometry()
    image_date = image.date()
    projection = image.select('B2').projection()
    crs = projection.crs()
    image2 = ee.Image(image.divide(10000))
    orig_image = ee.Image(image.select(s2_bands[0:6], common_bands[0:6]).divide(10000).copyProperties(image)).set('system:time_start', image.get('system:time_start'))

    cloud_bands = ['B1', 'B2', 'B4', 'B5', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
    control = image.select(cloud_bands).multiply(0.0001).reproject(crs, None, 10)
    cloud = tnn.predict(control, Sen2Cloud.Sen2Cloud_H1_scale,
                                    Sen2Cloud.Sen2Cloud_H1_offset,
                                    Sen2Cloud.Sen2Cloud_H2_scale,
                                    Sen2Cloud.Sen2Cloud_H2_offset,
                                    Sen2Cloud.Sen2Cloud_Out_scale,
                                    Sen2Cloud.Sen2Cloud_Out_offset).rename('cloud_prob')
    cloud_prob = ee.Image(1).divide(ee.Image(1).add((ee.Image(-1).multiply(cloud)).exp()))
    cloud = cloud_prob.gt(0.4).rename('cloud')
    dia_cloud  = cloud.focal_max(kernel= ee.Kernel.circle(radius=60, units='meters'))
    image = ee.Image(image.select(s2_bands[0:6], common_bands[0:6]).divide(10000).updateMask(dia_cloud.eq(0)).copyProperties(image)).set('system:time_start', image.get('system:time_start'))

    mcd43a1    = interp_mcd43a1.interp_mcd43a1(image)

    psf = ee.Kernel.gaussian(radius=1500, sigma= 300, units= 'meters', normalize= True)

    image = image.convolve(psf)
    saa  = ee.Image.constant(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE'))).rename('saa')
    sza  = ee.Image.constant(ee.Number(image.get('MEAN_SOLAR_ZENITH_ANGLE' ))).rename('sza')
    vaa  = ee.Image.constant(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B2'))).rename('vaa')
    vza  = ee.Image.constant(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B2' ))).rename('vza')

    image = image.addBands([saa, sza, vaa, vza]).addBands(mcd43a1.select(mcd43_names[0:18]))
    image = kernel.makeBRDFKernels(image)
    cos_raa = image.select('cos_raa')
    cos_sza = image.select('cos_sza')
    cos_vza = image.select('cos_vza')
    tco3 = ee.ImageCollection('TOMS/MERGED').filterDate(ee.Date(image_date).advance(-3, 'day'), ee.Date(image_date).advance(4, 'day')).mean().select('ozone').multiply(0.001).clip(geom)
    median = tco3.reduceRegion(ee.Reducer.percentile([50]), geom, 3000, None, None, False, 10e13).get('ozone')
    tco3 = tco3.unmask().where(tco3.mask().eq(0), ee.Number(median))


    elevation = ee.Image('USGS/SRTMGL1_003').select('elevation').multiply(0.001)
    image3 = ee.Image.cat([(image.select('BLUE').subtract(1.381298790609708504e-02)).divide(1.006170958348520772e+00), 
                                (image.select('GREEN').subtract(0.00190909687075)).divide(1.00087659541449), 
                                mcd43a1.select('BRDF_Albedo_Parameters_Band3_iso'),
                                mcd43a1.select('BRDF_Albedo_Parameters_Band3_vol'),
                                mcd43a1.select('BRDF_Albedo_Parameters_Band3_geo'),
                                mcd43a1.select('BRDF_Albedo_Parameters_Band4_iso'),
                                mcd43a1.select('BRDF_Albedo_Parameters_Band4_vol'),
                                mcd43a1.select('BRDF_Albedo_Parameters_Band4_geo'),
                                cos_sza, cos_vza, cos_raa, tco3.rename('tco3'), elevation.rename('ele')]).reproject(crs, None, 3000).float().clip(geom).float()
                        
    boa_diff  = image.select('simu_boa_b5').subtract(image.select('SWIR1')).rename('diff')
    pp = boa_diff.reduceRegion(ee.Reducer.percentile([15, 50, 85]), geom, 3000, None, None, False, 10e13)
    diff_min  = ee.Number(pp.get('diff_p15'))
    diff_max  = ee.Number(pp.get('diff_p85'))
    stab_mask = (boa_diff.lte(diff_max)).And(boa_diff.gte(diff_min))

    boa_diff  = image.select('simu_boa_b6').subtract(image.select('SWIR2')).rename('diff')
    pp = boa_diff.reduceRegion(ee.Reducer.percentile([15, 50, 85]), geom, 3000, None, None, False, 10e13)
    diff_min  = ee.Number(pp.get('diff_p15'))
    diff_max  = ee.Number(pp.get('diff_p85'))
    stab_mask = (boa_diff.lte(diff_max)).And(boa_diff.gte(diff_min)).And(stab_mask)

    median = ee.Number(pp.get('diff_p50'))
    diff_min  = ee.Number(median).subtract(0.02)
    diff_max  = ee.Number(median).add(0.02)
    stab_mask = (boa_diff.lte(diff_max)).And(boa_diff.gte(diff_min)).And(boa_diff.abs().lt(0.03))

    BLUE = image3.select('BLUE')
    pp   = BLUE.reduceRegion(ee.Reducer.percentile([10, 90]), geom, 3000, None, None, False, 10e13)
    BLUE_min  = ee.Number(pp.get('BLUE_p10'))
    BLUE_max  = ee.Number(pp.get('BLUE_p90'))
    BLUE_mask = (image3.select('BLUE').lte(BLUE_max)).And(image3.select('BLUE').gte(BLUE_min))
    SWIR2 = image.select('SWIR2')
    pp = SWIR2.reduceRegion(ee.Reducer.percentile([10, 90]), geom, 3000, None, None, False, 10e13)
    SWIR_min  = ee.Number(pp.get('SWIR2_p10'))
    SWIR_max  = ee.Number(pp.get('SWIR2_p90'))
    SWIR_mask = (SWIR2.lte(SWIR_max)).And(SWIR2.gte(SWIR_min))
    image3 = image3.updateMask(BLUE_mask.And(stab_mask).And(image3.mask()).And(SWIR_mask))
    aot = tnn.predict(image3, inverse_6S_AOT.Inverse_6S_AOT_brdf_H1_scale,
                                inverse_6S_AOT.Inverse_6S_AOT_brdf_H1_offset,
                                inverse_6S_AOT.Inverse_6S_AOT_brdf_H2_scale,
                                inverse_6S_AOT.Inverse_6S_AOT_brdf_H2_offset,
                                inverse_6S_AOT.Inverse_6S_AOT_brdf_Out_scale,
                                inverse_6S_AOT.Inverse_6S_AOT_brdf_Out_offset).rename('aot')

    valid_pix = aot.mask().reduceRegion(reducer= ee.Reducer.sum(), geometry= geom, scale= 3000, maxPixels=10e13)
    prior_aot = ee.ImageCollection('MODIS/006/MOD08_M3').select('Aerosol_Optical_Depth_Land_Ocean_Mean_Mean').filter(ee.Filter.date(ee.Date(image_date).advance(-1, 'year'), image_date)).first().multiply(0.001).clip(geom).rename('aot')
    aot = ee.Image(ee.Algorithms.If(ee.Number(valid_pix.get('aot')).lt(5), prior_aot, aot))

    pp = aot.reduceRegion(ee.Reducer.percentile([10, 50, 90]), geom, 3000, None, None, False, 10e13)
    aot_min  = ee.Number(pp.get('aot_p10')).subtract(0.001)
    aot_max  = ee.Number(pp.get('aot_p90')).add(0.001)
    median   = ee.Number(pp.get('aot_p50'))
    aot_mask = (aot.lte(aot_max)).And(aot.gte(aot_min))
    aot = aot.updateMask(aot_mask.And(aot.mask()))
    filled_aot = aot.unmask().where(aot.mask().eq(0), ee.Number(median))
    kernelSize = 10
    ker = ee.Kernel.square(kernelSize * 500, "meters")
    aot     = filled_aot.reduceNeighborhood(reducer=ee.Reducer.mean(), kernel= ker, inputWeight='mask', skipMasked=False, optimization='boxcar').rename('aot')
    aot = aot.max(0.01)
    control = ee.Image.cat([ image2.select('B9'), 
                            image2.select('B8A'), 
                            cos_sza, cos_vza, cos_raa, aot, tco3, elevation]).updateMask(dia_cloud.eq(0).And(image2.select('B8A').gt(0.12))).reproject(crs, None, 3000).float().float()
    tcwv = tnn.predict(control,inverse_S2_TCWV.Inverse_S2_TCWV_H1_scale,
                            inverse_S2_TCWV.Inverse_S2_TCWV_H1_offset,
                            inverse_S2_TCWV.Inverse_S2_TCWV_H2_scale,
                            inverse_S2_TCWV.Inverse_S2_TCWV_H2_offset,
                            inverse_S2_TCWV.Inverse_S2_TCWV_Out_scale,
                            inverse_S2_TCWV.Inverse_S2_TCWV_Out_offset).rename('tcwv')

    pp = tcwv.reduceRegion(ee.Reducer.percentile([10, 50, 90]), geom, 3000, None, None, False, 10e13)
    tcwv_min  = ee.Number(pp.get('tcwv_p10'))
    tcwv_max  = ee.Number(pp.get('tcwv_p90'))
    median    = ee.Number(pp.get('tcwv_p50'))
    tcwv_mask = (tcwv.lte(tcwv_max)).And(tcwv.gte(tcwv_min))
    tcwv = tcwv.updateMask(tcwv_mask.And(tcwv.mask()))
    filled_tcwv = tcwv.unmask().where(tcwv.mask().eq(0), ee.Number(median))
    kernelSize = 6
    ker = ee.Kernel.square(kernelSize * 500, "meters")
    tcwv     = filled_tcwv.reduceNeighborhood(reducer= ee.Reducer.mean(), kernel= ker, inputWeight= 'mask', skipMasked=False, optimization='boxcar').rename('tcwv')
    tcwv = tcwv.max(0.01)
    image = orig_image
    saa  = ee.Image.constant(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE'))).rename('saa')
    sza  = ee.Image.constant(ee.Number(image.get('MEAN_SOLAR_ZENITH_ANGLE' ))).rename('sza')
    vaa  = ee.Image.constant(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B2'))).rename('vaa')
    vza  = ee.Image.constant(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B2' ))).rename('vza')
    raa  = vaa.subtract(saa)
    deg2rad = ee.Number(math.pi).divide(ee.Number(180.0))
    cos_sza = (sza.multiply(deg2rad)).cos()
    cos_vza = (vza.multiply(deg2rad)).cos()
    cos_raa = (raa.multiply(deg2rad)).cos()
    ele = ee.Image('USGS/SRTMGL1_003').select('elevation').multiply(0.001)
    res = 3000
    image4 = ee.Image.cat([cos_sza, cos_vza, cos_raa, aot, tcwv, tco3, ele]).reproject(crs, None, res).float().clip(geom)

    name     = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A',  'B11', 'B12']
    new_name = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A',  'B11', 'B12']

    boa = s2b_ac.S2B_AC(image4, image2.select(name, new_name))
    return ee.Image(boa.select(new_name, name)).addBands([aot.rename('AOT'), tcwv.rename('WVP'), cloud.rename('MSK_CLDPRB'),tco3.rename('TCO3')])
