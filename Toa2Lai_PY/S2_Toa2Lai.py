import ee, math, gc
from Toa2Lai_PY import NN_cloud as Sen2Cloud
from Toa2Lai_PY import mcd19_prior as mcd19_prior
from Toa2Lai_PY import NN_prosail as NN_prosail
from Toa2Lai_PY import get_xps as xps
from Toa2Lai_PY import mcd19_prior_v006 as mcd19_prior_v006



def get_ozone_prior(image):
    image_date = image.date()
    tco3 = ee.ImageCollection('TOMS/MERGED').filterDate(ee.Date(image_date).advance(-3, 'day'), ee.Date(image_date).advance( 3, 'day')).mean()
    start_year = ee.Date(image_date).advance(-3, 'year').get('year')
    end_year   = ee.Date(image_date).advance( 3, 'year').get('year')
    month      = ee.Date(image_date).get('month')

    monthly_climatology = ee.ImageCollection('TOMS/MERGED').filter(ee.Filter.calendarRange(start_year, end_year, 'year')).filter(ee.Filter.calendarRange(month , month, 'month')).median().multiply(0.001)

    ozone = ee.Image(ee.Algorithms.If(tco3.bandNames(), tco3.select('ozone').multiply(0.001), monthly_climatology.select('ozone')))

    tco3 = ozone.unmask().where(ozone.mask().eq(0), monthly_climatology)

    ker = ee.Kernel.square(200, "pixels")
    tco3_mean = tco3.reduceNeighborhood(reducer=ee.Reducer.mean(), kernel=ker, inputWeight='mask', skipMasked=False, optimization='boxcar').rename('tco3')
    tco3 = tco3.unmask().where(tco3.mask().eq(0), tco3_mean)
    return tco3


def s2_cloud(image):

    image1 = image.multiply(0.0001).float()
    cloud_bands = ['B1', 'B2', 'B4', 'B5', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
    control = image1.select(cloud_bands)
    cloud_prob = Sen2Cloud.sen2cloud_prob(control)

    NDSI = image1.normalizedDifference(['B3', 'B11']).rename('NDSI')
    snow = NDSI.gt(0.45).And(image1.select('B8' ).gt(0.1500)).And(image1.select('B12').lt(0.12)).rename('snow')
    cloud = cloud_prob.select('cloud_prob').gt(0.6).And(snow.eq(0)).rename('cloud_mask')

    kernel1 = ee.Kernel.circle(radius=10, units='meters')
    kernel2 = ee.Kernel.circle(radius=10, units='meters')
    cloud = cloud.focal_min(kernel=kernel1, iterations=2).focal_max(kernel=kernel2, iterations=10)

    dataset = ee.Image('USGS/USGS/GMTED2010_FULL')
    elevation = dataset.select('be75').clip(image.geometry())

    thresh = ee.Image(-1.31071432e-05).multiply(elevation).add(elevation.pow(2).multiply(9.92142860e-09)).add(3.80892904e-03)
    thresh = thresh.add(0.01)
    b10 = image1.select('B10')
    cirrus = b10.gt(thresh).rename('cirrus')
    cloud = cloud.eq(1).Or(cirrus.eq(1))

    return image.updateMask(cloud.eq(0))

def s2_add_angs(image):
    saa = ee.Image(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE'))).rename('saa')
    sza = ee.Image(ee.Number(image.get('MEAN_SOLAR_ZENITH_ANGLE' ))).rename('sza')
    vza = ee.List([ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B1')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B2')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B3')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B4')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B5')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B6')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B7')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B8')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B8A')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B9')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B10')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B11')),
                        ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B12')),
                        ])
                    
    vaa = ee.List([ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B1')), 
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B2')),
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B3')), 
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B4')), 
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B5')), 
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B6')), 
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B7')), 
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8')), 
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8A')), 
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B9')), 
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B10')), 
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B11')), 
                        ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B12'))
                        ])
    vza = ee.Image.constant(vza.reduce(ee.Reducer.mean())).rename('vza')
    vaa = ee.Image.constant(vaa.reduce(ee.Reducer.mean())).rename('vaa')

    image = image.addBands([sza.float(), saa.float(), vza.float(), vaa.float()])
    return image



def do_ac(xap, xbp, xcp, image, band):
    band = ee.String(band)
    toa = image.select(band).divide(10000)
    y = toa.multiply(xap.select(band)).subtract(xbp.select(band))
    sur = y.divide(y.multiply(xcp.select(band)).add(1))
    return sur


def get_sur(image):

    image = s2_add_angs(image)

    geom  = image.geometry()
    image_date = image.date()
    projection = image.select('B2').projection()
    crs = projection.crs()

    vza = image.select('vza')
    vaa = image.select('vaa')
    sza = image.select('sza')
    saa = image.select('saa')
    raa  = ee.Image(vaa).subtract(saa)
    zero_image = image.select('B2').multiply(0)

    deg2rad = ee.Number(math.pi).divide(ee.Number(180.0))
    cos_sza = (sza.multiply(deg2rad)).cos().rename('cos_sza')
    cos_vza = (vza.multiply(deg2rad)).cos().rename('cos_vza')
    cos_raa = (raa.multiply(deg2rad)).cos().rename('cos_raa')


    tcwv = mcd19_prior.get_wv_prior(image).divide(8)
    aot  = mcd19_prior.get_aot_prior(image).divide(3)

    gc.collect()

    tcwv_v006 = mcd19_prior_v006.get_wv_prior_v006(image).divide(8)
    aot_v006  = mcd19_prior_v006.get_aot_prior_v006(image).divide(3)

    gc.collect()

    aot = ee.Image([aot, aot_v006]).reduce(ee.Reducer.firstNonNull()).rename('aot')
    tcwv = ee.Image([tcwv, tcwv_v006]).reduce(ee.Reducer.firstNonNull()).rename('tcwv')

    gc.collect()

    tco3 = get_ozone_prior(image).rename('tco3')
    ele = ee.Image('USGS/USGS/GMTED2010_FULL').select('be75').multiply(0.00001).rename('ele')
    inp = ee.Image.cat([cos_sza, cos_vza, cos_raa, aot, tcwv, tco3, ele]).reproject(crs, None, 500)

    gc.collect()
    xap = xps.get_xap(inp)
    gc.collect()
    xbp = xps.get_xbp(inp)
    gc.collect()
    xcp = xps.get_xcp(inp)
    gc.collect()

    xap = xap.multiply(ee.Image([21.432022999999997, 19.288629999999998, 23.570948, 14.282488, 12.966325000000001, 11.919912, 9.017822, 8.096992, 10.290405, 5.0013630000000004, 6.812678])) \
                    .add(ee.Image([0.437156, 0.446027, 0.544881, 0.529572, 0.510017, 0.513989, 0.519447, 0.554209, 0.540255, 0.780309, 0.835003]))
    xbp = xbp.multiply(ee.Image([8.068029, 6.088372, 4.780875, 3.54993, 3.238609, 3.042271, 2.80878, 2.4387809999999996, 2.599056, 0.895949, 0.5270859999999999])) \
                    .add(ee.Image([0.023714, 0.015876, 0.00949, 0.004827, 0.00382, 0.003152, 0.002511, 0.001704, 0.002037, 0.000146, 4.3e-05]))
    xcp = xcp.multiply(ee.Image([0.25760399999999994, 0.22974, 0.230782, 0.222948, 0.218427, 0.211096, 0.201904, 0.185895, 0.19265, 0.09445300000000001, 0.050554] )) \
                    .add(ee.Image([0.073024, 0.050022, 0.030634, 0.015767, 0.012464, 0.010273, 0.008168, 0.005514, 0.006588, 0.000191, 4.2e-05]))

    B1 = do_ac(xap, xbp, xcp, image, 'B1')
    B2 = do_ac(xap, xbp, xcp, image, 'B2')
    B3 = do_ac(xap, xbp, xcp, image, 'B3')
    B4 = do_ac(xap, xbp, xcp, image, 'B4')
    B5 = do_ac(xap, xbp, xcp, image, 'B5')
    B6 = do_ac(xap, xbp, xcp, image, 'B6')
    B7 = do_ac(xap, xbp, xcp, image, 'B7')
    B8 = do_ac(xap, xbp, xcp, image, 'B8')
    B8A = do_ac(xap, xbp, xcp, image, 'B8A')
    B11 = do_ac(xap, xbp, xcp, image, 'B11')
    B12 = do_ac(xap, xbp, xcp, image, 'B12')
    sur = ee.Image([B1, B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12]).multiply(10000).int()

    inps = ee.Image.cat([B2, B3, B4, B5, B6, B7, B8A, B11, B12, cos_sza, cos_vza, cos_raa])

    lai = NN_prosail.nn_lai(inps)

    gc.collect()

    lai = (lai.pow(0.01).multiply(100).subtract(100)).multiply(-2).multiply(1000).int()

    sur = sur.addBands(lai).max(ee.Image(0))

    return ee.Image(sur.copyProperties(image).set('system:footprint',  image.get('system:footprint'))
                                            .set('system:index',      image.get('system:index'))
                                            .set('system:time_end',   image.get('system:time_end'))
                                            .set('system:time_start', image.get('system:time_start')))

# exports.inv_prosail = get_sur
# exports.s2_cloud = s2_cloud
