import ee, math
from SIAC_PY import TWO_NN as tnn
from SIAC_PY import Sen2Cloud as Sen2Cloud
# Sen2Cloud = require('users/marcyinfeng/SIAC/:/Sen2Cloud')


def get_cloud_prob(image):
    projection = image.select('B2').projection()
    crs = projection.crs()
    cloud_bands = ['B1', 'B2', 'B4', 'B5', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
    control = image.select(cloud_bands).multiply(0.0001).reproject(crs, None, 10)
    cloud = tnn.predict(control, Sen2Cloud.Sen2Cloud_H1_scale,
                            Sen2Cloud.Sen2Cloud_H1_offset,
                            Sen2Cloud.Sen2Cloud_H2_scale,
                            Sen2Cloud.Sen2Cloud_H2_offset,
                            Sen2Cloud.Sen2Cloud_Out_scale,
                            Sen2Cloud.Sen2Cloud_Out_offset).rename('cloud_prob')
    cloud_prob = ee.Image(1).divide(ee.Image(1).add((ee.Image(-1).multiply(cloud)).exp()))
    return image.addBands(cloud_prob.rename('cloud_prob'))

s2_bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'B10', 'B8A', 'B7']
l8_bands = ['B2', 'B3', 'B4', 'B5', 'B6',  'B7',  'B9',  'B10', 'B8']
l7_bands = ['B1', 'B2', 'B3', 'B4', 'B5',  'B7']
common_bands = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'CIRRUS', 'EB1', 'EB2'] 

def PCP_Sentinel2(image):
    NDSI = image.normalizedDifference(['GREEN', 'SWIR1']).rename('NDSI')
    NDVI = image.normalizedDifference(['NIR',     'RED']).rename('NDVI')
    basic_test = image.select('SWIR2').gt(0.03).And(NDSI.lt(0.8)).And(NDVI.lt(0.8))
    mean_vis   = (image.select('BLUE').add(image.select('GREEN')).add(image.select('RED'))).divide(3)
    whiteness = ((image.select('BLUE' ).subtract(mean_vis)).divide(mean_vis)).abs().add(((image.select('GREEN').subtract(mean_vis)).divide(mean_vis)).abs()).add(((image.select('RED').subtract(mean_vis)).divide(mean_vis)).abs()).rename('WHITENESS')
    whiteness_test =  whiteness.lt(0.7)
    HOT_test = (image.select('BLUE' ).subtract(image.select('RED').multiply(0.5)).subtract(0.08)).gt(0)
    NIR_SWIR1_Test =  (image.select('NIR').divide(image.select('SWIR1'))).gt(0.75)
    Water_Test = (NDVI.lt(0.01).And(image.select('NIR').lt(0.11))).Or(NDVI.lt( 0.1).And(image.select('NIR').lt(0.05))).rename('WATER_TEST')
    NIR = image.select('NIR')
    NIR = NIR.where(NIR.lt(0.12), 0.12)
    pcp = basic_test.And(whiteness_test).And(HOT_test).And(NIR_SWIR1_Test).rename('PCP')
    return image.addBands(pcp).addBands(Water_Test).addBands(whiteness).addBands(NDVI).addBands(NDSI)

def CDI(image):
    R_8A_8 = image.select('NIR').divide(image.select('EB1')).rename('R_8A_8')
    R_8A_7 = image.select('EB2').divide(image.select('EB1')).rename('R_8A_7')

    V_8A_8 = R_8A_8.reduceNeighborhood(reducer=ee.Reducer.stdDev(), kernel= ee.Kernel.circle(radius=140, units='meters'))
    V_8A_7 = R_8A_7.reduceNeighborhood(reducer= ee.Reducer.stdDev(), kernel= ee.Kernel.circle(radius=140, units='meters'))

    cdi = (V_8A_7.subtract(V_8A_8)).divide(V_8A_7.add(V_8A_8)).rename('CDI')
    return image.addBands(cdi)

def CDI_mask(image):
    image = CDI(image)
    CDI_Cloud = (image.select('PCP').And(image.select('CDI').lt(-0.25))).rename('CDI_Cloud')
    return image.addBands([CDI_Cloud]).updateMask(CDI_Cloud.eq(0)) 

def s2_scale(image):
    return image.divide(10000).cast({'BLUE'  : 'float'}) \
                            .cast({'GREEN' : 'float'}) \
                            .cast({'RED'   : 'float'}) \
                            .cast({'NIR'   : 'float'}) \
                            .cast({'SWIR1' : 'float'}) \
                            .cast({'SWIR2' : 'float'}) \
                            .cast({'CIRRUS': 'float'}) \
                            .cast({'EB1'   : 'float'}) \
                            .cast({'EB2'   : 'float'}) \
                            .clip(image.geometry()) \
                            .copyProperties(image) \
                            .set({'system:time_start': image.get('system:time_start')})


def L8_QA_clouds(image):
    qa = image.select('BQA')

    cloudBitMask  = 1 << 4
    cirrusBitMask = 1 << 6

    mask = (qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))).eq(1).rename('QA_cloud')
    return image.updateMask(mask)

def get_shadow(image):
    NIR  = image.select('NIR').multiply(10000).int()
    SWIR = image.select('SWIR1').multiply(10000).int()
    NIR_backg = NIR.reduceRegion(reducer= ee.Reducer.percentile([17.5]), geometry= image.geometry(), scale=30, maxPixels= 1e9).get('NIR')
                
    SWIR_backg = SWIR.reduceRegion(reducer= ee.Reducer.percentile([17.5]), geometry= image.geometry(), scale=30, maxPixels= 1e9).get('SWIR1')
    NIR_backg = ee.Number(ee.Algorithms.If(NIR_backg&1, NIR_backg, 0))

    SWIR_backg = ee.Number(ee.Algorithms.If(SWIR_backg&1, SWIR_backg, 0))
    NIR =  NIR.where( NIR.mask().eq(0),  NIR_backg)
    NIR = ee.Algorithms.FMask.fillMinima( NIR)
    NIR = NIR.subtract(image.select('NIR'))
    SWIR = SWIR.where(SWIR.mask().eq(0), SWIR_backg)
    SWIR = ee.Algorithms.FMask.fillMinima(SWIR)
    SWIR = SWIR.subtract(image.select('SWIR1'))
    shadow_prob = NIR.min(SWIR).rename('shadow_prob')
    return image.addBands(shadow_prob)

def s2_add_angs(image):
    saa = ee.Image(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE'))).rename('saa')
    sza = ee.Image(ee.Number(image.get('MEAN_SOLAR_ZENITH_ANGLE' ))).rename('sza')
    vza = ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B2')) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B3'))) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B4'))) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B8'))) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B8A'))) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B11'))) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B12'))) \
                    
    vaa = ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B2')) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B3'))) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B4'))) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8'))) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8A'))) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B11'))) \
        .add(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B12')))
    vza = ee.Image(vza.divide(7)).rename('vza')
    vaa = ee.Image(vaa.divide(7)).rename('vaa')
    image = image.addBands([sza.float(), saa.float(), vza.float(), vaa.float()])
    return image

def projectShadows(cloudMask, cloudHeights, sza, saa):
    nominalScale = cloudMask.projection().nominalScale()
    azR  = ee.Number(saa).add(180).multiply(math.pi).divide(180.0)
    zenR = ee.Number(sza).multiply(math.pi).divide(180.0)
    def mapCloudHeights(cloudHeight):
        cloudHeight = ee.Number(cloudHeight)
        shadowCastedDistance = zenR.tan().multiply(cloudHeight)
        x = azR.sin().multiply(shadowCastedDistance).divide(nominalScale)
        y = azR.cos().multiply(shadowCastedDistance).divide(nominalScale)
        return cloudMask.changeProj(cloudMask.projection(), cloudMask.projection().translate(x, y))
    
    shadows = list(map(mapCloudHeights, cloudHeights))
    shadowMask = ee.ImageCollection.fromImages(shadows).max()
    return shadowMask

def get_cloud_prob(image):

    cloud_bands = ['B1', 'B2', 'B4', 'B5', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
    control = image.select(cloud_bands).multiply(0.0001)
    cloud = tnn.predict(control, Sen2Cloud.Sen2Cloud_H1_scale,
                        Sen2Cloud.Sen2Cloud_H1_offset,
                        Sen2Cloud.Sen2Cloud_H2_scale,
                        Sen2Cloud.Sen2Cloud_H2_offset,
                        Sen2Cloud.Sen2Cloud_Out_scale,
                        Sen2Cloud.Sen2Cloud_Out_offset).rename('cloud_prob')
    cloud_prob = ee.Image(1).divide(ee.Image(1).add((ee.Image(-1).multiply(cloud)).exp()))
    return cloud_prob.rename('cloud_prob')

def get_cloud(image):
    projection = image.select('B2').projection()
    crs = projection.crs()
    cloud_prob = get_cloud_prob(image)
    T_start = image.date().advance(-20, 'day')
    T_end   = image.date().advance( 20, 'day')
    roi = image.geometry()
    l8_dataset = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').filterDate(T_start, T_end).filterBounds(roi).map(L8_QA_clouds).select(l8_bands, common_bands)
    
    s2_dataset = ee.ImageCollection('COPERNICUS/S2') \
                        .filterDate(T_start, T_end) \
                        .filterBounds(roi) \
                        .filter(ee.Filter.eq('MGRS_TILE', image.get('MGRS_TILE'))).select(s2_bands, common_bands).map(s2_scale).map(PCP_Sentinel2).map(CDI_mask)
    this_image = image.select(s2_bands, common_bands)
    this_image     = ee.Image(s2_scale(this_image))
    this_image     = PCP_Sentinel2(this_image)
    this_image     = CDI(this_image)
    cdi       = this_image.select('CDI').unmask()
    NDSI      = ee.Image(this_image.select('NDSI')).unmask()
    PCP       = ee.Image(this_image.select('PCP')).unmask()
    CDI_cloud = ee.Image(CDI_mask(this_image).select('CDI_Cloud')).unmask()
    snow       = NDSI.gt(0.15).And(image.select('B11').gt(0.1100)).And(image.select('B8A'  ).gt(0.1000)).rename('snow')
    sel_bands = common_bands.slice(0,6)
    merged_datasets = s2_dataset.select(sel_bands).merge(l8_dataset.select(sel_bands))
    image = s2_add_angs(image)
    fitted = merged_datasets.median()
    angs = image.select(['sza', 'saa', 'vza', 'vaa'])
    fitted   = fitted.select(common_bands.slice(0, 6)).clip(image.geometry())
    image = image.select(s2_bands.slice(0,6), common_bands.slice(0,6)).divide(10000)
    diff = image.select(common_bands.slice(0,6)).subtract(fitted).select(common_bands.slice(0,6), ['diff_BLUE', 'diff_GREEN', 'diff_RED', 'diff_NIR', 'diff_SWIR1', 'diff_SWIR2'])
                    
    image    = image.addBands(fitted.select(common_bands.slice(0,4), ['fitted_BLUE', 'fitted_GREEN', 'fitted_RED', 'fitted_NIR']))

    bad_pixel = fitted.lt(0.0001).reduce(ee.Reducer.allNonZero())
    diff_mean = diff.select('diff_RED').add(diff.select('diff_GREEN')).add(diff.select('diff_BLUE' )).divide(3)

    diff_mean2 = diff.select('diff_NIR').add(diff.select('diff_SWIR1')).add(diff.select('diff_SWIR2' )).divide(3)

    diff_mean3 = diff.select('diff_NIR') \
                            .add(diff.select('diff_SWIR1')) \
                            .add(diff.select('diff_SWIR2')) \
                            .add(diff.select('diff_RED')) \
                            .add(diff.select('diff_GREEN')) \
                            .add(diff.select('diff_BLUE' )).rename('all_bands').reproject(crs, None, 10)

    kernel1 = ee.Kernel.circle(radius= 2, units= 'pixels')
    kernel2 = ee.Kernel.circle(radius= 4, units= 'pixels')
    snow_thresh_1 = 0.05
    land_cloud = (diff_mean.gt(0.025)).And(NDSI.lte(snow_thresh_1)).focal_min(kernel= kernel1).focal_max(kernel= kernel2)
    bad_pixel_cloud = bad_pixel.eq(1).And(CDI_cloud)
    land_cloud = (diff_mean.gt(0.025)).And(snow.eq(0)).And(bad_pixel.eq(0))
    potencial_cloud1 = (diff_mean.gt(0.05)).And(NDSI.lt(0.5)).And(bad_pixel.eq(0))
    potencial_cloud2 = (diff_mean.gt(0.1)).And(NDSI.lt(0.5)).And(bad_pixel.eq(0))
    snow_cloud = (diff_mean.gt(0.025).Or(diff_mean.lt(-0.025))).And(snow.eq(1)).And(CDI_cloud)
                                    
    cloud      = land_cloud.Or(snow_cloud) \
                            .Or(potencial_cloud1) \
                            .Or(potencial_cloud2) \
                            .Or(bad_pixel_cloud).focal_min(kernel= kernel1).focal_max(kernel= kernel2).rename('cloud')
    image = get_shadow(image)
    
    saa = ee.Number(this_image.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
    sza = ee.Number(this_image.get('MEAN_SOLAR_ZENITH_ANGLE' ))
    cloudHeights = ee.List.sequence(200, 5000, 100);
    shadowMask = projectShadows(cloud, cloudHeights, sza, saa, 10)
    shadowMask = shadowMask.And(diff_mean.lt(-0.01)).And(image.select('shadow_prob').gt(200)).rename('shadow')
    
    shadow = diff_mean.lt(-0.005).And(image.select('shadow_prob').gt(200)).And(cloud.eq(0)).focal_min(kernel= kernel1).focal_max(kernel= kernel2).rename('shadow')
    return image.addBands([cloud, NDSI, PCP, CDI_cloud, shadow, cdi, diff])
