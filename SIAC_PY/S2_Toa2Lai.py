import ee, math
from pathlib import Path
from SIAC_PY import prosail as Prosail
from SIAC_PY import AC_processor as siac_processor
current_dir = Path(__file__).resolve().parent

# prosail = require('./Prosail')
# siac_processor = require('C/\Users\Dell\OneDrive\Desktop\PS_Work\wimex-im-cdsm-cesbio/SIAC///AC_processor')

def inv_prosail(image):
    image     = ee.Image(image)
    # geom  = image.geometry()
    # image_date = image.date()
    # projection = image.select('B2').projection()
    # crs = projection.crs()
    boa = siac_processor.get_boa(image)
    saa  = ee.Image.constant(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE'))).rename('saa')
    sza  = ee.Image.constant(ee.Number(image.get('MEAN_SOLAR_ZENITH_ANGLE' ))).rename('sza')
    vaa  = ee.Image.constant(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B2'))).rename('vaa')
    vza  = ee.Image.constant(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B2' ))).rename('vza')
    raa  = vaa.subtract(saa)
    deg2rad = ee.Number(math.pi).divide(ee.Number(180.0))
    cos_sza = (sza.multiply(deg2rad)).cos()
    cos_vza = (vza.multiply(deg2rad)).cos()
    cos_raa = (raa.multiply(deg2rad)).cos()
    name     = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A',  'B11', 'B12']
    new_name = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A',  'B11', 'B12']
    prosail_bands = ['B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B8A']
    prosail_inputs = ee.Image.cat([boa.select(name, new_name).select(prosail_bands), cos_sza, cos_vza, cos_raa])

    H1_scale = Prosail.Prosail_H1_scale
    H2_scale = Prosail.Prosail_H2_scale
    Out_scale = Prosail.Prosail_Out_scale

    H1_offset = Prosail.Prosail_H1_offset
    H2_offset = Prosail.Prosail_H2_offset
    Out_offset = Prosail.Prosail_Out_offset

    arrayImage1D = prosail_inputs.toArray()
    arrayImage2D = arrayImage1D.toArray(1)
    imageAxis = 0
    # bandAxis  = 1
    # arrayLength = arrayImage2D.arrayLength(imageAxis)
    arrayImage2D    = arrayImage2D
    h1  = arrayImage2D.arrayTranspose().matrixMultiply(ee.Image(ee.Array(H1_scale))).add(ee.Image(ee.Array(H1_offset)).toArray(1).arrayTranspose())
    in1 = h1.max(0)
    h2  = in1.matrixMultiply(ee.Image(ee.Array(H2_scale))).add(ee.Image(ee.Array(H2_offset)).toArray(1).arrayTranspose())
    in2 = h2.max(0)
    oup = in2.matrixMultiply(ee.Image(ee.Array(Out_scale)).toArray(1))
    out = oup.add(ee.Image(ee.Array(Out_offset[0])).toArray(1).arrayTranspose()).arrayProject([1]).arrayFlatten([['Lai']])
    Lai    = out.select('Lai').log().multiply(-2)
    # fname = ee.String(image.get('PRODUCT_ID'))
    # fname = image.get('PRODUCT_ID')
    return Lai.multiply(100).toInt()
    # ee.batch.Export.image.toDrive(collection=Lai.multiply(100).toInt(), description=str(fname), scale=10, folder='S2_SUR', fileFormat='GeoTIFF', region=geom, maxPixels=1e13)
    # ee.Export.image.toDrive({image: Lai.multiply(100).toInt(), description: fname + '_lai', scale: 10, fileFormat: 'GeoTIFF', folder:'S2_SUR', region: geom, maxPixels: 1e13})
