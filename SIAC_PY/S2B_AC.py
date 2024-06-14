import ee
from SIAC_PY import TWO_NN as tnn
from SIAC_PY import S2_RAD_Offsets as S2_offsets
from SIAC_PY import S2_RAD_Scales as S2_scales

# var tnn = require('users/marcyinfeng/SIAC/:/Two_NN')
# S2_offsets = require('users/marcyinfeng/SIAC/:/S2_RAD_Offsets')
# S2_scales  = require('users/marcyinfeng/SIAC/:/S2_RAD_Scales')



bands = ee.List(['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A',  'B11', 'B12'])

def predict(image, H1_scale, H1_offset, H2_scale, H2_offset, Out_scales, Out_offsets):
    arrayImage1D = image.toArray()
    arrayImage2D = arrayImage1D.toArray(1)
    imageAxis = 0
    bandAxis  = 1
    arrayLength = arrayImage2D.arrayLength(imageAxis)
    arrayImage2D    = arrayImage2D
    h1  = arrayImage2D.arrayTranspose().matrixMultiply(ee.Image(ee.Array(H1_scale))).add(ee.Image(ee.Array(H1_offset)).toArray(1).arrayTranspose())
    in1 = h1.max(0)
    h2  = in1.matrixMultiply(ee.Image(ee.Array(H2_scale))).add(ee.Image(ee.Array(H2_offset)).toArray(1).arrayTranspose())
    in2 = h2.max(0)
    out_scale_offsets = Out_scales.zip(Out_offsets)
    def get_out(out_scale_offset):
        Out_scale  = ee.List(out_scale_offset).get(0)
        Out_offset = ee.List(out_scale_offset).get(1)
        oup = in2.matrixMultiply(ee.Image(ee.Array(Out_scale)).toArray(1))
        out = oup.arrayProject([0]).arrayFlatten([['xx']]).add(ee.Image(ee.Number(ee.List(Out_offset).get(0))))
        return out
    ret = out_scale_offsets.map(get_out)
    return ret

def S2B_B01_AC(control, image):
    i=0
    xap = predict(control, S2_scales.S2B_xap_H1_scale,  S2_offsets.S2B_xap_H1_offset, 
                                S2_scales.S2B_xap_H2_scale,  S2_offsets.S2B_xap_H2_offset, 
                                S2_scales.S2B_xap_Out_scale, S2_offsets.S2B_xap_Out_offset)
                                
    xbp = predict(control, S2_scales.S2B_xbp_H1_scale,  S2_offsets.S2B_xbp_H1_offset, 
                                S2_scales.S2B_xbp_H2_scale,  S2_offsets.S2B_xbp_H2_offset, 
                                S2_scales.S2B_xbp_Out_scale, S2_offsets.S2B_xbp_Out_offset)
    
    xcp = predict(control, S2_scales.S2B_xcp_H1_scale,  S2_offsets.S2B_xcp_H1_offset, 
                                S2_scales.S2B_xcp_H2_scale,  S2_offsets.S2B_xcp_H2_offset, 
                                S2_scales.S2B_xcp_Out_scale, S2_offsets.S2B_xcp_Out_offset)
    
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa


def S2B_B02_AC(control, image):
    i=1
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                S2_offsets.S2B_xap_H1_offset.get(i),
                                S2_scales.S2B_xap_H2_scale.get(i),
                                S2_offsets.S2B_xap_H2_offset.get(i),
                                S2_scales.S2B_xap_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                S2_offsets.S2B_xbp_H1_offset.get(i),
                                S2_scales.S2B_xbp_H2_scale.get(i),
                                S2_offsets.S2B_xbp_H2_offset.get(i),
                                S2_scales.S2B_xbp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                S2_offsets.S2B_xcp_H1_offset.get(i),
                                S2_scales.S2B_xcp_H2_scale.get(i),
                                S2_offsets.S2B_xcp_H2_offset.get(i),
                                S2_scales.S2B_xcp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa


def S2B_B03_AC(control, image):
    i=2
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                    S2_offsets.S2B_xap_H1_offset.get(i),
                                    S2_scales.S2B_xap_H2_scale.get(i),
                                    S2_offsets.S2B_xap_H2_offset.get(i),
                                    S2_scales.S2B_xap_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                    S2_offsets.S2B_xbp_H1_offset.get(i),
                                    S2_scales.S2B_xbp_H2_scale.get(i),
                                    S2_offsets.S2B_xbp_H2_offset.get(i),
                                    S2_scales.S2B_xbp_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                    S2_offsets.S2B_xcp_H1_offset.get(i),
                                    S2_scales.S2B_xcp_H2_scale.get(i),
                                    S2_offsets.S2B_xcp_H2_offset.get(i),
                                    S2_scales.S2B_xcp_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                    
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa


def S2B_B04_AC(control, image):
    i=3
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                S2_offsets.S2B_xap_H1_offset.get(i),
                                S2_scales.S2B_xap_H2_scale.get(i),
                                S2_offsets.S2B_xap_H2_offset.get(i),
                                S2_scales.S2B_xap_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                S2_offsets.S2B_xbp_H1_offset.get(i),
                                S2_scales.S2B_xbp_H2_scale.get(i),
                                S2_offsets.S2B_xbp_H2_offset.get(i),
                                S2_scales.S2B_xbp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                S2_offsets.S2B_xcp_H1_offset.get(i),
                                S2_scales.S2B_xcp_H2_scale.get(i),
                                S2_offsets.S2B_xcp_H2_offset.get(i),
                                S2_scales.S2B_xcp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa



def S2B_B05_AC(control, image):
    i=4
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                S2_offsets.S2B_xap_H1_offset.get(i),
                                S2_scales.S2B_xap_H2_scale.get(i),
                                S2_offsets.S2B_xap_H2_offset.get(i),
                                S2_scales.S2B_xap_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                S2_offsets.S2B_xbp_H1_offset.get(i),
                                S2_scales.S2B_xbp_H2_scale.get(i),
                                S2_offsets.S2B_xbp_H2_offset.get(i),
                                S2_scales.S2B_xbp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                S2_offsets.S2B_xcp_H1_offset.get(i),
                                S2_scales.S2B_xcp_H2_scale.get(i),
                                S2_offsets.S2B_xcp_H2_offset.get(i),
                                S2_scales.S2B_xcp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa



def S2B_B06_AC(control, image):
    i=5
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                S2_offsets.S2B_xap_H1_offset.get(i),
                                S2_scales.S2B_xap_H2_scale.get(i),
                                S2_offsets.S2B_xap_H2_offset.get(i),
                                S2_scales.S2B_xap_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                S2_offsets.S2B_xbp_H1_offset.get(i),
                                S2_scales.S2B_xbp_H2_scale.get(i),
                                S2_offsets.S2B_xbp_H2_offset.get(i),
                                S2_scales.S2B_xbp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                S2_offsets.S2B_xcp_H1_offset.get(i),
                                S2_scales.S2B_xcp_H2_scale.get(i),
                                S2_offsets.S2B_xcp_H2_offset.get(i),
                                S2_scales.S2B_xcp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa



def S2B_B07_AC(control, image):
    i=6
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                S2_offsets.S2B_xap_H1_offset.get(i),
                                S2_scales.S2B_xap_H2_scale.get(i),
                                S2_offsets.S2B_xap_H2_offset.get(i),
                                S2_scales.S2B_xap_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                S2_offsets.S2B_xbp_H1_offset.get(i),
                                S2_scales.S2B_xbp_H2_scale.get(i),
                                S2_offsets.S2B_xbp_H2_offset.get(i),
                                S2_scales.S2B_xbp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                S2_offsets.S2B_xcp_H1_offset.get(i),
                                S2_scales.S2B_xcp_H2_scale.get(i),
                                S2_offsets.S2B_xcp_H2_offset.get(i),
                                S2_scales.S2B_xcp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa


def S2B_B08_AC(control, image):
    i=7
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                    S2_offsets.S2B_xap_H1_offset.get(i),
                                    S2_scales.S2B_xap_H2_scale.get(i),
                                    S2_offsets.S2B_xap_H2_offset.get(i),
                                    S2_scales.S2B_xap_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                    S2_offsets.S2B_xbp_H1_offset.get(i),
                                    S2_scales.S2B_xbp_H2_scale.get(i),
                                    S2_offsets.S2B_xbp_H2_offset.get(i),
                                    S2_scales.S2B_xbp_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                    S2_offsets.S2B_xcp_H1_offset.get(i),
                                    S2_scales.S2B_xcp_H2_scale.get(i),
                                    S2_offsets.S2B_xcp_H2_offset.get(i),
                                    S2_scales.S2B_xcp_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                    
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa


def S2B_B8A_AC(control, image):
    i=8
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                    S2_offsets.S2B_xap_H1_offset.get(i),
                                    S2_scales.S2B_xap_H2_scale.get(i),
                                    S2_offsets.S2B_xap_H2_offset.get(i),
                                    S2_scales.S2B_xap_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                    S2_offsets.S2B_xbp_H1_offset.get(i),
                                    S2_scales.S2B_xbp_H2_scale.get(i),
                                    S2_offsets.S2B_xbp_H2_offset.get(i),
                                    S2_scales.S2B_xbp_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                    S2_offsets.S2B_xcp_H1_offset.get(i),
                                    S2_scales.S2B_xcp_H2_scale.get(i),
                                    S2_offsets.S2B_xcp_H2_offset.get(i),
                                    S2_scales.S2B_xcp_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                    
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa



def S2B_B09_AC(control, image):
    i=9
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                S2_offsets.S2B_xap_H1_offset.get(i),
                                S2_scales.S2B_xap_H2_scale.get(i),
                                S2_offsets.S2B_xap_H2_offset.get(i),
                                S2_scales.S2B_xap_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                S2_offsets.S2B_xbp_H1_offset.get(i),
                                S2_scales.S2B_xbp_H2_scale.get(i),
                                S2_offsets.S2B_xbp_H2_offset.get(i),
                                S2_scales.S2B_xbp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                S2_offsets.S2B_xcp_H1_offset.get(i),
                                S2_scales.S2B_xcp_H2_scale.get(i),
                                S2_offsets.S2B_xcp_H2_offset.get(i),
                                S2_scales.S2B_xcp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa


def S2B_B10_AC(control, image):
    i=10
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                S2_offsets.S2B_xap_H1_offset.get(i),
                                S2_scales.S2B_xap_H2_scale.get(i),
                                S2_offsets.S2B_xap_H2_offset.get(i),
                                S2_scales.S2B_xap_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                S2_offsets.S2B_xbp_H1_offset.get(i),
                                S2_scales.S2B_xbp_H2_scale.get(i),
                                S2_offsets.S2B_xbp_H2_offset.get(i),
                                S2_scales.S2B_xbp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                S2_offsets.S2B_xcp_H1_offset.get(i),
                                S2_scales.S2B_xcp_H2_scale.get(i),
                                S2_offsets.S2B_xcp_H2_offset.get(i),
                                S2_scales.S2B_xcp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa


def S2B_B11_AC(control, image):
    i=11
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                    S2_offsets.S2B_xap_H1_offset.get(i),
                                    S2_scales.S2B_xap_H2_scale.get(i),
                                    S2_offsets.S2B_xap_H2_offset.get(i),
                                    S2_scales.S2B_xap_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                    S2_offsets.S2B_xbp_H1_offset.get(i),
                                    S2_scales.S2B_xbp_H2_scale.get(i),
                                    S2_offsets.S2B_xbp_H2_offset.get(i),
                                    S2_scales.S2B_xbp_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                    S2_offsets.S2B_xcp_H1_offset.get(i),
                                    S2_scales.S2B_xcp_H2_scale.get(i),
                                    S2_offsets.S2B_xcp_H2_offset.get(i),
                                    S2_scales.S2B_xcp_Out_scale.get(i),
                                    ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                    
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa


def S2B_B12_AC(control, image):
    i=12
    xap = tnn.predict(control,S2_scales.S2B_xap_H1_scale.get(i),
                                S2_offsets.S2B_xap_H1_offset.get(i),
                                S2_scales.S2B_xap_H2_scale.get(i),
                                S2_offsets.S2B_xap_H2_offset.get(i),
                                S2_scales.S2B_xap_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xap_Out_offset.get(i)).get(0))

    xbp = tnn.predict(control,S2_scales.S2B_xbp_H1_scale.get(i),
                                S2_offsets.S2B_xbp_H1_offset.get(i),
                                S2_scales.S2B_xbp_H2_scale.get(i),
                                S2_offsets.S2B_xbp_H2_offset.get(i),
                                S2_scales.S2B_xbp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xbp_Out_offset.get(i)).get(0))
                                
    xcp = tnn.predict(control,S2_scales.S2B_xcp_H1_scale.get(i),
                                S2_offsets.S2B_xcp_H1_offset.get(i),
                                S2_scales.S2B_xcp_H2_scale.get(i),
                                S2_offsets.S2B_xcp_H2_offset.get(i),
                                S2_scales.S2B_xcp_Out_scale.get(i),
                                ee.List(S2_offsets.S2B_xcp_Out_offset.get(i)).get(0))
                                
    y = image.select([bands.get(i)]).multiply(xap).subtract(xbp); 
    boa = y.divide(xcp.multiply(y).add(ee.Image(1))).rename([bands.get(i)])
    return boa


def SIAC_S2B_10m(control_vars, image):
    B01 = S2B_B01_AC(control_vars, image)
    B02 = S2B_B02_AC(control_vars, image)
    B03 = S2B_B03_AC(control_vars, image)
    B04 = S2B_B04_AC(control_vars, image)
    B05 = S2B_B05_AC(control_vars, image)
    B06 = S2B_B06_AC(control_vars, image)
    B07 = S2B_B07_AC(control_vars, image)
    B08 = S2B_B08_AC(control_vars, image)
    B8A = S2B_B8A_AC(control_vars, image)
    B09 = S2B_B09_AC(control_vars, image)
    B10 = S2B_B10_AC(control_vars, image)
    B11 = S2B_B11_AC(control_vars, image)
    B12 = S2B_B12_AC(control_vars, image)
    return ee.Image.cat([B02, B03, B04, B08])


def SIAC_S2B_B05_B07_20m(control_vars, image):
    B01 = S2B_B01_AC(control_vars, image)
    B02 = S2B_B02_AC(control_vars, image)
    B03 = S2B_B03_AC(control_vars, image)
    B04 = S2B_B04_AC(control_vars, image)
    B05 = S2B_B05_AC(control_vars, image)
    B06 = S2B_B06_AC(control_vars, image)
    B07 = S2B_B07_AC(control_vars, image)
    B08 = S2B_B08_AC(control_vars, image)
    B8A = S2B_B8A_AC(control_vars, image)
    B09 = S2B_B09_AC(control_vars, image)
    B10 = S2B_B10_AC(control_vars, image)
    B11 = S2B_B11_AC(control_vars, image)
    B12 = S2B_B12_AC(control_vars, image)
    return ee.Image.cat([B05, B06, B07])

def SIAC_S2B_B8A_B12_20m(control_vars, image):
    B01 = S2B_B01_AC(control_vars, image)
    B02 = S2B_B02_AC(control_vars, image)
    B03 = S2B_B03_AC(control_vars, image)
    B04 = S2B_B04_AC(control_vars, image)
    B05 = S2B_B05_AC(control_vars, image)
    B06 = S2B_B06_AC(control_vars, image)
    B07 = S2B_B07_AC(control_vars, image)
    B08 = S2B_B08_AC(control_vars, image)
    B8A = S2B_B8A_AC(control_vars, image)
    B09 = S2B_B09_AC(control_vars, image)
    B10 = S2B_B10_AC(control_vars, image)
    B11 = S2B_B11_AC(control_vars, image)
    B12 = S2B_B12_AC(control_vars, image)
    return ee.Image.cat([B8A, B11, B12])

def SIAC_S2B_60m(control_vars, image):
    B01 = S2B_B01_AC(control_vars, image)
    B02 = S2B_B02_AC(control_vars, image)
    B03 = S2B_B03_AC(control_vars, image)
    B04 = S2B_B04_AC(control_vars, image)
    B05 = S2B_B05_AC(control_vars, image)
    B06 = S2B_B06_AC(control_vars, image)
    B07 = S2B_B07_AC(control_vars, image)
    B08 = S2B_B08_AC(control_vars, image)
    B8A = S2B_B8A_AC(control_vars, image)
    B09 = S2B_B09_AC(control_vars, image)
    B10 = S2B_B10_AC(control_vars, image)
    B11 = S2B_B11_AC(control_vars, image)
    B12 = S2B_B12_AC(control_vars, image)
    return ee.Image.cat([B01, B09, B10])

def S2B_AC(control, image):
    xap = predict(control, S2_scales.S2B_xap_H1_scale,  S2_offsets.S2B_xap_H1_offset, 
                                S2_scales.S2B_xap_H2_scale,  S2_offsets.S2B_xap_H2_offset, 
                                S2_scales.S2B_xap_Out_scale, S2_offsets.S2B_xap_Out_offset)
                                
    xbp = predict(control, S2_scales.S2B_xbp_H1_scale,  S2_offsets.S2B_xbp_H1_offset, 
                                S2_scales.S2B_xbp_H2_scale,  S2_offsets.S2B_xbp_H2_offset, 
                                S2_scales.S2B_xbp_Out_scale, S2_offsets.S2B_xbp_Out_offset)

    xcp = predict(control, S2_scales.S2B_xcp_H1_scale,  S2_offsets.S2B_xcp_H1_offset, 
                                S2_scales.S2B_xcp_H2_scale,  S2_offsets.S2B_xcp_H2_offset, 
                                S2_scales.S2B_xcp_Out_scale, S2_offsets.S2B_xcp_Out_offset)
    ret   = []
    for i in range(11):
        y = image.select([bands.get(i)]).multiply(ee.Image(xap.get(i))).subtract(ee.Image(xbp.get(i)))
        boa = y.divide(ee.Image(xcp.get(i)).multiply(y).add(ee.Image(1))).rename([bands.get(i)])
        ret.append(boa)
    return ee.Image.cat(ret)

