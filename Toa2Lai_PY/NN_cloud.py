import ee
from Toa2Lai_PY import NN_prosail as NN_prosail
# NN_prosail  = require('users/marcyinfeng/Algorithms/:/NN_prosail')

def sen2cloud_prob(image):
    ret = ee.List(NN_prosail.get_nn(ee.Image('users/marcyinfeng/update_s2_cloud_new_new'), 10, 10))
    w0 = ee.List(ret.get(0))
    b0 = ee.List(ret.get(1))
    w1 = ee.List(ret.get(2))
    b1 = ee.List(ret.get(3))
    w2 = ee.List(ret.get(4))
    b2 = ee.List(ret.get(5))
    w3 = ee.List(ret.get(6))
    b3 = ee.List(ret.get(7))

    arrayImage1D = image.toArray()
    arrayImage2D = arrayImage1D.toArray(1);
    imageAxis = 0
    bandAxis  = 1
    arrayLength = arrayImage2D.arrayLength(imageAxis)
    h1  = arrayImage2D.arrayTranspose().matrixMultiply(ee.Image(ee.Array(w0))).add(ee.Image(ee.Array(b0)).toArray(1).arrayTranspose())
    in1 = h1.max(0)

    h2  = in1.matrixMultiply(ee.Image(ee.Array(w1))).add(ee.Image(ee.Array(b1)).toArray(1).arrayTranspose())
    in2 = h2.max(0)

    h3  = in2.matrixMultiply(ee.Image(ee.Array(w2))).add(ee.Image(ee.Array(b2)).toArray(1).arrayTranspose())
    in3 = h3.max(0)

    h4  = in3.matrixMultiply(ee.Image(ee.Array(w3)).arrayTranspose()).add(ee.Image(ee.Array(b3)).toArray(1).arrayTranspose())
    prox_exp = ee.Image(2).pow(h4.divide(ee.Image(-0.6931)))
    cloud_prob = ee.Image(1).divide(ee.Image(1).add(prox_exp))
    cloud_prob = cloud_prob.arrayProject([0]).arrayFlatten([['cloud_prob']])
    cloud_seed = cloud_prob.gt(0.44)
    potential_cloud = cloud_prob.gt(0.2)
    kernel = ee.Kernel.circle(radius=10, units='meters')
    cloud = cloud_seed.focal_min(kernel=kernel, iterations=2).focal_max(kernel=kernel, iterations=5).eq(1).And(potential_cloud.eq(1)).rename('cloud')
    return cloud_prob.addBands(cloud)

# exports.sen2cloud_prob  = sen2cloud_prob
