import ee

def predict(image, H1_scale, H1_offset, H2_scale, H2_offset, Out_scale, Out_offset):
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
    oup = in2.matrixMultiply(ee.Image(ee.Array(Out_scale)))
    out = oup.arrayProject([0]).arrayFlatten([['xx']]).add(ee.Image(ee.Number(Out_offset)))

    return out