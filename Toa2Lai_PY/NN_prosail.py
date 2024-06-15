import ee

# loadNN = require('users/marcyinfeng/utils/:/load_NN')

def get_nn(nn_image, nodes, nxs):
    nn_image = ee.Image(nn_image)
    nodes = ee.Number(nodes)
    nxs = ee.Number(nxs)
    dns = ee.List(nn_image.reduceRegion(ee.Reducer.toList(), nn_image.geometry()).get('b1'))

    start = ee.Number(0)
    ws = nxs.multiply(nodes)
    w0 = dns.slice(0, ws)
    bs = nodes
    end = start.add(ws).add(bs)
    b0 = dns.slice(start.add(ws), end)

    start = end
    ws = nodes.multiply(nodes)
    w1 = dns.slice(start, start.add(ws))
    bs = nodes
    end = start.add(ws).add(bs)
    b1 = dns.slice(start.add(ws), end)

    start = end
    ws = nodes.multiply(nodes)
    w2 = dns.slice(start, start.add(ws))
    bs = nodes
    end = start.add(ws).add(bs)
    b2 = dns.slice(start.add(ws), end)

    start = end
    ws = nodes
    w3 = dns.slice(start, start.add(ws))
    bs = nodes
    end = start.add(ws).add(bs)
    b3 = dns.slice(start.add(ws), end)


    def accum(item, list):
        start = ee.Number(item).multiply(nodes)
        end   = (ee.Number(item).add(1)).multiply(nodes)
        return ee.List(list).add(w0.slice(start, end))

    arr = ee.List([])
    w0 = ee.List.sequence(0, nxs.subtract(1), 1).iterate(accum, arr)

    def accum(item, list):
        start = ee.Number(item).multiply(nodes)
        end   = (ee.Number(item).add(1)).multiply(nodes)
        return ee.List(list).add(w1.slice(start, end))

    arr = ee.List([])
    w1 = ee.List.sequence(0, nodes.subtract(1), 1).iterate(accum, arr)

    def accum(item, list):
        start = ee.Number(item).multiply(nodes)
        end   = (ee.Number(item).add(1)).multiply(nodes)
        return ee.List(list).add(w2.slice(start, end))

    arr = ee.List([])
    w2 = ee.List.sequence(0, nodes.subtract(1), 1).iterate(accum, arr)

    w3 = ee.List([w3])
    b3 = ee.List([b3.get(0)])

    return [w0, b0, w1, b1, w2, b2, w3, b3]

def nn(image, nn_image, nodes, nx):

    ret = ee.List(get_nn(nn_image, nodes, nx))
    w0 = ee.List(ret.get(0))
    b0 = ee.List(ret.get(1))
    w1 = ee.List(ret.get(2))
    b1 = ee.List(ret.get(3))
    w2 = ee.List(ret.get(4))
    b2 = ee.List(ret.get(5))
    w3 = ee.List(ret.get(6))
    b3 = ee.List(ret.get(7))

    arrayImage1D = image.toArray()
    arrayImage2D = arrayImage1D.toArray(1)
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
    out = h4.arrayProject([0]).arrayFlatten([['b1']])
    return out    

def nn_lai(image):
    out = nn(image, ee.Image('users/marcyinfeng/prosail_lai'), 15, 12).rename('lai')
    return out
def nn_lai_accurate(image):
    nnImage = ee.Image('users/marcyinfeng/NN_models/prosail_lai_accurate_2')
    out = loadNN.nn(image, nnImage, 12, 1, 4).rename('lai')
    return out

def nn_cab(image):
    nnImage = ee.Image('users/marcyinfeng/prosail_cab_new')
    image = ee.Image(image).select([0,1,2,3,4,5,6,9,10,11])
    out = loadNN.nn(image, nnImage, 10, 1, 4).rename('cab')

    return out

def nn_cw(image):
    out = nn(image, ee.Image('users/marcyinfeng/prosail_cw'), 15, 12).rename('cw')
    return out

def nn_cbrown(image):
    out = nn(image, ee.Image('users/marcyinfeng/prosail_cbrown'), 15, 12).rename('cbrown')
    return out

def nn_cbrown_accurate(image):
    nnImage = ee.Image('users/marcyinfeng/NN_models/nnCbrown_accurate')
    out = loadNN.nn(image, nnImage,  13, 1, 4).rename('cbrown')
    return out

# exports.get_nn    = get_nn
# exports.nn_lai    = nn_lai
# exports.nn_cab    = nn_cab
# exports.nn_cw     = nn_cw
# exports.nn_cbrown = nn_cbrown
# exports.nn_lai_accurate = nn_lai_accurate