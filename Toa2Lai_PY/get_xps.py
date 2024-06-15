import ee

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
    ws = nodes.multiply(nodes)
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

    def accum(item, list):
        start = ee.Number(item).multiply(nodes)
        end   = (ee.Number(item).add(1)).multiply(nodes)
        return ee.List(list).add(w3.slice(start, end))
    
    arr = ee.List([])
    w3 = ee.List.sequence(0, nodes.subtract(1), 1).iterate(accum, arr)

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
    h4  = in3.matrixMultiply(ee.Image(ee.Array(w3))).add(ee.Image(ee.Array(b3)).toArray(1).arrayTranspose())
    out = h4.arrayProject([1]).arrayFlatten([['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']])
    return out    

def get_xap(image):
    out = nn(image, ee.Image('users/marcyinfeng/6s_xap'), 11, 7)
    return out

def get_xbp(image):
    out = nn(image, ee.Image('users/marcyinfeng/6s_xbp'), 11, 7)
    return out

def get_xcp(image):
    out = nn(image, ee.Image('users/marcyinfeng/6s_xcp'), 11, 7)
    return out


# exports.get_xap = get_xap
# exports.get_xbp = get_xbp
# exports.get_xcp = get_xcp








