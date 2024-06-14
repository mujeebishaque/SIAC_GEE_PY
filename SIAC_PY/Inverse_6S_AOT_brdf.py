import ee

def get_3D(image, nb, nx, ny):
    array   = ee.List(image.reduceRegion(ee.Reducer.toList(), image.geometry()).get('b1'))
    sub_len = array.size().divide(nb)
    ret = []
    for j in range(nb):
        list_2d = []
        for i in range(nx):
            start = sub_len.multiply(j).add(i*ny)
            end   = sub_len.multiply(j).add((i+1)*ny)
            sub = array.slice(start, end)
            list_2d.append(sub)
        ret.append(list_2d)
    return ee.List(ret).reverse()

def get_2D(image, nb, ny):
    array   = ee.List(image.reduceRegion(ee.Reducer.toList(), image.geometry()).get('b1'))
    sub_len = array.size().divide(nb)
    list_2d = []
    for j in range(nb):
        start = sub_len.multiply(j)
        end   = sub_len.multiply(j).add(ny)
        sub = array.slice(start, end)
        list_2d.append(sub)
    return ee.List(list_2d).reverse()



Inverse_6S_AOT_brdf_H1_scale  = get_3D(ee.Image('users/marcyinfeng/Inverse_6S_AOT_brdf_H1_scale'),  1, 13, 64).get(0)
Inverse_6S_AOT_brdf_H2_scale  = get_3D(ee.Image('users/marcyinfeng/Inverse_6S_AOT_brdf_H2_scale'),  1, 64, 64).get(0)
Inverse_6S_AOT_brdf_Out_scale = get_3D(ee.Image('users/marcyinfeng/Inverse_6S_AOT_brdf_Out_scale'), 1, 64, 1).get(0)

Inverse_6S_AOT_brdf_H1_offset  = get_2D(ee.Image('users/marcyinfeng/Inverse_6S_AOT_brdf_H1_offset'),  1, 64).get(0)
Inverse_6S_AOT_brdf_H2_offset  = get_2D(ee.Image('users/marcyinfeng/Inverse_6S_AOT_brdf_H2_offset'),  1, 64).get(0)
Inverse_6S_AOT_brdf_Out_offset = get_2D(ee.Image('users/marcyinfeng/Inverse_6S_AOT_brdf_Out_offset'), 1, 1).get(0)
