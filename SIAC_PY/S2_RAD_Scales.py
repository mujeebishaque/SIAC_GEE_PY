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

def parse_out(image):
    dict = ee.Dictionary(image.reduceRegion(ee.Reducer.toList(), image.geometry()))
    ret  = ee.List([dict.get('b1'), dict.get('b2'), dict.get('b3'), dict.get('b4'),  dict.get('b5'),  dict.get('b6'), 
                        dict.get('b7'), dict.get('b8'), dict.get('b9'), dict.get('b10'), dict.get('b11')])
    return ret


S2A_xap_H1_scale  =  get_3D(ee.Image('users/marcyinfeng/S2A_xap_H1_scale'),  1,  7, 128).get(0)
S2A_xap_H2_scale  =  get_3D(ee.Image('users/marcyinfeng/S2A_xap_H2_scale'),  1, 128, 128).get(0)
S2A_xap_Out_scale = parse_out(ee.Image('users/marcyinfeng/S2A_shared_xap_Out_scale'))


S2A_xbp_H1_scale  =  get_3D(ee.Image('users/marcyinfeng/S2A_xbp_H1_scale'),  1,  7, 128).get(0)
S2A_xbp_H2_scale  =  get_3D(ee.Image('users/marcyinfeng/S2A_xbp_H2_scale'),  1, 128, 128).get(0)
S2A_xbp_Out_scale = parse_out(ee.Image('users/marcyinfeng/S2A_shared_xbp_Out_scale'))

S2A_xcp_H1_scale  =  get_3D(ee.Image('users/marcyinfeng/S2A_xcp_H1_scale'),  1,  7, 128).get(0)
S2A_xcp_H2_scale  =  get_3D(ee.Image('users/marcyinfeng/S2A_xcp_H2_scale'),  1, 128, 128).get(0)
S2A_xcp_Out_scale = parse_out(ee.Image('users/marcyinfeng/S2A_shared_xcp_Out_scale'))

S2B_xap_H1_scale  =  get_3D(ee.Image('users/marcyinfeng/S2B_xap_H1_scale'),  1,  7, 128).get(0)
S2B_xap_H2_scale  =  get_3D(ee.Image('users/marcyinfeng/S2B_xap_H2_scale'),  1, 128, 128).get(0)
S2B_xap_Out_scale = parse_out(ee.Image('users/marcyinfeng/S2B_shared_xap_Out_scale'))

S2B_xbp_H1_scale  =  get_3D(ee.Image('users/marcyinfeng/S2B_xbp_H1_scale'),  1,  7, 128).get(0)
S2B_xbp_H2_scale  =  get_3D(ee.Image('users/marcyinfeng/S2B_xbp_H2_scale'),  1, 128, 128).get(0)
S2B_xbp_Out_scale = parse_out(ee.Image('users/marcyinfeng/S2B_shared_xbp_Out_scale'))

S2B_xcp_H1_scale  =  get_3D(ee.Image('users/marcyinfeng/S2B_xcp_H1_scale'),  1,  7, 128).get(0)
S2B_xcp_H2_scale  =  get_3D(ee.Image('users/marcyinfeng/S2B_xcp_H2_scale'),  1, 128, 128).get(0)
S2B_xcp_Out_scale = parse_out(ee.Image('users/marcyinfeng/S2B_shared_xcp_Out_scale'))
