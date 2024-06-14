import ee

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

def parse_out(image):
    dict = ee.Dictionary(image.reduceRegion(ee.Reducer.toList(), image.geometry()))
    ret  = ee.List([dict.get('b1'), dict.get('b2'), dict.get('b3'), dict.get('b4'),  dict.get('b5'),  dict.get('b6'), 
                        dict.get('b7'), dict.get('b8'), dict.get('b9'), dict.get('b10'), dict.get('b11')])
    return ret

S2A_xap_H1_offset  = get_2D(ee.Image('users/marcyinfeng/S2A_xap_H1_offset'),  1, 128).get(0)
S2A_xap_H2_offset  = get_2D(ee.Image('users/marcyinfeng/S2A_xap_H2_offset'),  1, 128).get(0)
S2A_xap_Out_offset = parse_out(ee.Image('users/marcyinfeng/S2A_shared_xap_Out_offset'))


S2A_xbp_H1_offset  = get_2D(ee.Image('users/marcyinfeng/S2A_xbp_H1_offset'),  1, 128).get(0)
S2A_xbp_H2_offset  = get_2D(ee.Image('users/marcyinfeng/S2A_xbp_H2_offset'),  1, 128).get(0)
S2A_xbp_Out_offset = parse_out(ee.Image('users/marcyinfeng/S2A_shared_xbp_Out_offset'))

S2A_xcp_H1_offset  = get_2D(ee.Image('users/marcyinfeng/S2A_xcp_H1_offset'),  1, 128).get(0)
S2A_xcp_H2_offset  = get_2D(ee.Image('users/marcyinfeng/S2A_xcp_H2_offset'),  1, 128).get(0)
S2A_xcp_Out_offset = parse_out(ee.Image('users/marcyinfeng/S2A_shared_xcp_Out_offset'))

S2B_xap_H1_offset  = get_2D(ee.Image('users/marcyinfeng/S2B_xap_H1_offset'),  1, 128).get(0)
S2B_xap_H2_offset  = get_2D(ee.Image('users/marcyinfeng/S2B_xap_H2_offset'),  1, 128).get(0)
S2B_xap_Out_offset = parse_out(ee.Image('users/marcyinfeng/S2B_shared_xap_Out_offset'))

S2B_xbp_H1_offset  = get_2D(ee.Image('users/marcyinfeng/S2B_xbp_H1_offset'),  1, 128).get(0)
S2B_xbp_H2_offset  = get_2D(ee.Image('users/marcyinfeng/S2B_xbp_H2_offset'),  1, 128).get(0)
S2B_xbp_Out_offset = parse_out(ee.Image('users/marcyinfeng/S2B_shared_xbp_Out_offset'))

S2B_xcp_H1_offset  = get_2D(ee.Image('users/marcyinfeng/S2B_xcp_H1_offset'),  1, 128).get(0)
S2B_xcp_H2_offset  = get_2D(ee.Image('users/marcyinfeng/S2B_xcp_H2_offset'),  1, 128).get(0)
S2B_xcp_Out_offset = parse_out(ee.Image('users/marcyinfeng/S2B_shared_xcp_Out_offset'))
