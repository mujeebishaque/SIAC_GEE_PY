from pathlib import Path
import os, ee


current_dir = Path(__file__).resolve().parent.parent

service_account = 'x@x.iam.gserviceaccount.com'
credentials = ee.ServiceAccountCredentials(service_account, str(current_dir) + os.sep + 'ee-x.json')
ee.Authenticate()
ee.Initialize(credentials)

# _image = ee.Image('users/marcyinfeng/Prosail_H1_offset')
# print(_image.getInfo())

def get_2D(image, nb, ny):
    array   = ee.List(image.reduceRegion(ee.Reducer.toList(), image.geometry()).get('b1'))
    sub_len = array.size().divide(nb)
    list_2d = []
    for j in  range(nb):    
        start = sub_len.multiply(j)
        end   = sub_len.multiply(j).add(ny)
        sub = array.slice(start, end)
        list_2d.append(sub)
    return ee.List(list_2d).reverse()


def get_3D(image, nb, nx, ny):
    array = ee.List(image.reduceRegion(ee.Reducer.toList(), image.geometry()).get('b1'))
    sub_len = array.size().divide(nb)
    ret = []
    for j in range(nb):
        list_2d = []
        for i in range(nx):
            start = sub_len.multiply(j).add(i*ny)
            end = sub_len.multiply(j).add((i+1)*ny)
            sub = array.slice(start, end)
            list_2d.append(sub)
        ret.append(list_2d)
    return ee.List(ret).reverse()

def parse_out(image):
    dict = image.reduceRegion(ee.Reducer.toList(), image.geometry())
    ret  = [dict.get('b1')]
    return ret

Prosail_H1_offset  = get_2D(ee.Image('users/marcyinfeng/Prosail_H1_offset'),  1, 256).get(0)
Prosail_H2_offset  = get_2D(ee.Image('users/marcyinfeng/Prosail_H2_offset'),  1, 256).get(0)
Prosail_Out_offset = parse_out(ee.Image('users/marcyinfeng/Prosail_Out_offset'))

Prosail_H1_scale  = get_3D(ee.Image('users/marcyinfeng/Prosail_H1_scale'),  1,  10, 256).get(0)
Prosail_H2_scale  = get_3D(ee.Image('users/marcyinfeng/Prosail_H2_scale'),  1, 256, 256).get(0)
Prosail_Out_scale = get_3D(ee.Image('users/marcyinfeng/Prosail_Out_scale'), 1, 256, 1).get(0)