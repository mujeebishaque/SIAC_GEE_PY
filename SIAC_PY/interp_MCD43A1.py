import ee

mcd43_names = ['BRDF_Albedo_Parameters_Band1_iso', 
                'BRDF_Albedo_Parameters_Band2_iso', 
                'BRDF_Albedo_Parameters_Band3_iso', 
                'BRDF_Albedo_Parameters_Band4_iso', 
                'BRDF_Albedo_Parameters_Band6_iso', 
                'BRDF_Albedo_Parameters_Band7_iso',
                'BRDF_Albedo_Parameters_Band1_vol', 
                'BRDF_Albedo_Parameters_Band2_vol', 
                'BRDF_Albedo_Parameters_Band3_vol', 
                'BRDF_Albedo_Parameters_Band4_vol',
                'BRDF_Albedo_Parameters_Band6_vol', 
                'BRDF_Albedo_Parameters_Band7_vol', 
                'BRDF_Albedo_Parameters_Band1_geo', 
                'BRDF_Albedo_Parameters_Band2_geo', 
                'BRDF_Albedo_Parameters_Band3_geo', 
                'BRDF_Albedo_Parameters_Band4_geo', 
                'BRDF_Albedo_Parameters_Band6_geo', 
                'BRDF_Albedo_Parameters_Band7_geo']

def interp_mcd43a1(image):
    T_start = image.date().advance(-25, 'day')
    T_end   = image.date().advance( 25, 'day')
    geom    = image.geometry()
    date    = image.date()
    
    def update_mask(image):
        image = image.multiply(0.001).set({'system:time_start': image.get('system:time_start')}).copyProperties(image)
        return image
    
    mcd43a1 = ee.ImageCollection('MODIS/006/MCD43A1').filterDate(T_start, T_end).select(mcd43_names).map(update_mask)
    
    fitted_mcd43a1 = mcd43a1.median()
    return fitted_mcd43a1.select(mcd43_names[0:18])