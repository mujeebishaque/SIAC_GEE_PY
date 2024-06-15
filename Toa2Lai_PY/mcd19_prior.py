import ee

def get_best_aot(image):
    qa = image.select('AOD_QA')
    qa_mask = (1 << 8) + (1 << 9) + (1 << 10) + (1 << 11)
    aod_qa = qa.bitwiseAnd(qa_mask).rightShift(8)
    best_aot = aod_qa.eq(0)
    return image.updateMask(best_aot)

def get_wv_prior(image):
    image_date = ee.Date(image.get('system:time_start'))
    mcd19_wv = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES').filterDate(image_date.advance(-60, 'minute'), image_date.advance(60, 'minute')).map(get_best_aot).select('Column_WV').mean()

    bad_image = ee.Image(0)
    bad_image = bad_image.updateMask(bad_image.neq(0)).float()
    mcd19_wv = ee.Algorithms.If(mcd19_wv.bandNames(), mcd19_wv, bad_image)
    mcd19_wv = ee.Image(mcd19_wv)

    mcd19_wv = mcd19_wv.updateMask(mcd19_wv.gt(0))

    daily_mcd19_wv = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES').filterDate(image_date, image_date.advance(1, 'day')).map(get_best_aot).select('Column_WV').mean()

    bad_image = ee.Image(0)
    bad_image = bad_image.updateMask(bad_image.neq(0)).float()
    daily_mcd19_wv = ee.Algorithms.If(daily_mcd19_wv.bandNames(), daily_mcd19_wv, bad_image)
    daily_mcd19_wv = ee.Image(daily_mcd19_wv)

    dailymcd19_wv = daily_mcd19_wv.updateMask(daily_mcd19_wv.gt(0))
    ker = ee.Kernel.square(10, "pixels")
    mcd19_wv_mean = mcd19_wv.reduceNeighborhood(reducer=ee.Reducer.mean(), kernel=ker, inputWeight='mask', skipMasked=False, optimization='boxcar').rename('wv')

    mcd19_wv = mcd19_wv.unmask().where(mcd19_wv.mask().eq(0), mcd19_wv_mean)
    mcd19_wv = mcd19_wv.where(mcd19_wv.lte(0), daily_mcd19_wv)
    mcd19_wv = mcd19_wv.updateMask(mcd19_wv.gt(0)).divide(1000)
    ker = ee.Kernel.square(200, "pixels")
    wv_prior_mean = mcd19_wv.reduceNeighborhood(reducer=ee.Reducer.mean(), kernel=ker, inputWeight='mask', skipMasked=False, optimization='boxcar').rename('wv')

    wv_prior = mcd19_wv.unmask().where(mcd19_wv.mask().eq(0), wv_prior_mean)

    ker = ee.Kernel.square(10, "pixels")
    wv_prior = wv_prior.reduceNeighborhood(reducer=ee.Reducer.mean(), kernel=ker, inputWeight='mask', skipMasked=False, optimization='boxcar').rename('wv')

    return wv_prior.updateMask(wv_prior.gt(0)).clip(image.geometry()).rename('wv_prior')


def get_aot_prior(image):
    image_date = ee.Date(image.get('system:time_start'))
    mcd19_aot = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES').filterDate(image_date.advance(-60, 'minute'), image_date.advance(60, 'minute')).map(get_best_aot).select('Optical_Depth_055').mean()

    bad_image = ee.Image(0)
    bad_image = bad_image.updateMask(bad_image.neq(0)).float()
    mcd19_aot = ee.Algorithms.If(mcd19_aot.bandNames(), mcd19_aot, bad_image)
    mcd19_aot = ee.Image(mcd19_aot)

    mcd19_aot = mcd19_aot.updateMask(mcd19_aot.gt(0))

    daily_mcd19_aot = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES').filterDate(image_date, image_date.advance(1, 'day')).map(get_best_aot).select('Optical_Depth_055').mean()
    bad_image = ee.Image(0)
    bad_image = bad_image.updateMask(bad_image.neq(0)).float()
    daily_mcd19_aot = ee.Algorithms.If(daily_mcd19_aot.bandNames(), daily_mcd19_aot, bad_image)
    daily_mcd19_aot = ee.Image(daily_mcd19_aot)

    dailymcd19_aot = daily_mcd19_aot.updateMask(daily_mcd19_aot.gt(0))
    ker = ee.Kernel.square(10, "pixels")
    mcd19_aot_mean = mcd19_aot.reduceNeighborhood(reducer= ee.Reducer.mean(), kernel= ker, inputWeight= 'mask', skipMasked=False, optimization='boxcar').rename('aot')

    mcd19_aot = mcd19_aot.unmask().where(mcd19_aot.mask().eq(0), mcd19_aot_mean)
    mcd19_aot = mcd19_aot.where(mcd19_aot.lte(0), daily_mcd19_aot)
    mcd19_aot = mcd19_aot.updateMask(mcd19_aot.gt(0)).divide(1000)

    ker = ee.Kernel.square(200, "pixels")
    aot_prior_mean = mcd19_aot.reduceNeighborhood(reducer= ee.Reducer.mean(), kernel= ker, inputWeight= 'mask', skipMasked=False, optimization='boxcar').rename('aot')

    aot_prior = mcd19_aot.unmask().where(mcd19_aot.mask().eq(0), aot_prior_mean)

    ker = ee.Kernel.square(10, "pixels")
    aot_prior = aot_prior.reduceNeighborhood(reducer=ee.Reducer.mean(), kernel= ker, inputWeight= 'mask', skipMasked=False, optimization='boxcar').rename('aot')
    return aot_prior.updateMask(aot_prior.gt(0)).clip(image.geometry()).rename('aot_prior')


# exports.get_aot_prior = get_aot_prior
# exports.get_wv_prior = get_wv_prior
