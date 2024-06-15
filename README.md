# [SIAC on GEE](https://code.earthengine.google.com/30176ec6495fd7d1d91f633ff34f1bbb)
### Feng Yin
### Department of Geography, UCL
### ucfafyi@ucl.ac.uk

[![DOI](https://zenodo.org/badge/185251518.svg)](https://zenodo.org/badge/latestdoi/185251518)
[![DOI](https://zenodo.org/badge/117815245.svg)](https://zenodo.org/badge/latestdoi/117815245)


Python implemetation of inv_prosail function. 

#### Use like this:
```
def atmoscorrect(s2_img):
    return S2_Toa2Lai.inv_prosail(s2_img).set('imageId', s2_img.id())
```

### Reference:

Yin, F., Lewis, P. E., Gomez-Dans, J., & Wu, Q. (2019, February 21). A sensor-invariant atmospheric correction method: application to Sentinel-2/MSI and Landsat 8/OLI. https://doi.org/10.31223/osf.io/ps957

An example UI:

![SIAC_GEE_UI](SIAC_GEE_UI.png)  

### LICENSE
GNU GENERAL PUBLIC LICENSE V3
