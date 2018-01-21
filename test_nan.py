import numpy as np
import astropy.io.fits as pf

archivo_fit = pf.open('RCrA_19Jan.fit')
datos = archivo_fit[0].data
print(np.argwhere(np.isnan(datos)))
