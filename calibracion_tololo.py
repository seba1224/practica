import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt

cubo_thin = pf.open('integ_RCrA_19Jan.fit')
cubo_tololo = pf.open('integ_RCrA_tololo.fits')
cubo_3d_tololo = pf.open('DHT03_RCrA_raw.fits')
datos_integ_tololo = cubo_tololo[0].data
head_tololo = cubo_tololo[0].header
