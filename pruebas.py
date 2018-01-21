import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt


def values(h, j):
    """Entrega los valores correctos de los ejes de un archivo fits"""
    N = h['NAXIS' + str(j)]
    val = np.zeros(N)
    for i in range(0, N):
        val[i] = (i + 1 - float(h['CRPIX' + str(j)])) * float(h['CDELT' + str(j)]) + float(h['CRVAL' + str(j)])
    return val


cubo = pf.open('RCrA_19Jan.fit')
datos = cubo[0].data  # vel, lat, lon ...cuidado que los indices de la lat estan invertidos
# osea el 0 es el -19.00
header = cubo[0].header

lon = values(header, 1)  # l
lat = values(header, 2)  # b
vel = values(header, 3)  # vel
delta_vel = np.abs(vel[1]-vel[2])

# El canal 122 esta malo
datos[122, :, :] = (datos[121, :, :]+datos[123, :, :])

temp_rms = np.zeros([9, 9])
temp_rms2 = np.zeros([9, 9])  # no considera la linea de emision para cal el rms
for i in range(0, 9):
    for j in range(0, 9):
        # Aca tengo la duda si se considera la linea de emision pa calcular el RMS
        temp_rms[i, j] = np.sqrt(1.0/len(vel)*np.sum(datos[:, i, j]**2))
        temp_rms2[i, j] = np.sqrt(1.0/len(vel)*(np.sum(datos[0:49, i, j]**2) +
                                  np.sum(datos[113:, i, j]**2)))
        j += 1
    j = 0
    i += 1
sigma = np.sqrt(len(vel))*temp_rms*delta_vel
sigma2 = np.sqrt(len(vel))*temp_rms2*delta_vel


integ_fits = pf.open('integ_RCrA_19Jan.fit')
peak_fits = pf.open('peak_RCrA_19Jan.fit')
integ_data = integ_fits[0].data*delta_vel
peak_data = peak_fits[0].data
index_emission = np.argwhere((integ_data-3*sigma2) > 0)
# index_emission = np.argwhere((integ_data-3*sigma) > 0)
# index_peak = np.argwhere((peak_data-3*temp_rms2) > 0)
# index_peak = np.argwhere((peak_data-3*temp_rms) > 0)
