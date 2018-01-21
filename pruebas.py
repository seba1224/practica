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
for i in range(0, 9):
    for j in range(0, 9):
        # Aca tengo la duda si se considera la linea de emision pa calcular el RMS
        temp_rms[i, j] = np.sqrt(1.0/len(vel)*np.sum(datos[:, i, j]**2))
        j += 1
    j = 0
    i += 1
sigma = np.sqrt(len(vel))*temp_rms*delta_vel
