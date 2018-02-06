import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt

cubo_thin = pf.open('integ_RCrA_19Jan.fit')
cubo_tololo = pf.open('integ_RCrA_tololo.fits')
cubo_3d_tololo = pf.open('DHT03_RCrA_raw.fits')
datos_integ_tololo = cubo_tololo[0].data
head_tololo = cubo_3d_tololo[0].header





def values(h, j):
    """Entrega los valores correctos de los ejes de un archivo fits"""
    N = h['NAXIS' + str(j)]
    val = np.zeros(N)
    for i in range(0, N):
        val[i] = (i + 1 - float(h['CRPIX' + str(j)])) * float(h['CDELT' + str(j)]) + float(h['CRVAL' + str(j)])
    return val


lat_tololo = values(head_tololo, 1)  # l
lon_tololo = values(head_tololo, 2)  # b
vel_tololo = values(head_tololo, 3)  # vel
delta_vel = np.abs(vel_tololo[1]-vel_tololo[2])


fig_1 = plt.figure()
ax1 = fig_1.add_subplot(111)
mapa_integ_thin = ax1.imshow(datos_integ_tololo, extent=[lat_tololo[0], lat_tololo[-1], lon_tololo[0], lon_tololo[-1]])
plt.show()

fig_2 = plt.figure()
ax2 = fig_2.add_subplot(111)
mapa_integ = ax2.imshow(datos_integ_tololo, origin='lower')
plt.show()
