import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from astropy import units as u
import scipy.stats as stats
import scipy.optimize
from astropy import constants as c


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
peak_fits = pf.open('peak_RCrA_19jan.fit')
integ_data = integ_fits[0].data
peak_data = peak_fits[0].data
integ_data_flat = integ_data.reshape(81)
index_emission = np.argwhere((integ_data_flat-3*sigma2.reshape(81)) > 0)

W_co = integ_data_flat[index_emission[1:]]  # El 1er dato es una linea de base que no es plana
x_co = 2*10**20  # cm**-2(K*km/s)**-1 con 30% de incerteza
n_H2 = x_co*W_co*(u.cm)**-2
d = 130*u.pc
ang = 8.8*2*u.arcmin.to(u.rad)
area = (d**2*(ang)**2).to(u.cm**2)
N_tot = (np.sum(n_H2)*area)
m_H = 1.00794*u.u
m_H2 = 2.72*m_H  # correccion por abundancia de helio(paper Garcia)
M_tot = (m_H2*N_tot).to(u.M_sun)

"""Calculo mediante teorema del virial"""
datos_flat = datos.reshape(257, 81)
spec = datos_flat[49:113, index_emission[1:]]


def gaussian(x, amp, cen, wid):
    output = amp * np.exp(-(x-cen)**2 / wid)
    return output


def Modelo_gauss(x, a, b, A, mu, sigma):
    """
    entrega valores de la funcion del modelo gaussiano evaluada en punto
    o set de puntos
    """
    f = (a*x + b) - A * stats.norm(loc=mu, scale=sigma).pdf(x)
    return f


def ordenar(x_data, y_data, funcion, p_optimo):
    """
    Retorna arreglo ordenado de datos y de valores obtenidos con una funcion
    modelo
    """
    xmin = np.min(x_data)
    xmax = np.max(x_data)
    y_modelo_ordenado = np.sort(funcion(np.linspace(xmin, xmax, 200),
                                *p_optimo))
    y_data_ordenado = np.sort(y_data)
    return y_data_ordenado, y_modelo_ordenado


def prob_acumulada(y_data, y_modelo_ordenado):
    """
    Crea y entrega funcion de distribucion acumulada para el modelo
    """
    y_data_ordenado = np.sort(y_data)
    dist_acumulada_modelo = np.array([np.sum(y_modelo_ordenado <= yy) for yy in
                                     y_data_ordenado]) / len(y_modelo_ordenado)
    return dist_acumulada_modelo


desv_fit = np.zeros(len(spec[1, :]))
for i in range(0, len(spec[1, :])):
    x = vel[49:113]
    y = spec[:, i]

    a, b = np.polyfit(x, y, 1)
    A = np.max(y)
    mu = np.mean(x)  # Estoy dudoso de esta wea
    sigma = np.std(x)
    p0 = np.array([a[0], b[0], A, mu, sigma])
    p_optimo_gauss, a_covarianza_gauss = scipy.optimize.curve_fit(Modelo_gauss,
                                                                  x, y[:, 0], p0)
    desv_fit[i] = p_optimo_gauss[4]


def virial_mass(desv_fit, k, r):
    """desv_fit en km/s, r en pc, k es el coef de la distribucion de densidad
    rho(r) = r**-k"""
    mass = 3*(5 - 2*k)*r*desv_fit**2/((3-k)*c.G.to(u.pc*u.km**2/(u.M_sun*u.s**2)))
    output = np.sum(mass)
    return output


desv = desv_fit*u.km/u.s/np.sqrt(3)
r = np.sqrt(area/np.pi).to(u.pc)


masa_virial0 = virial_mass(desv, 0, r)
masa_virial1 = virial_mass(desv, 1, r)
masa_virial2 = virial_mass(desv, 2, r)

"""graficos de ejemplo"""


x = vel[49:113]
y = spec[:, 20]

a, b = np.polyfit(x, y, 1)
A = np.max(y)
mu = np.mean(x)  # Estoy dudoso de esta wea
sigma = np.std(x)
p0 = np.array([a[0], b[0], A, mu, sigma])
p_optimo_gauss, a_covarianza_gauss = scipy.optimize.curve_fit(Modelo_gauss,
                                                              x, y[:, 0], p0)
desv_fit[i] = p_optimo_gauss[4]
y_data_ordenado = np.sort(y)
CDF_modelo_gauss, y_modelo_ordenado_gauss = ordenar(x, y,
                                                    Modelo_gauss,
                                                    p_optimo_gauss)
Dn_gauss, prob_gauss = stats.kstest(y_data_ordenado, prob_acumulada,
                                    args=(y_modelo_ordenado_gauss,))
plt.clf()
plt.plot(x, Modelo_gauss(x, *p_optimo_gauss), color='r', label='fitting')
plt.plot(x, y, label='medicion')
plt.axvline(p_optimo_gauss[3], color='orange')
plt.axvline(p_optimo_gauss[3]+np.sqrt(2*np.log(2))*p_optimo_gauss[4], color='black')
plt.axvline(p_optimo_gauss[3]-np.sqrt(2*np.log(2))*p_optimo_gauss[4], color='black')
plt.show()
