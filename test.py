"""Esta es una prueba, pa avanzar un poco..uso datos que se usaron en astro
experimental, tomo el punto 15, 150 porq es el que mas me tinca """

import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt

cubo = pf.open('Cubo_de_datos.fits')
datos = cubo[0].data
# Parece q el 1er eje es la latitud, la 2da es la longitud y la 3era es la vel
headers = cubo[0].header

#Para el fitting la que mas me tinco es la grafica con lat=15, lon=150
""" Otro que igual es bueno pa jugar es el 11,11,: del DHT03_RCrA_interp.fits"""
