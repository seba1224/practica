import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt


"""Comparacion de los espectros integrados, ver q la integracion q hizo MacFITS
es solo la suma, falta multiplicar por delta v (parece q no..)"""
delta_vel = 0.31659999999999755
cubo_thin = pf.open('integ_RCrA_19Jan.fit')
cubo_wide = pf.open('integ_RCrA_wide_19Jan.fit')
integ_thin = cubo_thin[0].data
integ_wide = cubo_wide[0].data
head_integ_thin = cubo_thin[0].header
head_integ_wide = cubo_wide[0].header
fig_1 = plt.figure()
ax1 = fig_1.add_subplot(111)
mapa_integ_thin = ax1.imshow(integ_thin, extent=[-1, 1, -17, -19])
ax1.set_title('Espectro integrado espectrometro de alta resolucion')
cbar1 = fig_1.colorbar(mapa_integ_thin)

fig_2 = plt.figure()
ax2 = fig_2.add_subplot(111)
mapa_integ_wide = ax2.imshow(integ_wide, extent=[-1, 1, -17, -19])
ax2.set_title('Espectro integrado, espectrometro estandar')
cbar2 = fig_2.colorbar(mapa_integ_wide)


integ_thin_flat = integ_thin.reshape(81)
index_sort = np.argsort(integ_thin_flat)
integ_wide_flat = integ_wide.reshape(81)
fig_3 = plt.figure()
ax3 = fig_3.add_subplot(111)
ax3.plot(integ_wide_flat[index_sort], integ_thin_flat[index_sort], 'bo',
         label='datos medidos')
ax3.set_title('Integrated Temperature')
ax3.set_xlabel('Normal')
ax3.set_ylabel('Alta resolucion')

A = np.vstack([integ_wide_flat[index_sort], np.ones(len(integ_thin_flat))]).T
reg = np.linalg.lstsq(A, integ_thin_flat[index_sort])
ax3.plot(integ_wide_flat[index_sort], integ_wide_flat[index_sort]*reg[0][0] +
         reg[0][1], 'r', label='y=x$\\bullet$'+str(round(reg[0][0], 3))+'+' +
         str(round(reg[0][1], 3)))
ax3.legend()
corr_integ = np.corrcoef(integ_wide_flat[index_sort], integ_thin_flat[index_sort])

"""Se repite el procedimiento para los datos de peak temperature"""

fit_peak_thin = pf.open('peak_RCrA_19jan.fit')
fit_peak_wide = pf.open('peak_RCrA_wide_Jan19.fit')
peak_thin = fit_peak_thin[0].data
peak_wide = fit_peak_wide[0].data


fig_4 = plt.figure()
ax4 = fig_4.add_subplot(111)
mapa_peak_thin = ax4.imshow(peak_thin, extent=[-1, 1, -17, -19])
ax4.set_title('Peak de temperatura, espectrometro de alta resolucion')
cbar4 = fig_4.colorbar(mapa_peak_thin)
fig_5 = plt.figure()
ax5 = fig_5.add_subplot(111)
mapa_peak_wide = ax5.imshow(peak_wide, extent=[-1, 1, -17, -19])
ax5.set_title('Peak de temperatura, espectrometro estandar')
cbar5 = fig_5.colorbar(mapa_peak_wide)


peak_thin_flat = peak_thin.reshape(81)
index_sort2 = np.argsort(peak_thin_flat)
peak_wide_flat = peak_wide.reshape(81)
fig_6 = plt.figure()
ax6 = fig_6.add_subplot(111)
ax6.plot(peak_wide_flat[index_sort2], peak_thin_flat[index_sort2], 'bo',
         label='datos medidos')
ax6.set_title('Peak Temperature')
ax6.set_xlabel('Normal')
ax6.set_ylabel('Alta resolucion')

A2 = np.vstack([peak_wide_flat[index_sort2], np.ones(len(peak_thin_flat))]).T
reg2 = np.linalg.lstsq(A2, peak_thin_flat[index_sort2])
ax6.plot(peak_wide_flat[index_sort2], peak_wide_flat[index_sort2]*reg2[0][0] +
         reg2[0][1], 'r', label='y=x$\\bullet$'+str(round(reg2[0][0], 3))+'+' +
         str(round(reg2[0][1], 3)))
ax6.legend()
plt.show()
corr_peak = np.corrcoef(peak_thin_flat[index_sort2], peak_wide_flat[index_sort2])
