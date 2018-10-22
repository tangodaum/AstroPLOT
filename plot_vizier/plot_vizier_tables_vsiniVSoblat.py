
############################################
#######
import matplotlib.pyplot as plt
import numpy as np
import pyhdust.phc as phc
import pyhdust as phd
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import csv
import os
import re
from matplotlib import rc
from astropy import constants
from PyAstronomy import pyasl
rc('font', size=16)


# ==============================================================================
# Constants
pc = constants.pc.cgs  # Parsec
G = constants.G.value  # SI
M_sun = constants.M_sun.value  # SI
R_sun = constants.R_sun.value  # SI
L_sun = constants.L_sun.value  # SI
sigma_sb = constants.sigma_sb.value  # SI


# ==============================================================================
def vel_crit(r_pole, mass):

    vel_crit = 1e-3 * np.sqrt((2. * G * mass * M_sun) / (3. * r_pole * R_sun))
    # km/s
    return vel_crit


# ==============================================================================
def Upsilon(oblateness):
    '''
    Upsilon = veq / vcrit - Given the oblateness,
    it calculates the rotation rate (Upsilon)
    '''

    Upsilon = np.sqrt(3. - 3. / oblateness)
    return Upsilon
# ==============================================================================


Vizier.ROW_LIMIT = -1
cat = 'J/other/A+ARV/20.51/table4'

catalogs = Vizier.get_catalogs(cat)
# catalog1 = catalogs[0]
# catalog2 = catalogs[1]
catalog1 = catalogs[0]

# Operating with the data
data = catalog1.as_array()

# Print available data
print(data.dtype)

# Taking data
vsini = data['vsini'].data
oblat = data['Oblat'].data  # oblat = (b/a − 1)
theta = data['theta'].data  # stellar angular size
sptyp = data['SpType'].data
hd = data['HD'].data

sptyp = sptyp.astype('str')
hd = hd.astype('str')

cmap = 'rainbow'
count = np.array([0, 1, 2, 3, 4, 5, 6])
cor = phc.gradColor(count, cmapn=cmap, min=0, max=6)

# Plotting
# O, B, A, F, G, K, M
plt.clf()
upsilons = []
for i in range(len(sptyp)):
    if 'O' in str(sptyp[i]):
        cor_temp = cor[0]
    elif 'B' in str(sptyp[i]):
        cor_temp = cor[1]
    elif 'A' in str(sptyp[i]):
        cor_temp = cor[2]
    elif 'F' in str(sptyp[i]):
        cor_temp = cor[3]
    elif 'G' in str(sptyp[i]):
        cor_temp = cor[4]
    elif 'K' in str(sptyp[i]):
        cor_temp = cor[5]
    elif 'M' in str(sptyp[i]):
        cor_temp = cor[6]

    if 'I' in str(sptyp[i]):
        marker = 's'
        label = 'Iab'
    if 'II' in str(sptyp[i]):
        marker = 'D'
        label = 'II'
    if 'III' in str(sptyp[i]):
        marker = 'v'
        label = 'III'
    if 'V' in str(sptyp[i]):
        marker = 'o'
        label = 'V'
    if 'IV' in str(sptyp[i]):
        marker = '^'
        label = 'IV'
    # print(marker)
    st = sptyp[i][0:2]
    lum_class = np.copy(label)
    lum_class = lum_class.item()
    # print(sptyp[i])
    # print(st, lum_class)
    sdj = pyasl.SpecTypeDeJager()
    llum, lteff = sdj.lumAndTeff(st, lum_class)
    teff = 10**lteff
    obl = oblat[i] + 1
    ups = Upsilon(oblateness=obl)
    upsilons.append(ups)
    # print(st, lum_class, mass)
    # r_pole = ?
    # vcrit = vel_crit(r_pole, mass)
    # plt.scatter(x=vsini, y=oblat, marker=mark, markersize=5, color=cor_temp)
    plt.scatter(x=vsini[i], y=obl, s=100, color=cor_temp, alpha=0.7,
                linewidths=1.5, edgecolors='black', marker=marker)
    # plt.scatter(x=ups, y=obl, s=100, color=cor_temp, alpha=0.7,
    #             linewidths=1.5, edgecolors='black', marker=marker)


for i in range(len(sptyp)):
    # print(hd[i])
    if '158427' == hd[i]:
        # print(vsini[i], oblat[i] + 1)
        plt.annotate(s=r'$\alpha$ Arae', xy=(vsini[i], oblat[i] + 1),
                     xycoords='data', xytext=(vsini[i], oblat[i] + 1 - 0.1),
                     arrowprops=dict(arrowstyle="<-",
                                     connectionstyle="angle, angleA=180," +
                                                     "angleB=270, rad=5"))

plt.xlabel(r'$v \sin \,i \, \rm [km/s]$', fontsize=18)
plt.ylabel(r'$R_{\rm eq} / R_{\rm pole}$', fontsize=18)

plt.scatter(-1, -1, alpha=0.5, s=100, marker='s', label='I',
            linewidths=1.5, edgecolors='black')
plt.scatter(-1, -1, alpha=0.5, s=100, marker='D', label='II',
            linewidths=1.5, edgecolors='black')
plt.scatter(-1, -1, alpha=0.5, s=100, marker='v', label='III',
            linewidths=1.5, edgecolors='black')
plt.scatter(-1, -1, alpha=0.5, s=100, marker='o', label='V',
            linewidths=1.5, edgecolors='black')
plt.scatter(-1, -1, alpha=0.5, s=100, marker='^', label='IV',
            linewidths=1.5, edgecolors='black')

plt.xlim(0, 600)
plt.ylim(1, 1.5)

norm = plt.Normalize(vmin=0, vmax=6)
s_m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
s_m.set_array([])
colourbar = plt.colorbar(s_m, ticks=list(count))
# colourbar.set_label('ST')
colourbar.set_ticklabels(['O', 'B', 'A', 'F', 'G', 'K', 'M'])
plt.tight_layout()
plt.minorticks_on()
plt.legend(fontsize=14, loc=4, frameon=False, scatterpoints=1)
plt.savefig('van_Belle_oblat.pdf')

plt.show()


# ------------------------------------------------------------------------------
count = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
cor = phc.gradColor(count, cmapn=cmap, min=0, max=9)

# Plotting
# O, B, A, F, G, K, M
plt.clf()

for i in range(len(sptyp)):
    # print(sptyp[i])
    if ('B' == str(sptyp[i][0])) and \
       ('V' in str(sptyp[i])) and ('IV' not in str(sptyp[i])):
        if 'B0' in str(sptyp[i]):
            cor_temp = cor[0]
            mass = 17.5  # msun
            radius = 7.4  # rsun
        elif 'B1' in str(sptyp[i]):
            cor_temp = cor[1]
            mass = 13.21  # msun
            radius = 6.42  # rsun
        elif 'B2' in str(sptyp[i]):
            cor_temp = cor[2]
            mass = 9.11  # msun
            radius = 5.33  # rsun
        elif 'B3' in str(sptyp[i]):
            cor_temp = cor[3]
            mass = 7.6  # msun
            radius = 4.8  # rsun
        elif 'B4' in str(sptyp[i]):
            cor_temp = cor[4]
            mass = 6.62  # msun
            radius = 4.32  # rsun
        elif 'B5' in str(sptyp[i]):
            cor_temp = cor[5]
            mass = 5.8  # msun
            radius = 3.9  # rsun
        elif 'B6' in str(sptyp[i]):
            cor_temp = cor[6]
            mass = 5.17  # msun
            radius = 3.56  # rsun
        elif 'B7' in str(sptyp[i]):
            cor_temp = cor[7]
            mass = 4.45  # msun
            radius = 3.28  # rsun
        elif 'B8' in str(sptyp[i]):
            cor_temp = cor[8]
            mass = 3.8  # msun
            radius = 3.  # rsun
        elif 'B9' in str(sptyp[i]):
            cor_temp = cor[9]
            mass = 3.25  # msun
            radius = 2.7  # rsun

        marker = 'o'
        label = 'V'

        obl = oblat[i] + 1
        vcrit = vel_crit(radius, mass)
        # print(vcrit)

        plt.scatter(x=vsini[i] / vcrit, y=obl, s=100, color=cor_temp,
                    alpha=0.7, linewidths=1.5, edgecolors='black',
                    marker=marker)

for i in range(len(sptyp)):
    if '158427' == hd[i]:
        plt.annotate(s=r'$\alpha$ Arae', xy=(vsini[i], oblat[i] + 1),
                     xycoords='data', xytext=(vsini[i], oblat[i] + 1 - 0.1),
                     arrowprops=dict(arrowstyle="<-",
                                     connectionstyle="angle, angleA=180," +
                                                     "angleB=270, rad=5"))

plt.xlabel(r'$v \sin \,i \, / v_{\rm crit}$', fontsize=18)
plt.ylabel(r'$R_{\rm eq} / R_{\rm pole}$', fontsize=18)

# plt.scatter(-1, -1, alpha=0.5, s=100, marker='s', label='I',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='D', label='II',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='v', label='III',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='o', label='V',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='^', label='IV',
#             linewidths=1.5, edgecolors='black')

plt.vlines(x=1, ymin=1., ymax=1.5, linestyles='--', alpha=0.5)
plt.xlim(0.3, 1.2)
plt.ylim(1, 1.5)

norm = plt.Normalize(vmin=0, vmax=9)
s_m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
s_m.set_array([])
colourbar = plt.colorbar(s_m, ticks=list(count))
# colourbar.set_label('ST')
colourbar.set_ticklabels(['B0', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7',
                          'B8', 'B9'])
plt.tight_layout()
plt.minorticks_on()
plt.legend(fontsize=14, loc=4, frameon=False, scatterpoints=1)
plt.savefig('van_Belle_oblat2.pdf')

plt.show()


