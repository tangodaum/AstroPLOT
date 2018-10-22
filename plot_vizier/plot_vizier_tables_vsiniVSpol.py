############################################
#######
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import numpy as np
import pyhdust.phc as phc
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import csv
import os
import re
from pyhdust import phc
from matplotlib import rc
from astropy import constants
rc('font', size=16)


# ==============================================================================
# Constants
pc = constants.pc.cgs  # Parsec
G = constants.G.value  # SI
M_sun = constants.M_sun.value  # SI
R_sun = constants.R_sun.value  # SI
L_sun = constants.L_sun.value  # SI
sigma_sb = constants.sigma_sb.value  # SI

read_again_from_simbad = False  # Do you whish to read everything again :-(?


# ==============================================================================
def vel_crit(r_pole, mass):

    vel_crit = 1e-3 * np.sqrt((2. * G * mass * M_sun) / (3. * r_pole * R_sun))
    # km/s
    return vel_crit

# ==============================================================================


Vizier.ROW_LIMIT = -1
cat = 'J/A+A/368/912/appen1'
# cat = 'J/A+A/368/912/appen2'

catalogs = Vizier.get_catalogs(cat)
# catalog1 = catalogs[0]
# catalog2 = catalogs[1]
catalog1 = catalogs[0]

# Operating with the data
data = catalog1.as_array()

# Print available data
print(data.dtype)

# Taking data
vsini = data['__vsini_'].data
evsini = data['e__vsini_'].data
pol_intr = data['Pintr'].data
sptyp = data['SpType'].data
hr = data['HR'].data

cmap = 'rainbow'
count = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
cor = phc.gradColor(count, cmapn=cmap, min=0, max=6)

# Plotting
# O, B, A, F, G, K, M
plt.clf()
for i in range(len(sptyp)):
    if 'O' in str(sptyp[i]):
        cor_temp = cor[0]
    elif 'B0' in str(sptyp[i]):
        cor_temp = cor[1]
    elif 'B1' in str(sptyp[i]):
        cor_temp = cor[2]
    elif 'B2' in str(sptyp[i]):
        cor_temp = cor[3]
    elif 'B3' in str(sptyp[i]):
        cor_temp = cor[4]
    elif 'B4' in str(sptyp[i]):
        cor_temp = cor[5]
    elif 'B5' in str(sptyp[i]):
        cor_temp = cor[6]
    elif 'B6' in str(sptyp[i]):
        cor_temp = cor[7]
    elif 'B7' in str(sptyp[i]):
        cor_temp = cor[8]
    elif 'B8' in str(sptyp[i]):
        cor_temp = cor[9]
    elif 'B9' in str(sptyp[i]):
        cor_temp = cor[10]

    if 'I' in str(sptyp[i]):
        marker = 's'
        label = 'I'
    if 'II' in str(sptyp[i]):
        marker = 'D'
        label = 'II'
    if 'III' in str(sptyp[i]):
        marker = 's'
        label = 'III'
    if 'V' in str(sptyp[i]):
        marker = 'o'
        label = 'V'
    if 'IV' in str(sptyp[i]):
        marker = '^'
        label = 'IV'
    # print(marker)

    # plt.scatter(x=vsini, y=oblat, marker=mark, markersize=5, color=cor_temp)
    plt.scatter(x=vsini[i], y=pol_intr[i], s=100, color=cor_temp, alpha=0.7,
                linewidths=1.5, edgecolors='black', marker=marker)

for i in range(len(sptyp)):
    if 6510 == hr[i]:
        plt.annotate(s=r'$\alpha$ Arae', xy=(vsini[i], pol_intr[i]),
                     xycoords='data', xytext=(400, 1.5),
                     arrowprops=dict(arrowstyle="<-",
                                     connectionstyle="angle, angleA=90," +
                                                     "angleB=0, rad=5"))

plt.xlabel(r'$v \sin \,i \, \rm [km/s]$', fontsize=18)
plt.ylabel(r'$p\,[\%]$', fontsize=18)

# plt.scatter(-1, -1, alpha=0.5, s=100, marker='s', label='I',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='D', label='II',
#             linewidths=1.5, edgecolors='black')
plt.scatter(-1, -1, alpha=0.5, s=100, marker='s', label='III',
            linewidths=1.5, edgecolors='black')
plt.scatter(-1, -1, alpha=0.5, s=100, marker='o', label='V',
            linewidths=1.5, edgecolors='black')
plt.scatter(-1, -1, alpha=0.5, s=100, marker='^', label='IV',
            linewidths=1.5, edgecolors='black')

plt.xlim(0, 500)
plt.ylim(0.1, 2.5)

norm = plt.Normalize(vmin=0, vmax=10)
s_m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
s_m.set_array([])
colourbar = plt.colorbar(s_m, ticks=list(count))
# colourbar.set_label('ST')
colourbar.set_ticklabels(['O', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                          'B7', 'B8', 'B9'])
plt.tight_layout()
plt.minorticks_on()
plt.legend(fontsize=14, loc=1, frameon=True, framealpha=0.5, scatterpoints=1)
plt.savefig('yudin_2001_fig.pdf')

plt.show()


# ------------------------------------------------------------------------------
count = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
cor = phc.gradColor(count, cmapn=cmap, min=0, max=9)

# Plotting
# O, B, A, F, G, K, M
plt.clf()

sptyp = sptyp.astype('str')
for i in range(len(sptyp)):
    print(sptyp[i])
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

        # obl = oblat[i] + 1
        vcrit = vel_crit(radius, mass)
        # print(vcrit)

        plt.scatter(x=vsini[i] / vcrit, y=pol_intr[i], s=100, color=cor_temp,
                    alpha=0.7, linewidths=1.5, edgecolors='black',
                    marker=marker)
        # print(vsini[i] / vcrit)
    # plt.scatter(x=vsini, y=oblat, marker=mark, markersize=5, color=cor_temp)
    # plt.scatter(x=vsini[i], y=pol_intr[i], s=100, color=cor_temp, alpha=0.7,
    #             linewidths=1.5, edgecolors='black', marker=marker)

for i in range(len(sptyp)):
    if 6510 == hr[i]:
        # print(i, sptyp[i])
        mass = 9.11  # msun
        radius = 5.33  # rsun
        vcrit = vel_crit(radius, mass)
        plt.annotate(s=r'$\alpha$ Arae', xy=(vsini[i] / vcrit, pol_intr[i]),
                     xycoords='data', xytext=(vsini[i] / vcrit + 0.1, 1.5),
                     arrowprops=dict(arrowstyle="<-",
                                     connectionstyle="angle, angleA=90," +
                                                     "angleB=0, rad=5"))

plt.xlabel(r'$v \sin \,i \, / v_{\rm crit}$', fontsize=18)
plt.ylabel(r'$p\,[\%]$', fontsize=18)

# plt.scatter(-1, -1, alpha=0.5, s=100, marker='s', label='I',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='D', label='II',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='s', label='III',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='o', label='V',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='^', label='IV',
#             linewidths=1.5, edgecolors='black')

plt.xlim(0, 1)
plt.ylim(0.1, 2.5)

norm = plt.Normalize(vmin=0, vmax=10)
s_m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
s_m.set_array([])
colourbar = plt.colorbar(s_m, ticks=list(count))
# colourbar.set_label('ST')
colourbar.set_ticklabels(['O', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                          'B7', 'B8', 'B9'])
plt.tight_layout()
plt.minorticks_on()
plt.legend(fontsize=14, loc=1, frameon=True, framealpha=0.5, scatterpoints=1)
plt.savefig('yudin_2001_fig2.pdf')

plt.show()


# ------------------------------------------------------------------------------
count = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
cor = phc.gradColor(count, cmapn=cmap, min=0, max=9)


# ==============================================================================
def dist_simbad(star_name):
    '''
    query sptype from SIMBAD
    '''

    customSimbad = Simbad()
    customSimbad.get_votable_fields()
    customSimbad.add_votable_fields('plx')

    # print(customSimbad.get_votable_fields())
    result_table = customSimbad.query_object(star_name)
    try:
        plx = result_table['PLX_VALUE'][0]
        distance_pc = (1000. / plx)  # pc
        print(distance_pc)
        # print(vsini)
        return distance_pc
    except:
        # distance_pc = -1000
        return np.nan


# ==============================================================================
# Taking distances

# Saving vsini
if read_again_from_simbad is True:
    dist = []

    for i in range(len(hr)):
        star = 'hr' + str(hr[i])
        dis = dist_simbad(star)
        dist.append(dis)
        print(dis, star)
else:
    dist = np.loadtxt('tables/dists_vsiniXpol.csv')


# Plotting
# O, B, A, F, G, K, M
plt.clf()

sptyp = sptyp.astype('str')
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

        # obl = oblat[i] + 1
        vcrit = vel_crit(radius, mass)
        # print(vcrit)

        plt.scatter(x=dist[i], y=pol_intr[i], s=100,
                    color=cor_temp, alpha=0.7, linewidths=1.5,
                    edgecolors='black', marker=marker)
        # print(vsini[i] / vcrit)
    # plt.scatter(x=vsini, y=oblat, marker=mark, markersize=5, color=cor_temp)
    # plt.scatter(x=vsini[i], y=pol_intr[i], s=100, color=cor_temp, alpha=0.7,
    #             linewidths=1.5, edgecolors='black', marker=marker)

for i in range(len(sptyp)):
    if 6510 == hr[i]:
        # print(i, sptyp[i])
        mass = 9.11  # msun
        radius = 5.33  # rsun
        vcrit = vel_crit(radius, mass)
        plt.annotate(s=r'$\alpha$ Arae', xy=(dist[i], pol_intr[i]),
                     xycoords='data', xytext=(dist[i] + 50, 1.5),
                     arrowprops=dict(arrowstyle="<-",
                                     connectionstyle="angle, angleA=90," +
                                                     "angleB=0, rad=5"))

# plt.xlabel(r'$v \sin \,i \, / v_{\rm crit}$', fontsize=18)
plt.xlabel(r'$distance [pc]$', fontsize=18)
plt.ylabel(r'$p\,[\%]$', fontsize=18)

# plt.scatter(-1, -1, alpha=0.5, s=100, marker='s', label='I',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='D', label='II',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='s', label='III',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='o', label='V',
#             linewidths=1.5, edgecolors='black')
# plt.scatter(-1, -1, alpha=0.5, s=100, marker='^', label='IV',
#             linewidths=1.5, edgecolors='black')

plt.xlim(0, 1600)
plt.ylim(-0.1, 2.)

norm = plt.Normalize(vmin=0, vmax=10)
s_m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
s_m.set_array([])
colourbar = plt.colorbar(s_m, ticks=list(count))
# colourbar.set_label('ST')
colourbar.set_ticklabels(['O', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                          'B7', 'B8', 'B9'])
plt.tight_layout()
plt.minorticks_on()
plt.legend(fontsize=14, loc=1, frameon=True, framealpha=0.5, scatterpoints=1)
plt.savefig('yudin_2001_fig3.pdf')

plt.show()