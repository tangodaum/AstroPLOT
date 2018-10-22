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

Vizier.ROW_LIMIT = -1
# cat = ['J/MNRAS/361/1055/table1', 'J/A+A/393/887/SMC-Be']
cat = 'J/ApJ/722/605/table4'  # Huang 2010
catalogs = Vizier.get_catalogs(cat)
# catalog1 = catalogs[0]
# catalog2 = catalogs[1]
catalog1 = catalogs[0]

# Operating with the data
data = catalog1.as_array()

# Print available data
print(data.dtype)

# Taking data
vsini = data['vsini']
evsini = data['e_vsini']
teff = data['Teff']
eteff = data['e_Teff']
logg = data['log_g_']
elogg = data['e_log_g_']
hd = data['HD']

# Plotting
plt.errorbar(x=vsini, xerr=evsini, y=logg, yerr=elogg, linestyle='')
plt.show()


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
rc('font', size=16)

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
vsini = data['vsini']
oblat = data['Oblat']  # oblat = (b/a − 1)
theta = data['theta']  # stellar angular size
sptyp = data['SpType']
hd = data['HD']

cmap = 'rainbow'
count = np.array([0, 1, 2, 3, 4, 5, 6])
cor = phc.gradColor(count, cmapn=cmap, min=0, max=6)

# Plotting
# O, B, A, F, G, K, M
plt.clf()
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
        label = 'I'
    if 'II' in str(sptyp[i]):
        marker = 'D'
        label = 'Ii'
    if 'III' in str(sptyp[i]):
        marker = 'v'
        label = 'III'
    if 'V' in str(sptyp[i]):
        marker = 'o'
        label = 'V'
    if 'IV' in str(sptyp[i]):
        marker = '^'
        label = 'IV'
    print(marker)

    # plt.scatter(x=vsini, y=oblat, marker=mark, markersize=5, color=cor_temp)
    plt.scatter(x=vsini[i], y=oblat[i], s=100, color=cor_temp, alpha=0.7,
                linewidths=1.5, edgecolors='black', marker=marker)

for i in range(len(sptyp)):
    if 158427 == hd[i]:
        plt.annotate(s=r'$\alpha$ Arae', xy=(vsini[i], oblat[i]),
                     xycoords='data', xytext=(vsini[i], oblat[i] - 0.1),
                     arrowprops=dict(arrowstyle="<-",
                                     connectionstyle="angle, angleA=180," +
                                                     "angleB=270, rad=5"))

plt.xlabel(r'$v \sin \,i \, \rm [km/s]$', fontsize=18)
plt.ylabel(r'$(b/a) - 1$', fontsize=18)

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
plt.ylim(0, 0.5)

norm = plt.Normalize(vmin=0, vmax=6)
s_m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
s_m.set_array([])
colourbar = plt.colorbar(s_m, ticks=list(count))
# colourbar.set_label('ST')
colourbar.set_ticklabels(['O', 'B', 'A', 'F', 'G', 'K', 'M'])
plt.tight_layout()
plt.minorticks_on()
plt.legend(fontsize=14, loc=4, frameon=False, scatterpoints=1)
plt.savefig('yudin_2001_fig.pdf')

plt.show()





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
rc('font', size=16)

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