import matplotlib.pyplot as plt
import pyhdust as phd
import numpy as np
import pyhdust.phc as phc
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from matplotlib import rc
from PyAstronomy import pyasl
import pandas
from collections import Counter
import collections
from pyhdust import phc
from astropy import constants
rc('font', size=16)

read_again_from_simbad = False  # Do you whish to read everything again :-(?


# C =	constant RV
# SB1 =	spectroscopic binary, primary spectrum only detected
# SB2 =	spectroscopic binary, 2 spectra detected
# U =	unknown

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
def vsini_from_simbad(star_name):
    '''
    query sptype from SIMBAD
    '''

    customSimbad = Simbad()
    customSimbad.get_votable_fields()
    customSimbad.add_votable_fields('plx', 'plx_error', 'rot', 'sptype',
                                    'measurements')

    # print(customSimbad.get_votable_fields())
    result_table = customSimbad.query_object(star_name)
    try:
        vsini = result_table['ROT_Vsini'][0]
        # print(vsini)
        return vsini
    except:
        vsini = -1000
        return vsini

# ==============================================================================


Vizier.ROW_LIMIT = -1
# cat = ['J/MNRAS/361/1055/table1', 'J/A+A/393/887/SMC-Be']
cat = '2012MNRAS.424.1925C'
catalogs = Vizier.get_catalogs(cat)
# catalog1 = catalogs[0]
# catalog2 = catalogs[1]
catalog1 = catalogs[0]

# Operating with the data
data = catalog1.as_array()

# Print available data
print(data.dtype)

# Taking data
sptyp = data['SpType'].data
name = data['Name'].data
simbad_name = data['SimbadName'].data
binarity = data['Mult'].data

name = name.astype('str')
simbad_name = simbad_name.astype('str')
sptyp = sptyp.astype('str')
binarity = binarity.astype('str')

# Saving vsini
if read_again_from_simbad is True:
    vsini = []
    for i in range(len(simbad_name)):
        # name[i] = name[i].replace(' ', '')
        # if i > 517:
        print(simbad_name[i], '{} %'.format((i / len(simbad_name) * 100)))
        vsini_tmp = vsini_from_simbad(simbad_name[i])
        vsini.append(vsini_tmp)
    np.savetxt('tables/table_vsini_binarity.csv', vsini)
else:
    vsini = np.loadtxt('tables/table_vsini_binarity.csv')

# Plotting
cmap = 'rainbow'
count = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

cor = phc.gradColor(count, cmapn=cmap, min=0, max=9)

# Plotting
# O, B, A, F, G, K, M
plt.clf()

ylimt = []
condition = True
suffix = 'SB1'  # 'PM', 'MS', 'GS', 'SG', 'SD', 'WD', 'CP'
suffix2 = 'SB2'
for i in range(len(sptyp)):
    if ('B' == str(sptyp[i][0])) and \
       ('V' in str(sptyp[i])) and ('IV' not in str(sptyp[i])):
        st = sptyp[i][0:2]
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

        biny = str(binarity[i])
        if 'SB2' in biny:
            # print('yes')
            marker = 's'
            label = 'SB2'
        elif 'SB1' in biny:
            marker = 'o'
            label = 'SB1'
        elif 'C' in biny:
            marker = 'v'
            label = 'C'
        elif 'U' in biny:
            marker = '^'
            label = 'U'

        lc = str(sptyp[i])
        if 'I' in lc:
            lum_class = 'I'
        if 'II' in lc:
            lum_class = 'II'
        if 'III' in lc:
            lum_class = 'III'
        if 'V' in lc:
            lum_class = 'V'
        if 'IV' in lc:
            lum_class = 'IV'

        marker = 'o'
        label = 'V'

        # obl = oblat[i] + 1
        vcrit = vel_crit(radius, mass)
        try:
            sdj = pyasl.SpecTypeDeJager()
            # print(st, lum_class)
            llum, lteff = sdj.lumAndTeff(st, lum_class)
            teff = 10**lteff

            if (suffix in biny or suffix2 in biny) and (condition is True):
                print('aqui')
                plt.scatter(x=vsini[i], y=teff / 1e4, s=100, color=cor_temp,
                            alpha=0.5, linewidths=1.5, edgecolors='black',
                            marker=marker)
            if condition is False:
                plt.scatter(x=vsini[i], y=teff / 1e4, s=100, color=cor_temp,
                            alpha=0.5, linewidths=1.5, edgecolors='black',
                            marker=marker)
        except:
            pass

plt.ylabel(r'$10^4 \times T_{\rm eff} \,  [K]$', fontsize=18)
plt.xlabel(r'$v \sin i \, \rm [km/s] $', fontsize=18)

# plt.xlim(0, np.nanmax(vsini))
# plt.yscale('log')
# plt.xscale('log')
plt.xlim(-40, 400)
plt.ylim(0.8, 5)

norm = plt.Normalize(vmin=0, vmax=10)
s_m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
s_m.set_array([])
colourbar = plt.colorbar(s_m, ticks=list(count))
# colourbar.set_label('ST')
colourbar.set_ticklabels(['O', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                          'B7', 'B8', 'B9'])

plt.tight_layout()
plt.minorticks_on()
plt.legend(fontsize=14, loc=4, frameon=False, scatterpoints=1)
plt.savefig('vsini_VS_teff' + suffix + '.pdf')

plt.show()


# ------------------------------------------------------------------------------
# Ploting histogram
st = []
for i in range(len(sptyp)):
    st.append(sptyp[i][0:2])
st = np.array(st)
index = np.where([binarity == 'SB2'])
index2 = np.where([binarity == 'SB1'])
index = np.concatenate([index[1], index2[1]])
# st = st[index[1]]
st = st[index]

plt.clf()
letter_counts = Counter(st)

letter_counts = collections.OrderedDict(sorted(letter_counts.items()))
letter_counts.pop("OC")
letter_counts.pop("ON")
letter_counts.pop("--")
letter_counts.move_to_end('B0')
letter_counts.move_to_end('B1')
letter_counts.move_to_end('B2')
letter_counts.move_to_end('B3')
letter_counts.move_to_end('B4')
letter_counts.move_to_end('B5')
letter_counts.move_to_end('B6')
letter_counts.move_to_end('B7')
letter_counts.move_to_end('B8')
letter_counts.move_to_end('B9')

df = pandas.DataFrame.from_dict(letter_counts, orient='index')
df.plot(kind='bar', sort_columns=False, alpha=0.5, legend=None, width=1)
# plt.legend()
# df.hist()
plt.savefig('histogram_SB12.pdf')
plt.show()
