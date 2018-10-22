import matplotlib.pyplot as plt
import numpy as np
import pyhdust.phc as phc
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from matplotlib import rc
from PyAstronomy import pyasl
from numpy import sqrt, pi, exp, linspace, loadtxt
from lmfit.models import LinearModel, LorentzianModel
from lmfit import Model

rc('font', size=16)

read_again_from_simbad = False  # Do you whish to read everything again :-(?


# ==============================================================================
def gaussian(x, amp, cen, wid):
    "1-d gaussian: gaussian(x, amp, cen, wid)"
    return (amp / (sqrt(2 * pi) * wid)) * exp(-(x - cen)**2 / (2 * wid**2))


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
cat = 'J/A+A/583/A115/catalog'  # Huang 2010
catalogs = Vizier.get_catalogs(cat)
# catalog1 = catalogs[0]
# catalog2 = catalogs[1]
catalog1 = catalogs[0]

# Operating with the data
data = catalog1.as_array()

# Print available data
# print(data.dtype)

# Taking data
sptyp = data['Class'].data
name = data['Name'].data
bzh = data['__Bz__H_'].data
bzm = data['__Bz__m_'].data
e_bzm = data['e__Bz__m_'].data
e_bzh = data['e__Bz__H_'].data
e_bz = data['e__Bz__T_'].data
bz = data['__Bz__T_'].data
nz = data['__Nz__T_'].data

bz = bz / 1000.

name = name.astype('str')
sptyp = sptyp.astype('str')
bz = np.abs(bz)

# Saving vsini
if read_again_from_simbad is True:
    vsini = []
    for i in range(len(name)):
        vsini_tmp = vsini_from_simbad(name[i])
        vsini.append(vsini_tmp)
else:
    vsini = np.loadtxt('tables/table_vsini.csv')

# Plotting
cmap = 'rainbow'
count = np.array([0, 1, 2, 3, 4, 5, 6])
cor = phc.gradColor(count, cmapn=cmap, min=0, max=6)

# Plotting
# O, B, A, F, G, K, M
plt.clf()

arr_teff = []
ylimt = []
arr_bz = []
suffix = 'MS'  # 'PM', 'MS', 'GS', 'SG', 'SD', 'WD', 'CP'
for i in range(len(sptyp)):
    # print(sptyp[i][3:5], i)
    st = sptyp[i][3:5]
    if 'O' in str(sptyp[i][3:5]):
        cor_temp = cor[0]
    elif 'B' in str(sptyp[i][3:5]):
        cor_temp = cor[1]
    elif 'A' in str(sptyp[i][3:5]):
        cor_temp = cor[2]
    elif 'F' in str(sptyp[i][3:5]):
        cor_temp = cor[3]
    elif 'G' in str(sptyp[i][3:5]):
        cor_temp = cor[4]
    elif 'K' in str(sptyp[i][3:5]):
        cor_temp = cor[5]
    elif 'M' in str(sptyp[i][3:5]):
        cor_temp = cor[6]

    lc = str(sptyp[i][0:2])
    if 'PM' in lc:
        marker = 's'
        label = 'PM'
        lum_class = 'V'
    elif 'MS' in lc:
        marker = 'o'
        label = 'MS'
        lum_class = 'V'
    elif 'GS' in lc:
        marker = 'v'
        label = 'GS'
        lum_class = 'III'
    elif 'SG' in lc:
        marker = '^'
        label = 'SG'
        lum_class = 'Ia'
    elif 'SD' in lc:
        # print('aqui')
        marker = '>'
        label = 'SD'
        lum_class = 'V'
    elif 'WD' in lc:
        marker = 'd'
        label = 'WD'
        lum_class = 'V'
    elif 'CP' in lc:
        marker = '<'
        label = 'CP'
        lum_class = 'V'

    try:
        sdj = pyasl.SpecTypeDeJager()

        llum, lteff = sdj.lumAndTeff(st, lum_class)
        teff = 10**lteff
        arr_teff.append(teff)
        arr_bz.append(bz[i])

        # print(teff)
        # print(marker)
        # if np.isnan(vsini[i]) == False and suffix in st:
        if suffix in lc:
            # print(sptyp[i], vsini[i], sptyp[i][0:2])
            # ylimt.append(bz[i])
            # print(sptyp[i], cor_temp)
            # plt.scatter(x=vsini[i], y=bz[i], s=100, color=cor_temp,
            #               alpha=0.5, linewidths=1.5, edgecolors='black',
            #               marker=marker)
            plt.scatter(x=teff, y=bz[i], s=100, color=cor_temp, alpha=0.5,
                        linewidths=1.5, edgecolors='black', marker=marker)
    except:
        pass
for i in range(len(sptyp)):
    if 'HD158427' == name[i]:
        llum, lteff = sdj.lumAndTeff('B2', 'V')
        teff = 10**lteff
        plt.annotate(s=r'$\alpha$ Arae', xy=(teff, bz[i]),
                     xycoords='data', xytext=(teff + 8000, bz[i] + 1.2),
                     arrowprops=dict(arrowstyle="<-",
                                     connectionstyle="angle, angleA=200," +
                                                     "angleB=270, rad=15"))

plt.ylabel(r'$|<B_{\rm Z}>| \,  [kG]$', fontsize=18)
# plt.xlabel(r'$v \sin i \, \rm [km/s] $', fontsize=18)
plt.xlabel(r'$T_{\rm eff} \, \rm [K] $', fontsize=18)

# Gaussian Fit
# gmodel = Model(gaussian)
# result = gmodel.fit(arr_bz, x=arr_teff, amp=5, cen=5, wid=1)
# plt.plot(arr_teff, result.init_fit, 'k--')
# plt.plot(arr_teff, result.best_fit, 'r-')
# print(result.fit_report())


# plt.xlim(0, np.nanmax(vsini))
# plt.yscale('log')
# plt.xscale('log')
plt.ylim(-0.5, 10)
plt.xlim(2e3, 50000)

norm = plt.Normalize(vmin=0, vmax=6)
s_m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
s_m.set_array([])
colourbar = plt.colorbar(s_m, ticks=list(count))
# colourbar.set_label('ST')
colourbar.set_ticklabels(['O', 'B', 'A', 'F', 'G', 'K', 'M'])
plt.tight_layout()
plt.minorticks_on()
plt.legend(fontsize=14, loc=4, frameon=False, scatterpoints=1)
plt.savefig('bz_VS_teff' + suffix + '.pdf')

plt.show()

# ------------------------------------------------------------------------------
# Second plot

for i in range(len(sptyp)):
    # print(sptyp[i][3:5], i)
    st = sptyp[i][3:5]
    if 'O' in str(sptyp[i][3:5]):
        cor_temp = cor[0]
    elif 'B' in str(sptyp[i][3:5]):
        cor_temp = cor[1]
    elif 'A' in str(sptyp[i][3:5]):
        cor_temp = cor[2]
    elif 'F' in str(sptyp[i][3:5]):
        cor_temp = cor[3]
    elif 'G' in str(sptyp[i][3:5]):
        cor_temp = cor[4]
    elif 'K' in str(sptyp[i][3:5]):
        cor_temp = cor[5]
    elif 'M' in str(sptyp[i][3:5]):
        cor_temp = cor[6]

    lc = str(sptyp[i][0:2])
    if 'PM' in lc:
        marker = 's'
        label = 'PM'
        lum_class = 'V'
    elif 'MS' in lc:
        marker = 'o'
        label = 'MS'
        lum_class = 'V'
    elif 'GS' in lc:
        marker = 'v'
        label = 'GS'
        lum_class = 'III'
    elif 'SG' in lc:
        marker = '^'
        label = 'SG'
        lum_class = 'Ia'
    elif 'SD' in lc:
        # print('aqui')
        marker = '>'
        label = 'SD'
        lum_class = 'V'
    elif 'WD' in lc:
        marker = 'd'
        label = 'WD'
        lum_class = 'V'
    elif 'CP' in lc:
        marker = '<'
        label = 'CP'
        lum_class = 'V'

    try:
        sdj = pyasl.SpecTypeDeJager()

        if suffix in lc:
            plt.scatter(x=vsini[i], y=bz[i], s=100, color=cor_temp, alpha=0.5,
                        linewidths=1.5, edgecolors='black', marker=marker)

    except:
        pass
for i in range(len(sptyp)):
    if 'HD158427' == name[i]:
        llum, lteff = sdj.lumAndTeff('B2', 'V')
        teff = 10**lteff
        plt.annotate(s=r'$\alpha$ Arae', xy=(teff, bz[i]),
                     xycoords='data', xytext=(teff + 8000, bz[i] + 1.2),
                     arrowprops=dict(arrowstyle="<-",
                                     connectionstyle="angle, angleA=200," +
                                                     "angleB=270, rad=15"))

plt.ylabel(r'$|<B_{\rm Z}>| \,  [kG]$', fontsize=18)
plt.xlabel(r'$v \sin i \, \rm [km/s] $', fontsize=18)

plt.ylim(-0.1, 3)
plt.xlim(-10, 400)

norm = plt.Normalize(vmin=0, vmax=6)
s_m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
s_m.set_array([])
colourbar = plt.colorbar(s_m, ticks=list(count))
# colourbar.set_label('ST')
colourbar.set_ticklabels(['O', 'B', 'A', 'F', 'G', 'K', 'M'])
plt.tight_layout()
plt.minorticks_on()
plt.legend(fontsize=14, loc=4, frameon=False, scatterpoints=1)
plt.savefig('bz_VS_vsini' + suffix + '.pdf')

plt.show()