import pyfits
import matplotlib.pylab as plt
import os
from fnmatch import fnmatch
from matplotlib import rcParams
import numpy as np
from matplotlib import rc
rc('font', size=16)
rcParams['font.family'] = 'Times'
# > linfit / linprof
# cont0, cont1 = (6554., 6572.5)
cont0, cont1 = (2.1691, 2.1717)

dcont = 0.0005
sigC = 0.9
lin0, lin1 = (cont0 + dcont, cont1 - dcont)
linc = 2.17028

# READ CRIRES
common_folder = '/Users/tangodown/Dropbox/2_Artigos/tex_aara/spectroscopy/'
folder_data = common_folder + 'fits_data/crires/'

pattern = "*.fits"

list = []

for path, subdirs, files in os.walk(folder_data):
    for name in files:
        if fnmatch(name, pattern):
            print(os.path.join(path, name))
            list.append(os.path.join(path, name))
plt.clf()


# ==============================================================================
def smooth_spectrum(wave, flux, doplot=None):
    '''
    Smooth the spectrum by convolving with a (normalized) Hann
    window of 400 points.

    :param wave: Array with the wavelenght (numpy array)
    :param flux: Array with the fluxes (numpy array)
    :param doplot: Would you like to see the plot? (boolean)
    :return: smoothed flux (array)
    '''

    kernel = np.hanning(100)
    kernel = kernel / kernel.sum()
    smoothed = np.convolve(kernel, flux, mode='SAME')

    if doplot is True:
        plt.plot(wave, flux, label='original')
        plt.plot(wave, smoothed, label='smoothed', linewidth=3, color='orange')
        plt.yscale('log')
        plt.tight_layout()
        plt.legend()

    return smoothed


# ==============================================================================
def getCont(x, data, sigC, cont):
    ind = np.where((x < cont + dcont) & (x > cont - dcont))
    dvals = data[ind]
    med = np.median(dvals)
    if len(dvals) > 2:
        ind = np.where(abs(med - dvals) < sigC * np.std(dvals))
        dvals = dvals[ind]
    return np.average(dvals)


# ==============================================================================
def doNorm(lb, spec):
    spa = getCont(lb, spec, sigC, cont0)
    spb = getCont(lb, spec, sigC, cont1)
    a2 = (spb - spa) / (cont1 - cont0)
    a1 = spa - a2 * cont0
    return spec / (a1 + a2 * lb)


# ==============================================================================
def plot_read_fits(file, cor):
    hdulist = pyfits.open(file)
    data1 = hdulist[1].data
    wave1 = 1e-3 * data1['Wavelength']
    flux1 = data1['Extracted_OPT']
    mjd1 = hdulist[0].header['MJD-OBS']

    hdulist = pyfits.open(file)
    data2 = hdulist[2].data
    wave2 = 1e-3 * data2['Wavelength']
    flux2 = data2['Extracted_OPT']
    mjd2 = hdulist[0].header['MJD-OBS']

    hdulist = pyfits.open(file)
    data3 = hdulist[3].data
    wave3 = 1e-3 * data3['Wavelength']
    flux3 = data3['Extracted_OPT']
    mjd3 = hdulist[0].header['MJD-OBS']

    hdulist = pyfits.open(file)
    data4 = hdulist[4].data
    wave4 = 1e-3 * data4['Wavelength']
    flux4 = data4['Extracted_OPT']
    mjd4 = hdulist[0].header['MJD-OBS']
    # plt.plot(x=wave1, y=flux1, color=cor, label=mjd1)
    # plt.plot(x=wave2, y=flux2, color=cor, label=mjd1)
    # plt.plot(x=wave3, y=flux3, color=cor, label=mjd1)
    # plt.plot(x=wave4, y=flux4, color=cor, label=mjd1)
    # plt.show()

    return mjd1, wave1, flux1, mjd2, wave2, flux2, mjd3, wave3, flux3,\
        mjd4, wave4, flux4

# cor = ['red', 'blue', 'green']
cor = 'blue'

# ==============================================================================
# 1
# mjd1, wave1, flux1, wave2, flux2, wave3, flux3, wave4, flux4 =\
#     plot_read_fits(file=list[0], cor=cor)

# plt.plot(wave1, flux1, label=mjd1)
# plt.plot(wave2, flux2, label=mjd1)
# plt.plot(wave3, flux3, label=mjd1)
# plt.plot(wave4, flux4, label=mjd1)

# # 2
# mjd1, wave1, flux1, wave2, flux2, wave3, flux3, wave4, flux4 =\
#     plot_read_fits(file=list[1], cor=cor)

# plt.plot(wave1, flux1, label=mjd1)
# plt.plot(wave2, flux2, label=mjd1)
# plt.plot(wave3, flux3, label=mjd1)
# plt.plot(wave4, flux4, label=mjd1)

cor = ['red', 'black', 'red']
for i in range(len(list)):
    mjd1, wave1, flux1, mjd2, wave2, flux2, mjd3, wave3, flux3,\
        mjd4, wave4, flux4 =\
        plot_read_fits(file=list[i], cor=cor)

    # plt.scatter(wave2, flux2, marker='o', color='black', alpha=0.6)

    flux1 = doNorm(wave1, flux1)
    flux2 = doNorm(wave2, flux2)
    flux3 = doNorm(wave3, flux3)

    flux1 = smooth_spectrum(wave1, flux1)
    flux2_smooth = smooth_spectrum(wave2, flux2)
    flux3 = smooth_spectrum(wave3, flux3)


    # if np.average(wave2) > 2.1691 and np.average(wave2) < 2.1717:
    # plt.plot(wave1, flux1, label=mjd1)
    if mjd2 != 56113.27761326:
        print(mjd2)
        plt.scatter(wave2, flux2, marker='o',
                    alpha=0.6, s=15, label='MJD = {:.2f}'.format(mjd2),
                    color=cor[i])
        # plt.plot(wave2, flux2_smooth, label='{:.2f}'.format(mjd2),
        #          linewidth=4, alpha=0.5, color=cor[i])
    # plt.plot(wave3, flux3, label=mjd3)
    # plt.plot(wave4, flux4, label=mjd4)

fontsize = 16
plt.ylabel(r'$F_\lambda / F_{\rm c}$', fontsize=fontsize)
plt.xlabel(r'$\lambda \, [\mu m]$', fontsize=fontsize)
# plt.xticks(np.arange(2.1691, 2.1717, 0.001))
plt.xticks(np.arange(2.162, 2.169, 0.002), (2.162, 2.164, 2.166, 2.168))
plt.legend(fontsize=12, loc=1, frameon=False, markerscale=1, numpoints=1)

plt.text(x=2.1623, y=1.9, s='FeII', fontsize=18)
plt.hlines(y=1, xmin=2.162, xmax=2.170, linestyles='--', alpha=0.5)

plt.minorticks_on()
plt.xlim(2.162, 2.169)
plt.ylim(0.8, 2)
plt.savefig('results/crires_aara.pdf')
# plt.autoscale(enable=True)
plt.show()

