# ==============================================================================
import atpy
import matplotlib.pylab as plt
import numpy as np
import matplotlib.gridspec as gridspec
from pyhdust import phc
from matplotlib.ticker import FormatStrFormatter

common_folder = '/Users/tangodown/Dropbox/2_Artigos/tex_aara/spectroscopy/'
folder_data = 'fits_data/iso/'
vospec_file = common_folder + folder_data + 'vospec_aara_iso_esa_archive.xml'
folder_results = common_folder + 'results/'


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
def print_data_header(vospec_file):
    t2 = atpy.Table(vospec_file, tid=1)

    t2.columns

    return t2


# ==============================================================================
def region_around_line(w, flux, cont):
    '''cut out and normalize flux around a line

    Parameters
    ----------
    w : 1 dim np.ndarray
    array of wavelengths
    flux : np.ndarray of shape (N, len(w))
    array of flux values for different spectra in the series
    cont : list of lists
    wavelengths for continuum normalization [[low1,up1],[low2, up2]]
    that described two areas on both sides of the line
    '''
    # index is true in the region where we fit the polynomial

    indcont = ((w > cont[0][0]) & (w < cont[0][1]) |
               ((w > cont[1][0]) & (w < cont[1][1])))
    # index of the region we want to return
    indrange = (w > cont[0][0]) & (w < cont[1][1])
    # make a flux array of shape
    # (number of spectra, number of points in indrange)
    f = flux[indcont]
    wfit = w[indcont]

    linecoeff = np.polyfit(wfit, f, 2)
    f = flux / np.polyval(linecoeff, w)

    return w[indrange], f


# ==============================================================================
def plot_specific_lines(wave, flux, indx):
    flx = flux[indx]
    wvg = wave[indx]
    imax = np.where(flx == np.max(flx))
    imax = imax[0][0]

    w1 = float(wvg[0])
    w2 = float(wvg[2])  # wvg[imax] - 40
    w3 = float(wvg[-2])  # float(wvg[10])wvg[imax] + 40
    w4 = float(wvg[-1])

    cont = [[w1, w2], [w3, w4]]
    wvg, flx, = region_around_line(wvg, flx, cont)
    flx = np.delete(flx, [0, len(flx) - 1])
    return flx, wvg


# ==============================================================================
def plot_subplot(x, y, pos, ylabel, cor, loc, label, gs):

    sheet = plt.subplot(gs1[pos[0], pos[1]])
    plt.scatter(x, y, marker='o', color=cor, alpha=0.6)

    smoothed = smooth_spectrum(wave=x, flux=y)
    plt.plot(x, smoothed, color='red', linewidth=4, alpha=0.5, label=label)
    plt.setp(sheet.get_xticklabels(), visible=True)
    if pos[0] == 0 and pos[1] == 0:
        plt.xlim(min(x) + 0.0045, max(x) - 0.0045)
    elif pos[0] == 1 and pos[1] == 0:
        plt.xlim(min(x) + 0.01, max(x) - 0.01)
    elif pos[0] == 0 and pos[1] == 1:
        plt.xlim(min(x) + 0.001, max(x) - 0.001)
    elif pos[0] == 1 and pos[1] == 1:
        plt.xlim(min(x) + 0.008, max(x) - 0.008)

    plt.ylim(min(y) - 0.25, max(y) + 0.25)
    plt.yticks(np.linspace(min(y) - 0.05, max(y) + 0.05, 4), fontsize=11)

    plt.hlines(1, xmin=min(x) - 100,
               xmax=max(x) + 100, colors='k', linestyles='--')
    plt.ylabel(ylabel)
    # import matplotlib
    plt.rcParams['legend.handlelength'] = 0
    plt.rcParams['legend.numpoints'] = 1
    plt.legend(fontsize=10, loc=loc, frameon=False, markerscale=0)
    plt.minorticks_on()
    plt.locator_params(tight=True, nbins=6)
    sheet.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    sheet.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    return


# ==============================================================================
# Plotting data
t = print_data_header(vospec_file=vospec_file)

# VOSpec
flux = t['Flux0'][:]  # erg/cm2/s/A
wave = t['SpectralAxis0'][:]  # Angstrom

# smoothed = smooth_spectrum(wave=wave, flux=flux)

# plt.plot(wave, wave * flux, '-')
# plt.plot(wave, wave * smoothed, '-', color='red')
# plt.xscale('linear')
# plt.yscale('log')


# ------------------------------------------------------------------------------
# limits = [[26197, 26440], [28935, 29050], [30390, 30575], [40203, 40340],
#           [40520, 40645], [46536, 46630], [74595, 74770], [75028, 75277],
#           [223333, 224998]]

# Lines:
# He I - 4713 AA
# He I - 4.651 AA
# H I - 7.460 AA
# He I - 7.5182 AA


limits = [[30270, 30500], [40000, 40400], [46467, 46797], [74260, 75386]]

indx1 = np.where((wave > limits[0][0]) & (wave < limits[0][1]))
indx2 = np.where((wave > limits[1][0]) & (wave < limits[1][1]))
indx3 = np.where((wave > limits[2][0]) & (wave < limits[2][1]))
indx4 = np.where((wave > limits[3][0]) & (wave < limits[3][1]))

flx1, wvl1 = plot_specific_lines(wave=wave, flux=flux, indx=indx1)
flx2, wvl2 = plot_specific_lines(wave=wave, flux=flux, indx=indx2)
flx3, wvl3 = plot_specific_lines(wave=wave, flux=flux, indx=indx3)
flx4, wvl4 = plot_specific_lines(wave=wave, flux=flux, indx=indx4)

cmap = 'rainbow'
cor = phc.gradColor(np.arange(9), cmapn=cmap)
cor = 'black'

gs1 = gridspec.GridSpec(2, 2)

ylabel = r'$F_\lambda / F_{\rm c}$'
plot_subplot(x=wvl1 * 1e-4, y=flx1, pos=[0, 0], ylabel=ylabel, cor=cor, loc=1,
             label=r'He I - 3.037 $\mu m$', gs=gs1)
plot_subplot(x=wvl2 * 1e-4, y=flx2, pos=[1, 0], ylabel=ylabel, cor=cor, loc=2,
             label=r'Fe II - 4.021 $\mu m$', gs=gs1)
plt.xlabel(r'$\lambda \, [\mu m]$')

plot_subplot(x=wvl3 * 1e-4, y=flx3, pos=[0, 1], ylabel=ylabel, cor=cor, loc=1,
             label=r'He I - 4.651 $\mu m$', gs=gs1)
plot_subplot(x=wvl4 * 1e-4, y=flx4, pos=[1, 1], ylabel=ylabel, cor=cor, loc=3,
             label=r'H I - 7.460 $\mu m$, He I - 7.505 $\mu m$', gs=gs1)

plt.xlabel(r'$\lambda \, [\mu m]$')
plt.yscale('linear')
# plt.xscale('log')
plt.tight_layout()

plt.savefig(folder_results + 'iso_specs.pdf')
plt.show()




