#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Version 1.01
# Birth date: 2014 05 28
# Last modification: 2016 12 26

# coding: utf-8
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import math
import pyhdust as phd
from pyhdust import phc
import datetime as dt
from matplotlib import rc
import os
import jdcal
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import pyhdust.spectools as spt
from scipy.optimize import curve_fit
from astropy import constants
import csv
from matplotlib.ticker import FormatStrFormatter
from pyhdust import stats
import matplotlib.gridspec as gridspec
from astropy.wcs import WCS
from astropy.io import fits
import fitsio
from fnmatch import fnmatch
from pyhdust import spectools

rc('font', size=12)

# ==============================================================================
# Define these parameters
# c = phc.c.SI * 1e-3  # km/s
sigC = .9
# pc = 3.086e18
# Lsun = 3.9e33
pc = constants.pc.cgs.value
L_sun = constants.L_sun.cgs.value  # erg/ss
c = constants.c.si.value
dist = 82 * pc
Lum = 1179.70 * L_sun
norma = Lum / (4. * np.pi * dist**2)
plot_hdust_models = False

# file_star = 'fullsed_mod01_noCS_Be_M06.57_W0.76_b0.18_rp04.57_L02180_Ell.sed2'
file_star = 'fullsed_mod01_noCS_Be_M06.85_W0.96_b0.18_rp03.81_L01941_Ell.sed2'
computer = 'tangodown'
common_folder_data = '/Users/' + computer +\
    '/Dropbox/2_Artigos/tex_aara/photometry/'
# model_nodisk = common_folder + 'hdust_models/' + 'fullsed_mod01a-phot.sed2'
# model_disk = common_folder + 'hdust_models/' + 'fullsed_mod01a.sed2'

model_nodisk = common_folder_data + 'hdust_models/' + file_star
# model_disk = common_folder + 'hdust_models/' + file_disk
model_disk = common_folder_data + 'hdust_models/'
folder_tab = common_folder_data + 'tables/'


# ==============================================================================
# Basic Parameters
common_folder = '/Users/tangodown/Dropbox/2_Artigos/tex_aara/spectroscopy/'
folder_results = common_folder + 'results/'
folder_data_heros_blue = common_folder + 'fits_data/Heros/BLUE/'
folder_data_heros_red = common_folder + 'fits_data/Heros/RED/'
bess_table = common_folder + 'fits_data/bess/'
uves_data = common_folder + 'fits_data/uves_reduced/'
feros_data = common_folder + 'fits_data/feros/'
harps_data = common_folder + 'fits_data/harps/'
opd_data = common_folder + 'fits_data/opd_data/'

folder_table = common_folder + 'tables/'

if os.path.exists(folder_data_heros_red + '.DS_Store'):
    os.system('rm ' + folder_data_heros_red + '.DS_Store')

if os.path.exists(folder_data_heros_blue + '.DS_Store'):
    os.system('rm ' + folder_data_heros_blue + '.DS_Store')

if os.path.exists(bess_table + '.DS_Store'):
    os.system('rm ' + bess_table + '.DS_Store')

if os.path.exists(uves_data + '.DS_Store'):
    os.system('rm ' + uves_data + '.DS_Store')

if os.path.exists(harps_data + '.DS_Store'):
    os.system('rm ' + harps_data + '.DS_Store')

if os.path.exists(opd_data + '.DS_Store'):
    os.system('rm ' + opd_data + '.DS_Store')

colors = ['Black', 'Blue', 'Green', 'red', 'orange', 'brown', 'purple',
          'gray', 'dodgerblue', 'lightgreen', 'tomato', 'yellow', 'peru',
          'MediumVioletRed', 'LightSteelBlue', 'cyan', 'darkred', 'olive']

specparams = ['MJD', u'Equiv. Width ($\\AA$)', u'E/C',
              'Peak Separation (km/s)', 'Depth Central Reversal', 'V/R']
specparamsF = ['MJD', 'EqWidth', 'EC', 'PkSep', 'DpCtRev', 'VR']
npar = len(specparams)


# ==============================================================================
# Creating folders
if os.path.exists(folder_results) is False:
    os.system('mkdir ' + folder_results)

if os.path.exists(folder_table) is False:
    os.system('mkdir ' + folder_table)


# ==============================================================================
def create_txt_file(x, y, file_name):
    '''
    Create a txt file.

    :param x: array with n elements (array)
    :param y: array with n elements (array)
    :param file_name: file's name (string)
    :return: txt file
    '''

    writer = open(file_name, 'w')
    writer.write('#MJD, Pol, Error\n')
    writer.close()

    with open(file_name, 'a') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(x, y))

    return


# ==============================================================================
def near_pos(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx


# ==============================================================================
# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A * np.exp(-(x - mu)**2 / (2. * sigma**2))


# ==============================================================================
def date_to_jd(year, month, day):
    """
    jdutil.py ; day == FLOAT!!!
    """
    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month

    # this checks where we are in relation to October 15, 1582, the beginning
    # of the Gregorian calendar.
    if ((year < 1582) or
        (year == 1582 and month < 10) or
       (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)

    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)

    D = math.trunc(30.6001 * (monthp + 1))

    jd = B + C + D + day + 1720994.5

    return jd


# ==============================================================================
def hmsm_to_days(hour=0, min=0, sec=0, micro=0):
    """
    Convert hours, minutes, seconds, and microseconds to fractional days.
    """
    days = sec + (micro / 1.e6)
    days = min + (days / 60.)
    days = hour + (days / 60.)
    return days / 24.


# ==============================================================================
def jd_to_mjd(jd):
    return jd - 2400000.5


# ==============================================================================
def getPksp(x, data):
    i = near_pos(x, linc)
    vmax = (data[:i]).argmax()
    vmax = x[vmax]
    rmax = (data[i:]).argmax()
    rmax = x[rmax + i]
    return (rmax - vmax) / linc * c, rmax, vmax


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
def getEWVR(x, data, clevel):
    ind = np.where((x > lin0) & (x < linc))
    datai = data[ind]
    xi = x[ind]
    V = (clevel - np.sum(datai) / len(datai)) * (xi[-1] - xi[0])
    ind = np.where((x > linc) & (x < lin1))
    datai = data[ind]
    xi = x[ind]
    R = (clevel - np.sum(datai) / len(datai)) * (xi[-1] - xi[0])
    return V + R, V / R


# ==============================================================================
def getMJD(hdr):
    if ('MJD-OBS' in hdr) is False:
        date = hdr['DATE-OBS']
        if date[2] == '/':
            date = date[6:] + '-' + date[3: 5] + '-' + date[:2]
        if len(date) == 19:
            fmt = '%Y-%m-%dT%H:%M:%S'
        elif len(date) > 19:
            fmt = '%Y-%m-%dT%H:%M:%S.%f'
        else:
            if 'UT' in hdr:
                hour = hdr['UT']
            else:
                hour = hdr['UT-START']
            if hour.find(' ') > -1:
                fmt = '%Y-%m-%d' + 'T %H:%M:%S.%f'
            else:
                fmt = '%Y-%m-%d' + 'T%H:%M:%S.%f'
            date = date + 'T' + hour
        mjd = dt.datetime.strptime(date, fmt)
        mjd = date_to_jd(mjd.year, mjd.month, mjd.day +
                         hmsm_to_days(mjd.hour, mjd.minute,
                                      mjd.second, mjd.microsecond))
        mjd = jd_to_mjd(mjd)
    else:
        mjd = hdr['MJD-OBS']
    return mjd


# ==============================================================================
def getSpecParams(lbc, fitsname, type_data):
    try:
        if type_data == 'HEROS':
            hdr = pyfits.getheader(fitsname)
            data0 = pyfits.getdata(fitsname)
            x0 = hdr['CDELT1'] * (np.arange(len(data0)) -
                                  hdr['CRPIX1']) + hdr['CRVAL1']
            mjd = getMJD(hdr)
            ind = np.where((x0 < cont1) & (x0 > cont0))
            data = data0[ind]
            x = x0[ind]
        elif type_data == 'BESS':
            x, data, mjd = read_fits(fitsname)
            x = 1e4 * x
        elif type_data == 'PUCHEROS':
            x, data, mjd = read_pucheros(fitsname)
            x = 1e4 * x
        elif type_data == 'UVES':
            x, data, mjd = read_uves(fitsname)
        elif type_data == 'FEROS':
            x, data, mjd = read_feros(fitsname)
        elif type_data == 'HARPS':
            x, data, mjd = read_harps(fitsname)
        elif type_data == 'OPD':
            x, data, mjd = read_opd(fitsname)

        # Calc velocitires
        print(np.mean(x), np.mean(lbc))
        vels = wave2doppler(x, lbc)

# ------------------------------------------------------------------------------
        # check if is normalized; if not, do it:
        if np.median(data) > 2:
            data = spectools.normalize_spec(x, data)

        # clevel = (getCont(x, data, sigC, cont0) +
        #           getCont(x, data, sigC, cont1)) / 2.

        EC, vetop = spectools.ECcalc(vels, data)

        if EC > 1.1:
            pksp, rmax, vmax = spectools.PScalc(vels, data)
        else:
            pksp, rmax, vmax = np.NaN

        # pksp, rmax, vmax = getPksp(x, data)

        depth = np.max(data) - np.min(data)

        # EW = spectools.EWcalc(vels, data, clevel)
        # VR = spectools.VRcalc(vels, data)
        VR, EW = spectools.VREWcalc(vels, data)

        return mjd, EW, EC, pksp, depth, VR
    except:
        return


# ==============================================================================
def create_list_files(list_name, folder, folder_table):
    '''
    Creates a list of the files inside a given folder.

    :param list_name: list's name (string)
    :param folder: files' folder (string)
    :return: creates a txt file, with the files' paths
    '''

    a = open(folder_table + list_name + ".txt", "w")
    for path, subdirs, files in os.walk(folder):
        for filename in files:
            f = os.path.join(path, filename)
            a.write(str(f) + os.linesep)
    return


# ==============================================================================
def read_list_files_all(table_name, folder_table):
    '''
    Read list of files in a table, and returns all fits file in an array.

    :param folder_table: table's folder (string)
    :param table_name: Table's name (string)
    :return: list of files (txt file)
    '''

    file_data = folder_table + table_name

    files = open(file_data, 'r')
    lines = files.readlines()

    list_files = []
    for i in range(len(lines)):
        list_files.append(lines[i][:-1])

    files.close()

    return list_files


# ==============================================================================
def create_opd_list(root, folder):
    pattern = "*halpha.fits"
    files_list = []

    for path, subdirs, files in os.walk(root):
        for name in files:
            if fnmatch(name, pattern):
                print(os.path.join(path, name))
                files_list.append(os.path.join(path, name))
    return files_list


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
def splot(star, suffix, ref0, ref1, files, fmt, shift, cmap,
          fluxs_hdust, velcs_hdust):

    res = 25
    lbdas = []
    specs = []
    avgdif = []
    year = []
    month = []
    day = []
    nleg = []
    mjds = []
    global cor

    # for filei in files:
    for i in range(len(files)):
        filei = files[i]
        print(filei, fmt[i])
        if fmt[i] == 'HEROS' or fmt[i] == 'BESS':
            lbda, spec, mjd = read_fits(filei)
        elif fmt[i] == 1:
            lbda, spec, mjd = read_fits_narval(filei)
        elif fmt[i] == 2:
            lbda, spec, mjd = read_fits_narval_sextension(filei)
        elif fmt[i] == 3:
            lbda, spec, mjd = read_fits_dao(filei)
        elif fmt[i] == 4:
            lbda, spec, mjd = read_harps(filei)
        elif fmt[i] == 'PUCHEROS':
            lbda, spec, mjd = read_pucheros(filei)
        elif fmt[i] == 'UVES':
            lbda, spec, mjd = read_uves(filei)
        elif fmt[i] == 'FEROS':
            lbda, spec, mjd = read_feros(filei)
        elif fmt[i] == 'HARPS':
            lbda, spec, mjd = read_harps(filei)
        elif fmt[i] == 'OPD':
            lbda, spec, mjd = read_opd(filei)

        gcal = jdcal.jd2gcal(jdcal.MJD_0, mjd)
        year += [gcal[0]]
        month += [gcal[1]]
        day += [gcal[2]]
        nleg += [mjd]

        lbda, spec = cut_spec(lbda, spec, ref0, ref1)
        spec = list(spec)
        # print(np.mean(lbda))
        lbda = lbda + shift[i]
        lbda = list(lbda)
        # lenlb = len(lbda) / 2
        # R = (ref1 + ref0) / 2 / (lbda[lenlb + 1] - lbda[lenlb])
        # print('# Opening {}... R={}'.format(filei, R))
        lbdas += [lbda]
        print(mjd)
        mjds.append(mjd)
        specs += [spec]

    v1 = min(nleg)
    v2 = max(nleg)

    # print(len(mjds))
    # cor = phc.gradColor(nleg, cmapn=cmap, min=v1, max=v2)
    cor = phc.gradColor(mjds, cmapn=cmap, min=v1, max=v2)

    # print(np.shape(specs), np.shape(lbdas), np.shape(nleg), np.shape(cor))
    # nleg, lbdas, specs, cor = zip(*sorted(zip(nleg, lbdas, specs, cor)))
    temp = np.arange(len(cor))
    nleg, lbdas, specs, temp = zip(*sorted(zip(nleg, lbdas, specs, temp)))
    # nleg, cor, specs = zip(*sorted(zip(nleg, lbdas, specs)))
    cor2 = []
    for h in temp:
        cor2.append(cor[h])
    cor = np.copy(cor2)
    # cor = phc.gradColor(mjds, cmapn=cmap, min=v1, max=v2)
    # p = mjds.argsort()
    # nleg = nleg[p]
    # lbdas = lbdas[p]
    # specs = specs[p]
    # cor = cor[p]

    lbdref = np.linspace(ref0, ref1, res)
    ipspecs = interpolspline_specs(lbdref, lbdas, specs)
    ind = np.where((lbdref <= ref1 - .001) & (lbdref >= ref0 + .001))

    plt.clf()

    for i in range(len(files)):
        # print(np.mean(lbdas[i]))
        avgdif += [np.mean(ipspecs[i][ind] - ipspecs[0][ind])]
        plt.plot(lbdas[i], specs[i] - avgdif[i] + 0.3 * i, color=cor[i],
                 label=('%s-%s-%s' % (year[i], month[i], day[i])))
        # print('tamanho lbda', len(lbdas[i]))

    # Colour bar
    # norm = plt.Normalize(vmin=np.min(nleg), vmax=np.max(nleg))
    norm = plt.Normalize(vmin=v1, vmax=v2)
    s_m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    colourbar = plt.colorbar(s_m)
    colourbar.set_label('MJD')

    # General properties
    plt.ylim(0.0, 2.5)
    plt.xlim(ref0, ref1)
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    plt.ylabel(r'F$_\lambda$ $/$ F$_c$')
    plt.minorticks_on()
    plt.rc('font', family='serif', size=13)
    plt.autoscale(enable=True, axis='y')
    plt.savefig('results/' + star + '_{}.eps'.format(suffix), dpi=300)

    # Plotting in the velocity space
    plt.clf()
    if suffix == 'halpha':
        central_wave = 0.6563
        len_flux = 700

    if suffix == 'hbeta':
        central_wave = 0.4861
        len_flux = 600

    # print(alf)
    avg_vel = np.zeros(len_flux)
    avg_flux = np.zeros(len_flux)
    cont = 0
    for i in range(len(files)):
        res = np.abs(lbdas[i][1] - lbdas[i][0])
        # print(np.mean(lbdas[i]))
        vel, flux = fluxXvel(wave=lbdas[i], flux=specs[i],
                             central_wave=central_wave, delta=res)
        # vel, flux = spt.lineProf(lbdas[i], specs[i], lbc=central_wave)
        # print(np.mean(vel))
        if len(flux) == len_flux:
            avg_flux = avg_flux + flux
            avg_vel = avg_vel + vel
            cont = cont + 1
        flux = np.array(flux)
        plt.plot(vel, flux + np.random.random(), color='black',
                 alpha=np.random.random())
    avg_flux = avg_flux / cont
    avg_vel = avg_vel / cont
    plt.plot(avg_vel, avg_flux, '-', label='mean', color='red')

    # Plotting hdust models ****
    for i in range(len(fluxs_hdust)):
        plt.plot(velcs_hdust[i], fluxs_hdust[i])
        # print(np.mean(velcs_hdust[i]), np.mean(fluxs_hdust[i]))

    plt.legend()

    if suffix == 'halpha':
        plt.ylim(0.0, 5.5)
    if suffix == 'hbeta':
        plt.ylim(0.0, 2.5)

# ------------------------------------------------------------------------------
# ***** Plot hdust model
    if plot_hdust_models is True:
        os.system('rm ' + folder_tab + 'files_disks.txt')
        create_list_files(list_name='files_disks', folder=model_disk,
                          folder_table=folder_tab)
        file_disk = read_list_files_all(table_name='files_disks.txt',
                                        folder_table=folder_tab)

        i_angles = [90, 83.6, 77.2, 70.5, 63.6, 56.3, 48.2, 38.9, 27.3, 0.]
        chi2_red = []
        chi2red = np.zeros((10, len(file_disk)))
        for i in range(len(file_disk)):
            array = plot_hdust_model(fullsedfile=file_disk[i])
            nobs = len(array[:, :, 7])

            for j in range(nobs):
                angle = i_angles[j]
                lbdarr = array[j, :, 2] * 1.e4  # Angstrom
                flux = array[j, :, 3]  # ***
                res = np.abs(lbdarr[1] - lbdarr[0])
                vel, flux = fluxXvel(wave=lbdarr, flux=flux,
                                     central_wave=central_wave, delta=res)
                ef = avg_flux * 0.05  # ***
                for k in range(len(avg_vel)):
                    tmp = find_nearest(array=lbdarr, value=avg_vel[k])
                    index = tmp[1]
                    if angle >= 30. and angle <= 70.:
                        peso = 0.1
                    else:
                        peso = 1.0
                    # print(index, k, len(avg_flux))
                    chi2red[j][i] = (chi2red[j][i] +
                                     (avg_flux[k] - flux[index])**2 /
                                     (ef[k]) ** 2) * peso

        chi2min = np.min(chi2red)

        # encontrando o indice do minimo chi2
        index = np.argwhere(chi2red == chi2min)
        # print(index)
        index = index[0]

        # calculando o chi2 reduzido
        m = 2  # numero de graus de liberdade
        n = len(lbdarr)
        chi2_red_min = chi2min / (n - m - 1)
        chi2_red = chi2red / (n - m - 1)

        bestmodel = file_disk[index[1]]
        best_angle = index[0]

        array = plot_hdust_model(fullsedfile=bestmodel)
        lbdarr_best = array[best_angle, :, 2]
        flux_best = array[best_angle, :, 3]

        vel_best, flux_best = fluxXvel(wave=lbdarr_best, flux=flux_best,
                                       central_wave=central_wave, delta=res)

        # print(len(vel_best), len(flux_best))

# ------------------------------------------------------------------------------
        # print(min(avg_vel), max(avg_vel))
        plt.plot(vel_best, flux_best, linewidth=5)
        plt.xlim(min(avg_vel), max(avg_vel))
        plt.ylabel(r'F$_\lambda$ $/$ F$_c$')
        plt.xlabel(r'$V [km/s]$')
        plt.savefig('results/' + star + '_vel2_{}.pdf'.format(suffix), dpi=300)

    return lbdas, specs, nleg


# ==============================================================================
def flux_fit():
    # **** Voltar aqui
    # In this section, we look at a simple example of fitting a Gaussian to a
    # simulated dataset. We use the Gaussian1D and Trapezoid1D models and the
    # LevMarLSQFitter fitter to fit the data:

    import numpy as np
    from astropy.modeling import models, fitting
    import matplotlib.pyplot as plt

    # Generate fake data
    np.random.seed(0)
    x = np.linspace(-5., 5., 200)
    y = 3 * np.exp(-0.5 * (x - 1.3)**2 / 0.8**2)
    y += np.random.normal(0., 0.2, x.shape)

    # Fit the data using a box model
    t_init = models.Trapezoid1D(amplitude=1., x_0=0., width=1., slope=0.5)
    fit_t = fitting.LevMarLSQFitter()
    t = fit_t(t_init, x, y)

    # Fit the data using a Gaussian
    g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y)

    # Plot the data with the best-fit model
    plt.figure(figsize=(8, 5))
    plt.plot(x, y, 'ko')
    plt.plot(x, t(x), 'b-', lw=2, label='Trapezoid')
    plt.plot(x, g(x), 'r-', lw=2, label='Gaussian')
    plt.xlabel('Position')
    plt.ylabel('Flux')
    plt.legend(loc=2)
    # plt.show()

    return


# ==============================================================================
def wave2doppler(w, w0):
    """ Wavelength to vel., in km/s. `wl` and `lbc` units must be the same. """
    doppler = ((w - w0) / w0 * c) * 1e-5
    return doppler


# ==============================================================================
def fluxXvel(wave, flux, central_wave, delta):
    '''
    Function to plot a line profile.

    :param wave: Observed wavelenght of the line (numpy array)
    :param flux: Array with the fluxes (numpy array)
    :param flux_plus_i: Array with the fluxes add by a constant value
     (numpy array)
    :param delta: range to be considered around the lamb_obs (float)
    :param central_wave: Lab central wave (float)
    :param label: label to be plotted in the legend
    :param ax: subfigure name (ax, ay or az)
    :param fit: Do nothing (boolean)
    :return: velocity and associated flux (arrays)
    '''

    vel = []
    for i in range(len(wave)):
        # print(wave[i], central_wave, c)
        delta_lamb = wave[i] - central_wave
        vel_tmp = (c * 1e-3) * (delta_lamb / central_wave)
        vel.append(vel_tmp)  # km/s
        # print(wave[i], central_wave, c, vel_tmp)

    return vel, flux


# ==============================================================================
def calc_spec_params(lbc, fitslist, npar, colors, folder_results, type_data,
                     specparamsF, typ, wave):
    outnames = np.zeros(len(fitslist), dtype='|S128')
    outmtx = np.zeros((len(fitslist), npar))

    for i in range(len(fitslist)):
        fitsname = fitslist[i]
        # print(fitsname)
        outnames[i] = fitsname
        # print(fitsname, type_data[i])
        outmtx[i] = getSpecParams(lbc, fitsname, type_data[i])

    # Sort outmtx by JD
    idsort = outmtx[:, 0].argsort()
    outnames = outnames[idsort]
    outmtx = outmtx[idsort]

    for i in range(1, npar):
        plt.clf()
        plt.plot(outmtx[:, 0], outmtx[:, i], 'ro', color=colors[i])
        plt.xlabel(specparams[0])
        plt.ylabel(specparams[i])
        if False:
            ind = np.where((np.isnan(outmtx[:, i]) is False) &
                           (outmtx[:, 0] > 0))
            # print(outmtx[:, 0][ind])
            fit = np.polyfit(outmtx[:, 0][ind], outmtx[:, i][ind], 1)
            fit_fn = np.poly1d(fit)

            plt.title(u'H$' + wave + '$; a$_1$ = {:.2e}'.format(fit[0]))
            plt.plot(outmtx[:, 0][ind], fit_fn(outmtx[:, 0][ind]), '--k')
        else:
            plt.title(u'H$' + wave + '$')
        plt.rc('font', family='serif', size=13)
        plt.minorticks_on()
        create_txt_file(x=outmtx[:, 0], y=outmtx[:, i],
                        file_name=folder_results + typ + '_' +
                        specparamsF[i] + '.txt')

        plt.savefig(folder_results +
                    'fits_' + typ + '_{}.png'.format(specparamsF[i]))

    np.savetxt(folder_results + 'fitsnames_' + typ + '.txt',
               outnames, fmt='%s')
    np.savetxt(folder_results + 'fitsdata_' + typ + '.txt', outmtx)

    return


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
def read_fits(fname):
    # FEROS/HEROS and Bess
    # wavelength in Angstroms
    # normalized flux
    imfits = pyfits.open(fname)
    spec = imfits[0].data
    lbda = np.arange(len(spec)) * imfits[0].header['CDELT1'] +\
        imfits[0].header['CRVAL1']
    to_u = 1e-4
    hdr = pyfits.getheader(fname)
    mjd = getMJD(hdr)
    # print(mjd)
    return lbda * to_u, spec, mjd


# ==============================================================================
def read_fits_dao(fname):
    # DAO - Canadian
    # wavelength in Angstroms
    # normalized flux
    imfits = pyfits.open(fname)
    spec = imfits[0].data
    new_spec = []
    for i in range(4200 - 1):  # VOLTAR AQUI
        new_spec.extend(spec[i])
    lbda = np.arange(len(new_spec)) * imfits[0].header['DELTA_WL'] +\
        imfits[0].header['REFPIXEL']
    new_spec = np.array(new_spec)
    to_u = 1e-4
    hdr = pyfits.getheader(fname)
    mjd = getMJD(hdr)
    # print(mjd)
    return lbda * to_u, new_spec, mjd


# ==============================================================================
def read_fits_narval(fname):
    # NARVAL
    # wavelength in Angstroms
    # normalized flux
    imfits = pyfits.open(fname)
    spec = imfits[1].data.field('FLUX')
    spec_err = imfits[1].data.field('FLUX_ERR')
    lbda = hdulist2[1].data.field('AWAV')
    narval_date = imfits[0].header['TMID']
    to_u = 1e-4
    hdr = pyfits.getheader(fname)
    mjd = getMJD(hdr)
    # print(mjd)
    return lbda * to_u, spec, mjd


# ==============================================================================
def read_fits_narval_sextension(fname):
    # NARVAL with .s extension
    # wavelength in Angstroms
    # normalized flux
    imfits = np.loadtxt(fname, skiprows=2,
                        dtype=[('Wave', '|f7'), ('Flux', '|f7'),
                               ('Flux_err', '|f7')])
    spec = imfits['Flux']
    spec_err = imfits['Flux_err']
    lbda = imfits['Wave']  # Nanometer
    lbda = 10. * lbda  # Angstroms
    to_u = 1e-4
    mjd = [55873.5, 55873.5, 55873.5, 55948.5, 55948.5, 55948.5,
           55948.5, 55948.5]
    mjd = 55911.0
    return lbda * to_u, spec, mjd


# ==============================================================================
def read_harps(fname):
    to_u = 1e-4
    spec = fitsio.read(fname)
    header = fitsio.read_header(fname)
    lbda = np.arange(len(spec)) * header['CDELT1'] + header['CRVAL1']
    lbda = lbda * to_u  # to microns
    mjd = header['MJD-OBS']
    return lbda, spec, mjd


# ==============================================================================
def read_pucheros(fname):
    # Read pucheros
    # print(fname)
    data = np.loadtxt(fname, unpack=True)
    to_u = 1e-4
    lambd = data[0] * to_u  # to microns
    flux = data[1]  # normalized
    mjd = 57148.
    return lambd, flux, mjd


# ==============================================================================
def read_opd(fname):
    # file = 'aara.halpha.fits'
    spec = fitsio.read(fname)
    header = fitsio.read_header(fname)
    lbda = np.arange(len(spec)) * header['CDELT1'] + header['CRVAL1']
    jd = header['JD']
    mjd = jd - 2400000.5
    to_u = 1e-4
    lbda = to_u * lbda
    # plt.clf()
    # plt.plot(lbda, spec)
    # plt.show()
    # print(lbda, spec, mjd)
    return lbda, spec, mjd


# ==============================================================================
def read_uves(filename):
    '''Read a UVES spectrum from the ESO pipeline

    Parameters
    ----------
    filename : string
    name of the fits file with the data

    Returns
    -------
    wavelength : np.ndarray
    wavelength (in Ang)
    flux : np.ndarray
    flux (in erg/s/cm**2)
    date_obs : string
    time of observation
    '''

    try:
        suff = 'RED_SCI'
        if suff in filename:
            sp = fits.open(filename)
            header = sp[0].header
            wcs = WCS(header)
            index = np.arange(header['NAXIS1'])

            wavelength = wcs.wcs_pix2world(index[:, np.newaxis], 0)
            wavelength = wavelength.flatten()
            wavelength = wavelength * 1e-4
            flux = sp[0].data
            flux = smooth_spectrum(wavelength, flux)
            mjd = header['MJD-OBS']
            # cont = [[0.653, 0.654], [0.659, 0.660]]
            # wavelength, flux = region_around_line(wavelength, flux, cont)
            return 1e10 * wavelength, flux, mjd
    except:
        pass


# ==============================================================================
def read_feros(fname):
    # Read pucheros
    # print(fname)

    data = fitsio.read(fname)
    wave = data['WAVE'][0]
    flux = data['FLUX'][0]

    to_u = 1e-4
    lambd = wave * to_u  # to microns
    header = fitsio.read_header(fname)
    mjd = header['MJD-OBS']

    return lambd, flux, mjd


# ==============================================================================
def read_txt(table_name):
    '''
    Read a simple txt file.

    :param table_name: name of the table
    :param ncols: number of columns
    :return: x, y (arrays) or x, y, z (arrays)
    '''

    data = np.loadtxt(table_name, unpack=True)

    return data


# ==============================================================================
def normalize_range(lb, spec, a, b):
    a2 = (spec[b] - spec[a]) / (lb[b] - lb[a])
    a1 = spec[a] - a2 * lb[a]
    return spec / (a1 + a2 * lb)


# ==============================================================================
def cut_spec(lbda, spec, ref0, ref1):
    # print len(spec)
    ind = np.where((lbda > ref0) & (lbda < ref1))
    lbda = lbda[ind]
    spec = spec[ind]
    indices = np.argsort(lbda)
    lbda = lbda[indices]
    spec = spec[indices]
    ind = []
    for j in range(len(lbda) - 1):
        if lbda[j] - lbda[j + 1] >= 0:
            ind += [j]
    lbda = np.delete(lbda, ind)
    spec = np.delete(spec, ind)
    # print('# Deleted {} positions of {}'.format(len(ind), len(lbda)))
    # print len(spec)
    if max(spec) > 2.:
        spec = normalize_range(lbda, spec, 0, -1)
    return lbda, spec


# ==============================================================================
def interpolate_specs(lbd_ref, lbdas, specs):
    res = len(lbd_ref)
    ipspecs = np.zeros((len(lbdas), res))
    for i in range(len(lbdas)):
        # print('# Interpolating {}...'.format(i))
        f = interpolate.interp1d(lbdas[i], specs[i], kind='cubic')
        ipspecs[i] = f(lbd_ref)
    return ipspecs


# ==============================================================================
def interpolspline_specs(lbd_ref, lbdas, specs):
    res = len(lbd_ref)
    ipspecs = np.zeros((len(lbdas), res))
    for i in range(len(lbdas)):
        # print('# Interpolating spline {}...'.format(i))
        for j in range(len(lbdas[i]) - 1):
            if lbdas[i][j] - lbdas[i][j + 1] >= 0:
                # print j
                f = interpolate.splrep(lbdas[i], specs[i], s=0)
                ipspecs[i] = interpolate.splev(lbd_ref, f, der=0)
    return ipspecs


# ==============================================================================
def plot_hdust_model(fullsedfile):
    input_file = str(fullsedfile)
    array = phd.readfullsed2(input_file)
    return array


# ==============================================================================
def plot_3D_ribbon(x, y, z, cor, label, broken_axis):
    plt.rc('ytick', labelsize=10)
    rc('font', size=10)
    labelpad = 15

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    if broken_axis is False:
        fig = plt.figure()
        ax = Axes3D(fig)
        for i in range(len(z)):
            z_temp = np.array([z[i]] * len(x[i]))  # np.copy(z[i])
            ax.plot(xs=x[i], ys=z_temp, zs=y[i], zdir='z', c=cor[i])
            # ax.fill_between(x[i], 0, y[i], alpha=0.3)
        # ax.set_ylim(51300, 51400)

    ax.set_xlabel(r'$\lambda [\mu m]$', labelpad=labelpad)
    ax.set_zlabel(r'$F_\lambda / F_c$', labelpad=labelpad * 0.5)
    ax.set_ylabel('MJD', labelpad=labelpad)
    fig.show()
    fig.savefig('results/' + 'line_profiles_' + label + '.png')

    return


# ==============================================================================
def delta_v(vel, flux):
    '''
    Calcula o shift na velocidade, baseado no centro ajustado p gaussiana
    '''

    # from plot_corr import gauss
    # coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
    v_corte = np.zeros(2)
    i_corte = np.argsort(np.abs((flux - 1.) - (flux.max() - 1.) * .5))
    v_corte[0] = vel[i_corte[0]]

    v_aux = vel[i_corte[1]]
    count = 2
    while v_corte[0] * v_aux > 0:
        v_aux = vel[i_corte[count]]
        count += 1
    v_corte[1] = v_aux
    v_corte.sort()
    delt_v = 100.

    v_trechos = np.array([v_corte[0] - delt_v, v_corte[0],
                          v_corte[1], v_corte[1] + delt_v])

    asas = np.where((vel > v_trechos[0]) * (vel < v_trechos[1]) +
                    (vel > v_trechos[2]) * (vel < v_trechos[3]))

    # p0 is the initial guess for the fitting coefficients
    # (A, mu and sigma above)
    p0 = [flux.max() - 1., v_trechos.mean(), v_trechos[2] - v_trechos[1]]
    coeff, var_matrix = curve_fit(gauss, vel, flux - 1., p0=p0)

    return coeff[1]


# ==============================================================================
def find_nearest(array, value):
    '''
    Find the nearest value inside an array
    '''

    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


# ==============================================================================
# Creating list of files
create_list_files('list_files_red', folder=folder_data_heros_red,
                  folder_table=folder_table)
fitslist_red = read_list_files_all(table_name='list_files_red.txt',
                                   folder_table=folder_table)

create_list_files('list_files_blue', folder=folder_data_heros_blue,
                  folder_table=folder_table)
fitslist_blue = read_list_files_all(table_name='list_files_blue.txt',
                                    folder_table=folder_table)

# BESS data
create_list_files('list_files_bess', folder=bess_table,
                  folder_table=folder_table)
bess_data = read_list_files_all(table_name='list_files_bess.txt',
                                folder_table=folder_table)

# Pucheros data
pucheros_data = ['/Users/tangodown/Dropbox/2_Artigos/tex_aara/spectroscopy/' +
                 'fits_data/pucheros/' +
                 'HD158427_2015-05-06_08-34-38_final_corr_corrAIR.txt']

# UVES data
create_list_files('list_files_uves', folder=uves_data,
                  folder_table=folder_table)
uves = read_list_files_all(table_name='list_files_uves.txt',
                           folder_table=folder_table)
uves_red = []
uves_blue = []
for i in range(len(uves)):
    print(uves[i])
    try:
        wave, flux, mjd = read_uves(uves[i])
        if np.mean(wave) < 0.65:
            uves_blue.append(uves[i])
        else:
            uves_red.append(uves[i])
    except:
        pass

# FEROS data
create_list_files('list_files_feros', folder=feros_data,
                  folder_table=folder_table)
feros = read_list_files_all(table_name='list_files_feros.txt',
                            folder_table=folder_table)

# HARPS data
create_list_files('list_files_harps', folder=harps_data,
                  folder_table=folder_table)
harps = read_list_files_all(table_name='list_files_harps.txt',
                            folder_table=folder_table)

# OPD data
# create_list_files('list_files_opd', folder=opd_data,
#                   folder_table=folder_table)
# opd = read_list_files_all(table_name='list_files_opd.txt',
#                           folder_table=folder_table)
opd = create_opd_list(root=common_folder, folder=opd_data)

# ------------------------------------------------------------------------------
typ_obs_red = np.concatenate([['HEROS'] * len(fitslist_red),
                             ['BESS'] * len(bess_data),
                             ['PUCHEROS'] * len(pucheros_data),
                             ['UVES'] * len(uves_red),
                             ['FEROS'] * len(feros),
                             ['HARPS'] * len(harps),
                             ['OPD'] * len(opd)])
typ_obs_blue = np.concatenate([['HEROS'] * len(fitslist_blue),
                              ['BESS'] * len(bess_data),
                              ['PUCHEROS'] * len(pucheros_data),
                              ['UVES'] * len(uves_blue),
                              ['FEROS'] * len(feros),
                              ['HARPS'] * len(harps),
                              ['OPD'] * len(opd)])

total_data_red = np.concatenate([bess_data, fitslist_red, pucheros_data,
                                 uves_red, feros, harps, opd])
bess_blue = np.delete(bess_data, [1, 2])
total_data_blue = np.concatenate([bess_blue, fitslist_blue, pucheros_data,
                                  uves_blue, feros, harps, opd])
total_data_red = list(total_data_red)
total_data_blue = list(total_data_blue)


# ==============================================================================
# Calculating parameters

# BLUE
cont0, cont1 = (4850, 4876)
# cont0, cont1 = (0.4850, 0.4876)
dcont = 1.0
lin0, lin1 = (cont0 + dcont, cont1 - dcont)
linc = 4863
# linc = 0.4863
calc_spec_params(lbc=linc, fitslist=total_data_blue, npar=npar, colors=colors,
                 folder_results=folder_results, type_data=typ_obs_blue,
                 specparamsF=specparamsF, typ='blue', wave='\\beta')

# RED
cont0, cont1 = (6554., 6572.5)
# cont0, cont1 = (0.6554, 0.65725)
dcont = 1.0
lin0, lin1 = (cont0 + dcont, cont1 - dcont)
linc = 6562.8
# linc = 0.65628
# print(fitslist_red)
calc_spec_params(lbc=linc, fitslist=total_data_red, npar=npar, colors=colors,
                 folder_results=folder_results, type_data=typ_obs_red,
                 specparamsF=specparamsF, typ='red', wave='\\alpha')


# ==============================================================================
# Plotting line plot_line_profiles
plot_line_profiles = True
if plot_line_profiles is True:
    # Plotting Halpha profiles
    star = 'aara'
    ref0 = .6530
    ref1 = .6600
    # ref0 = .6520
    # ref1 = .6610

    suffix = 'halpha'
    cmap = 'inferno'

    # Plot
    lbd0 = 0.6563
    # c = 29979245800. / 1e5
    os.system('rm ' + folder_tab + 'files_disks.txt')
    create_list_files(list_name='files_disks', folder=model_disk,
                      folder_table=folder_tab)
    file_disk = read_list_files_all(table_name='files_disks.txt',
                                    folder_table=folder_tab)

    fluxs_hdust, velcs_hdust = [], []
    cor = phc.gradColor(np.arange(len(file_disk)), cmapn=cmap)
    for i in range(len(file_disk)):
        # print(file_disk[i])
        array = plot_hdust_model(fullsedfile=file_disk[i])
        obs = array[:, 0, 0]
        nobs = len(obs)
        for iobs in range(nobs):
            lbd_hdust = array[iobs, :, 2]
            # print(lbd_hdust)
            sed = array[iobs, :, 3]
            flux_hdust = norma * sed
            normflx = spt.linfit(lbd_hdust, flux_hdust)
            # print(np.mean(normflx))
            vel_hdust, flux_halpha = spt.lineProf(lbd_hdust, normflx, lbc=lbd0)
            # res = np.abs(lbd_hdust[1] - lbd_hdust[0])
            # vel_hdust, flux_halpha = fluxXvel(lbd_hdust, flux_hdust,
            #                                   lbd0, res)
            fluxs_hdust.append(flux_halpha)
            velcs_hdust.append(vel_hdust)
        try:
            deltav = delta_v(vel_hdust, flux_halpha)
            plt.plot(vel_hdust - deltav, flux_halpha,
                     color=cor[i])  # label=label[iobs]
        except:
            print('It was not possible for the file {}'.format(file_disk[i]))
    # plt.xlim(-400, 400)
    # plt.savefig('teste.pdf')
    print(75 * '=')
    print('\nhalfa')
    shift = [0] * len(total_data_red)

    lbdas, specs, mjds = splot(star, suffix, ref0, ref1, total_data_red,
                               typ_obs_red, shift, cmap, fluxs_hdust,
                               velcs_hdust)

    # Second option
    plot_3D_ribbon(x=lbdas, y=specs, z=mjds, cor=cor,
                   label=suffix, broken_axis=False)

# ------------------------------------------------------------------------------
    # Plotting Hbeta profiles
    print(75 * '=')
    print('\nhbeta')
    star = 'aara'
    ref0 = 0.4830
    ref1 = 0.4890
    suffix = 'hbeta'
    cmap = 'inferno'

    shift = [0] * len(total_data_blue)

    lbdas, specs, mjds = splot(star, suffix, ref0, ref1, total_data_blue,
                               typ_obs_blue, shift, cmap, fluxs_hdust,
                               velcs_hdust)
    # Second option
    plot_3D_ribbon(x=lbdas, y=specs, z=mjds, cor=cor,
                   label=suffix, broken_axis=False)


# ==============================================================================
def plot_subplot(x, y, pos, ylabel, cor, gs):

    sheet = plt.subplot(gs1[pos])
    plt.scatter(x, y, marker='o', color=cor, alpha=0.6)
    plt.setp(sheet.get_xticklabels(), visible=False)
    xlim = np.linspace(min(x) - 2, max(x) + 2, len(x))
    plt.xlim(min(x) - 2, max(x) + 2)
    plt.ylim(min(y) - 0.10, max(y) + 0.10)
    plt.yticks(np.linspace(min(y) - 0.05, max(y) + 0.05, 4), fontsize=10)
    sheet.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    avg = np.median(y)
    sigma_avg = stats.mad(y)
    plt.hlines(avg - sigma_avg, xmin=min(x) - 100,
               xmax=max(x) + 100, colors='k', linestyles='--')
    plt.hlines(avg + sigma_avg, xmin=min(x) - 100,
               xmax=max(x) + 100, colors='k', linestyles='--',
               label=ylabel + ' $\pm \, \sigma = {0:.2f} \pm {1:.2f} \%$'.
               format(avg, sigma_avg))
    plt.fill_between(xlim, (avg - sigma_avg) * np.ones(len(x)),
                     (avg + sigma_avg) * np.ones(len(x)),
                     facecolor='orange', alpha=0.3,
                     linestyles='--', color='gray')
    plt.ylabel(ylabel, fontsize=10)
    plt.legend(fontsize=10, loc=4, frameon=False)
    plt.minorticks_on()

    return


# ==============================================================================
# Secular evolution
mjd_1, pksp = read_txt(table_name=folder_results + 'red_PkSep.txt')
mjd_2, EW = read_txt(table_name=folder_results + 'red_EqWidth.txt')
mjd_3, EC = read_txt(table_name=folder_results + 'red_EC.txt')
mjd_4, VR = read_txt(table_name=folder_results + 'red_VR.txt')
mjd_5, DpCtRev = read_txt(table_name=folder_results + 'red_DpCtRev.txt')


# ------------------------------------------------------------------------------
# Plots
cmap = 'rainbow'
cor = phc.gradColor(np.arange(len(mjd_1)), cmapn=cmap)

gs1 = gridspec.GridSpec(5, 1)
gs1.update(hspace=0.00)

indx = np.where(mjd_1 < 51370)
mjd_1 = mjd_1[indx]
mjd_2 = mjd_2[indx]
mjd_3 = mjd_3[indx]
mjd_4 = mjd_4[indx]
mjd_5 = mjd_5[indx]
cor = phc.gradColor(np.arange(len(mjd_1)), cmapn=cmap)

pksp = pksp[indx]
EW = EW[indx]
EC = EC[indx]
VR = VR[indx]
DpCtRev = DpCtRev[indx]

plot_subplot(x=mjd_1, y=pksp / 1000, pos=0,
             ylabel=r'$Pk. sep. \, \mathrm{[km/s]}$',
             cor=cor, gs=gs1)
plot_subplot(x=mjd_2, y=EW, pos=1, ylabel=r'$EW \, \mathrm{[\AA]}$',
             cor=cor, gs=gs1)
plot_subplot(x=mjd_3, y=EC, pos=2, ylabel=r'$EC$', cor=cor, gs=gs1)
plot_subplot(x=mjd_4, y=VR, pos=3, ylabel=r'$V/R$', cor=cor, gs=gs1)
plot_subplot(x=mjd_5, y=DpCtRev, pos=4, ylabel=r'$DpCtRev$', cor=cor, gs=gs1)
plt.xlabel('MJD', fontsize=fontsize)
plt.yscale('linear')
plt.xscale('log')
plt.tight_layout()
plt.savefig(folder_results + 'spec_params.png')



# ----------
# plt.clf()
# f_diff = []
# for i in range(len(specs)):
#     # tmp = np.array(specs[i])
#     tmp = specs[i]
#     flux_avg = np.mean(tmp)
#     tmp2 = specs[i] - flux_avg * np.ones(len(specs[i]))
#     tmp2 = list(tmp2)
#     if len(tmp2) > 690 and len(tmp2) < 710:
#         f_diff.append(tmp2)
# plt.imshow(f_diff, aspect="auto", origin="lower")
# plt.show()


# X, Y = np.meshgrid(mjds, range(num_of_frequency_bins))

