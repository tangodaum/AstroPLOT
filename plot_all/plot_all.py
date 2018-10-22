#!/usr/bin/env python
'''
Plots stuff
'''

# ==============================================================================
# Import packages
import numpy as np
import xml.dom.minidom
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as pcolors
import matplotlib.cm as cmx
import jdcal
import jdutil as jd
import atpy
import pyfits
import pyhdust.spectools as spt
import pyhdust as phd
import datetime as dt
from PyAstronomy.pyasl import binningx0dt
from pyhdust import phc
from glob import glob
from astropy.time import Time
from scipy.optimize import curve_fit
from PyAstronomy import pyasl
from xml.dom.minidom import parse
import os


# ==============================================================================
# General options
show_figure = False
ymin = 0.0
ymax = 1.2
xmin = 3170
xmax = 9200
autoscale = False
legend_position = (0.25, 0.75)
star = '51oph'
suffix = 'hpol'
period = 'all'
individual = True
# fullsedfile = 'hdust_models/fullsed_mod02b.sed2'
fullsedfile = 'hdust_models/fullsed_mod01b.sed2'
folder_commom = '/Users/tangodown/Dropbox/polarization/'
folder_model = '/Users/Tangodown/Dropbox/plot_all/hdust_models/'
vo_folder = '/Users/tangodown/Dropbox/plot_all/data/'
tab_folder = folder_commom + 'tables/'
fig_folder = folder_commom + 'figures/'
table = '51_oph.txt'  # PLOT GOOD DATA
extension = '.pdf'

# ==============================================================================
# read HDUST model
# plot model
tab = phd.readfullsed2(folder_model + 'fullsed_mod01b.sed2')
lbd_hdust = tab[:, :, 2]
nlbd = len(lbd_hdust)
obs = tab[:, 0, 0]
nobs = len(obs)
incl = np.arccos(obs) * 180. / np.pi
sed = tab[:, :, 3]
pol_hdust = tab[:, :, 7] 


# physical parameters
pc = 3.086e18
Lsun = 3.9e33
plx = 12.48
dist = 124. / plx * pc
Lum = 231. * Lsun
norma = Lum / (4. * np.pi * dist**2)
# iobs = 0
# flux_hdust[i_observer, i_lambda]
flux_hdust = norma * sed[:, :] 

# avermelhar modelo
# tab = np.loadtxt('extinction_law.ascii')
# lbd_ext = tab[:, 0]
# kappa = tab[:, 1]
# k_v = 211.4
# EB_V = .12
# EB_V = .2
# A_v = 3.1 * EB_V
# A_lambda = A_v * (kappa / k_v)

# flux_new = np.zeros([nobs, nlbd])
# for iobs in xrange(nobs):
#     for ilbd in xrange(nlbd):
#         ilbd_near = np.where(np.abs(lbd_hdust[ilbd] - lbd_ext) == np.min(np.abs(lbd_hdust[ilbd] - lbd_ext)))[0][0]
#         flux_new[iobs, ilbd] = flux_hdust[iobs, ilbd] * 10.**(-.4 * A_lambda[ilbd_near])

# flux_old = flux_hdust
# flux_hdust = flux_new

# ==============================================================================
def plot_hdust_model(path, fullsedfile):
    input_file = str(path) + str(fullsedfile)
    array = phd.readfullsed2(input_file)
    return array


# ==============================================================================
def gentkdates(mjd0, mjd1, fact, step, dtstart=None):
    """ Generates round dates between > mjd0 and < mjd1 in a given step.
    Valid steps are:

        'd/D/dd/DD' for days;
        'm/M/mm/MM' for months;
        'y/Y/yy/YY/yyyy/YYYY' for years.

    dtstart (optional) is expected to be in datetime.datetime.date() format
    [i.e., datetime.date(yyyy, m, d)]

    fact must be an integer
    """

    # check sanity of dtstart
    if dtstart is None:
        dtstart = dt.datetime(*jdcal.jd2gcal(jdcal.MJD_0, mjd0)[:3]).date()
        mjdst = jdcal.gcal2jd(dtstart.year, dtstart.month, dtstart.day)[1]
    else:
        mjdst = jdcal.gcal2jd(dtstart.year, dtstart.month, dtstart.day)[1]
        if mjdst < mjd0 - 1 or mjdst > mjd1:
            print('# Warning! Invalid "dtstart". Using mjd0.')
            dtstart = dt.datetime(*jdcal.jd2gcal(jdcal.MJD_0, mjd0)[:3]).date()
    # define step 'position' and vector:
    basedata = [dtstart.year, dtstart.month, dtstart.day]
    dates = []
    mjd = mjdst
    if step.upper() in ['Y', 'YY', 'YYYY']:
        i = 0
        while mjd < mjd1 + 1:
            dates += [dt.datetime(*basedata).date()]
            basedata[i] += fact
            mjd = jdcal.gcal2jd(*basedata)[1]
    elif step.upper() in ['M', 'MM']:
        i = 1
        while mjd < mjd1 + 1:
            dates += [dt.datetime(*basedata).date()]
            basedata[i] += fact
            while basedata[i] > 12:
                basedata[0] += 1
                basedata[1] -= 12
            mjd = jdcal.gcal2jd(*basedata)[1]
    elif step.upper() in ['D', 'DD']:
        i = 2
        daysvec = np.arange(1, 29, fact)
        if basedata[i] not in daysvec:
            j = 0
            while daysvec[j + 1] < basedata[i]:
                j += 1
            daysvec += basedata[i] - daysvec[j]
            idx = np.where(daysvec < 29)
            daysvec = daysvec[idx]
        else:
            j = np.where(daysvec == basedata[i])[0]
        while mjd < mjd1 + 1:
            dates += [dt.datetime(*basedata).date()]
            j += 1
            if j == len(daysvec):
                j = 0
                basedata[1] += 1
                if basedata[1] == 13:
                    basedata[1] = 1
                    basedata[0] += 1
            basedata[i] = daysvec[j]
            mjd = jdcal.gcal2jd(*basedata)[1]
    else:
        print('# ERROR! Invalid step')
        raise SystemExit(1)
    return dates


# ==============================================================================
def getMJD(hdr):
    if ('MJD-OBS' in hdr) is False:
        date = hdr['DATE']
        if date[2] == '/':
            date = date[6:] + '-' + date[3:5] + '-' + date[:2]
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
        mjd = jd.date_to_jd(mjd.year, mjd.month, mjd.day +
                            jd.hmsm_to_days(mjd.hour, mjd.minute, mjd.second,
                                            mjd.microsecond))
        mjd = jd.jd_to_mjd(mjd)
    else:
        mjd = hdr['MJD-OBS']
    return mjd


# ==============================================================================
def read_fits(fname):
    # HPOL https://archive.stsci.edu/hpol/reading.html?print=1
    # wavelength in Angstroms

    imfits = pyfits.open(fname)
    err = imfits[1].data.field('ERROR')
    lbda = imfits[1].data.field('WAVELENGTH')
    Q = imfits[1].data.field('Q')
    U = imfits[1].data.field('U')

    err = np.reshape(err, err.shape[1])
    lbda = np.reshape(lbda, lbda.shape[1])
    Q = np.reshape(Q, Q.shape[1])
    U = np.reshape(U, U.shape[1])
    hdr = pyfits.getheader(fname)
    mjd = getMJD(hdr)

    # flux, lbdas, Q, U, err = read_lis_hpol_files(table=fname)

    # Convert Polarization Q,U to Pol,PA (degrees)
    # Temporarily renormalize values to avoid floating over/underflow

    Raddeg = .01745329

    Pol = []
    PA = []
    FF = []
    nlbda = []
    error = []
    avg_err = np.mean(err)

    for i in range(len(Q)):
        FF = np.abs(Q[i]) + np.abs(U[i])
        pol_0 = FF * np.sqrt((Q[i] / FF)**2 + (U[i] / FF)**2.)
        if fname[9:12] != 'ret':
            if FF != 0. and err[i] != 0. and lbda[i] >= 3600. and\
               lbda[i]<=9000. and err[i]<1*avg_err:
                Pol.append(FF * np.sqrt((Q[i] / FF)**2 + (U[i] / FF)**2.))
                PA.append(np.fabs(0.5 * np.arctan2(U[i], Q[i]) /\
                          Raddeg + 180.))
                nlbda.append(lbda[i])
                error.append(err[i])
        else:
            # print('estou aqui')
            if FF != 0. and err[i] != 0. and lbda[i] >= 3180. and\
               lbda[i] <= 7700. and err[i] < 1 * avg_err:
                Pol.append(FF * np.sqrt((Q[i] / FF)**2 + (U[i] / FF)**2.))
                PA.append(np.fabs(0.5 * np.arctan2(U[i], Q[i]) /\
                                  Raddeg + 180.))
                nlbda.append(lbda[i])
                error.append(err[i])
    return nlbda, Pol, PA, mjd, error

# 3180A-7700A
# ==============================================================================
def read_lis_hpol_files(table, typ):

    data = open(table)
    lines = data.readlines()

    if typ == 0:
        fluxes = lines[71:310]
        Wavelengths = lines[369:608]
        Q = lines[610:849]
        U = lines[851:1090]
        errors = lines[1092:1331]
    if typ == 1:
        fluxes = lines[63:267]
        Wavelengths = lines[322:526]
        Q = lines[528:732]
        U = lines[734:939]
        errors = lines[940:1144]

    # print(len(fluxes), len(Wavelengths), len(Q), len(U), len(errors))
    flux = []
    lbdas = []
    Q_param = []
    U_param = []
    errors_param = []
    for i in range(len(fluxes)):
        flux_line = fluxes[i][1:-1].split(' ')
        flux_line = float(flux_line[0])
        flux.append(flux_line)

        lbda_line = Wavelengths[i][1:-1].split(' ')
        lbda_line = float(lbda_line[0])
        lbdas.append(lbda_line)

        Q_line = Q[i][1:-1].split(' ')
        Q_line = float(Q_line[0])
        Q_param.append(Q_line)

        U_line = U[i][1:-1].split(' ')
        U_line = float(U_line[0])
        U_param.append(U_line)

        error_line = errors[i][1:-1].split(' ')
        error_line = float(error_line[0])
        errors_param.append(error_line)

    return flux, lbdas, Q_param, U_param, errors_param


# ==============================================================================
# To plot the second array files
def splot(files, fig_folder, fmt, shift, wave_t, Pol_t, error_t,
          lbd_t, mean_t, ef_t, scalarMap):

    lbdas = []
    pols = []
    pas = []
    year = []
    month = []
    day = []
    nleg = []
    err = []
    for filei in files:
        i = files.index(filei)
        print('Reading file: %s' % files[i])
        if fmt[i] == 0:
            nlbda, Pol, PA, mjd, error = read_fits(filei)
        elif fmt[i] == 1:
            nlbda, Pol, PA, mjd, error = read_fits(filei)

        gcal = jdcal.jd2gcal(jdcal.MJD_0, mjd)
        year += [gcal[0]]
        month += [gcal[1]]
        day += [gcal[2]]
        nleg += [mjd]
        lbdas += [nlbda]
        pols += [Pol]
        pas += [PA]
        err += [error]

    # fig = plt.figure(figsize=(8, 6))
    # ax = fig.add_subplot(111, axisbg='white')
    # ax = plt.subplot2grid((4, 4), (0, 0), colspan=3, rowspan=4)
    # plt.subplot2grid((4, 4), (0, 0), colspan=3, rowspan=4)
    colors = phc.gradColor(np.arange(len(files)), cmapn='nipy_spectral')
    for i in range(len(files)):
        lbd = np.array(lbdas[i][:])
        pol = np.array(pols[i][:])
        error = np.array(err[i][:])

        bin_y, bin_lbda = binningx0dt(x=lbd, y=pol, x0=np.min(lbdas[i]),
                                      nbins=15, yerr=error, useBinCenter=True,
                                      removeEmpty=True, removeNoError=True)

        bin_lbda = bin_y[::, 0]
        bin_pol = bin_y[::, 1]
        error = bin_y[::, 2]

        if month[i] / 10. < 1. and day[i] / 10. < 1.:
            plt.step(bin_lbda, bin_pol, color=colors[i], where='mid',
                    label=('HPOL: %s-0%s-0%s' % (year[i], month[i], day[i])),
                    linewidth=1.5)
            plt.errorbar(bin_lbda, bin_pol, yerr=error,
                        ecolor=colors[i], fmt=' ')
        if month[i] / 10. < 1. and day[i] / 10. > 1.:
            plt.step(bin_lbda, bin_pol, color=colors[i],
                    label=('HPOL: %s-0%s-%s' % (year[i], month[i], day[i])),
                    linewidth=1.5)
            plt.errorbar(bin_lbda, bin_pol, yerr=error,
                        ecolor=colors[i], fmt=' ')
        if month[i] / 10. > 1. and day[i] / 10. > 1.:
            plt.step(bin_lbda, bin_pol, color=colors[i],
                    label=('HPOL: %s-%s-%s' % (year[i], month[i], day[i])),
                    linewidth=1.5)
            plt.errorbar(bin_lbda, bin_pol, yerr=error,
                          ecolor=colors[i], fmt=' ')

    plt.errorbar(wave_t, Pol_t, error_t, marker='o', linestyle=' ', alpha=0.5,
            label='OPD-Data:\n2009/09/20 - 2014/03/09', color='blue')
    plt.errorbar(lbd_t, mean_t, ef_t, marker='o', color='red', alpha=0.5,
            linestyle=' ', label='OPD-Average Data')

    return


# ==============================================================================
# To plot individual fits
def splot_ind(file, fmt, shift, scalarMap):

    print('Reading file: %s' % file[0])

    if fmt == 0:
        nlbda, Pol, PA, mjd, error = read_fits(file[0])
    elif fmt == 1:
        nlbda, Pol, PA, mjd, error = read_fits(file[0])
    return nlbda, Pol, PA, mjd, error


# ==============================================================================
# To plot individual fits using splot_ind function
def plot_fits(file, fmt, shift, scalarMap):
    lbd, pol, PA, mjd, error = splot_ind(file, 0, 0, scalarMap)
    lbd = np.array(lbd)
    pol = np.array(pol)
    error = np.array(error)

    bin_y, bin_lbda = binningx0dt(lbd, pol, x0=np.min(lbd),
                                  nbins=20, yerr=error, useBinCenter=False)
    bin_lbda = bin_y[::, 0]
    bin_pol = bin_y[::, 1]
    error = bin_y[::, 2]

    return bin_lbda, bin_pol, error


# ==============================================================================
# SED

# plot VOSA
# dir0 = '../vosa/vosa_results_11567/'
# common_folder = '/Users/jami/Desktop/51_Oph_IAG/HDUST_res/'
# model_nodisk = common_folder + 'hdust_models/' + 'fullsed_mod01b.sed2'
# model_disk = common_folder + 'jami2/runs/hdust/51_OPH/' + 'fullsed_mod05b.sed2'
# Use a vo format xml (define here the just field "vospec_file_name")

# ------------------------------------------------------------------------------
# VOSpec
vospec_file_name = '51oph_vospec.xml'
fits_file2 = vo_folder + vospec_file_name  # VOSPEC

# To acess the field names and data
t2 = atpy.Table(fits_file2, tid=1)
t2.columns

flux2 = t2['Flux0'][:]  # erg/cm2/s/A
wave2 = t2['SpectralAxis0'][:]  # Angstrom

fig = plt.figure(figsize=(8, 6))
plt.subplot(221)
# ax = fig.add_subplot(111, axisbg='white')
# ax = plt.subplot2grid((4, 4), (0, 0), colspan=3, rowspan=4)
plt.plot(wave2, flux2, 'o', label='VOSpec', color='blue', markersize=4)


label = ['${:.0f}^{{\circ}}$'.format(incl[0]),
         '${:.0f}^{{\circ}}$'.format(incl[1]),
         '${:.0f}^{{\circ}}$'.format(incl[2]),
         '${:.0f}^{{\circ}}$'.format(incl[3]),
         '${:.0f}^{{\circ}}$'.format(incl[4])]

# plot HDUST model
linestyle = ['-', '--', ':', '-.', '-']
for iobs in range(nobs):
    plt.plot(lbd_hdust[iobs]*1e4, lbd_hdust[iobs]*flux_hdust[iobs],
             color='red', ls=linestyle[iobs], label=label[iobs])

plt.xlabel('$\lambda\, \mathrm{[\AA]}$', fontsize=16)
plt.ylabel('$\lambda F_{\lambda}\, \mathrm{[erg\, s^{-1}\, cm^{-2}]}$',
           fontsize=16)
plt.yscale('log')
plt.xlim(3000, 9000)
plt.tight_layout()
plt.legend(loc=0, fontsize=6)


# ==============================================================================
# (SED) IR

plt.subplot(222)
plt.plot(wave2, flux2, 'o', label='VOSpec', color='red', markersize=4)

#plot HDUST model
linestyle = ['-', '--', ':', '-.', '-']
label = ['${:.0f}^{{\circ}}$'.format(incl[0]),
         '${:.0f}^{{\circ}}$'.format(incl[1]),
         '${:.0f}^{{\circ}}$'.format(incl[2]),
         '${:.0f}^{{\circ}}$'.format(incl[3]),
         '${:.0f}^{{\circ}}$'.format(incl[4])]
for iobs in range(nobs):
    plt.plot(lbd_hdust[iobs], lbd_hdust[iobs]*flux_hdust[iobs],
             color='red', ls=linestyle[iobs], label=label[iobs])


plt.xlabel('$\lambda\, \mathrm{[\mu m]}$', fontsize=15)
plt.ylabel('$\lambda F_{\lambda}\, \mathrm{[erg\, s^{-1}\, cm^{-2}]}$',
           fontsize=16)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.xlim(1, 100)


# ==============================================================================
# Polarization

# plt.subplot(223)
plt.subplot(223, axisbg='white')
# plt.subplot2grid((4, 4), (0, 0), colspan=3, rowspan=4)
#	xc	bandwidth   (angstrom)			xmin	xmax
# Johnson-Cousins
#	U	3663	       650 		U	0.3338	0.3988
#	B	4361	       890		B	0.3916	0.4806
#	V	5448	       840	   	V	0.5028	0.5868
#	R	6407	       158		R	0.5617	0.7197
#	I	7980	       1540		I	0.721	0.875


# filt_arr = np.array(['u', 'b', 'v', 'r', 'i'])
# lbd_arr = np.array([3663., 4361., 5448., 6407., 7980.])
# cor = ['cyan', 'blue', 'green', 'red', 'magenta']

# # read data
# file_table = dir0 + '51oph_polarization.txt'
# data = np.loadtxt(file_table, dtype={'names': ('MJD', 'night', 'ccd', 'filt', \
#                   'calc', 'stdstars', 'dth', 'sigdth', 'P', 'Q', 'U', 'th', 'sigP',\
#                   'sigth', 'outfile', 'star', 'flag', 'tags'), 'formats': (np.float,\
#                   '|S10', '|S10', '|S10', np.float, '|S20', np.float, np.float, \
#                   np.float, np.float, np.float, np.float, np.float, np.float, \
#                   '|S20', np.int, '|S2', '|S30')})

# ordem_mjd = data['MJD'].argsort()
# MJD = data['MJD'][ordem_mjd]
# P = data['P'][ordem_mjd]
# Q = data['Q'][ordem_mjd]
# U = data['U'][ordem_mjd]
# th = data['th'][ordem_mjd]
# sigP = data['sigP'][ordem_mjd]
# sigth = data['sigth'][ordem_mjd]
# filt = data['filt'][ordem_mjd]

# for i in xrange(len(filt_arr)):
#     for j in xrange(len(P[filt == filt_arr[i]])):
#         plt.errorbar(lbd_arr[i], P[filt == filt_arr[i]][j],
#                      yerr=sigP[filt == filt_arr[i]][j], color=cor[i],
#                      ecolor=cor[i], marker='s')

# ------------------------------------------------------------------------------
# Plotting OPD Data
opd_pol = np.loadtxt(tab_folder + table,
                     dtype=[('Date', '|S7'), ('JD', '|f8'),
                            ('Filter', '|S1'), ('Q', '|f8'),
                            ('U', '|f8'), ('SIGMA', '|f8'), ('P', '|f8'),
                            ('Theta', '|f8'), ('SIGMAtheor', '|f8')])

Date = opd_pol['Date']
JD = opd_pol['JD']
Filter = opd_pol['Filter']
Q = opd_pol['Q']
U = opd_pol['U']
SIGMA = opd_pol['SIGMA']
P = opd_pol['P']
Theta = opd_pol['Theta']
SIGMAtheor = opd_pol['SIGMAtheor']

colormap = np.array(['violet', 'blue', 'green', 'red', 'orange'])

# fig = plt.figure()
# ax = fig.add_subplot(111)
file_name = '51 Oph'

filtros = np.unique(Filter)
print(filtros)
for filt in filtros:
    print('Filter: %s' % filt)
    Pol = []
    MJD = []
    error = []
    for i in range(len(Date)):
        if Filter[i] == filt:
            band = filt
            Pol.append(P[i])
            MJD.append(JD[i])
            error.append(SIGMA[i])
            if Filter[i] == filtros[0]:
                cor = colormap[0]
                lab = '$P_U$'
            elif Filter[i] == filtros[1]:
                cor = colormap[1]
                lab = '$P_B$'
            elif Filter[i] == filtros[2]:
                cor = colormap[2]
                lab = '$P_V$'
            elif Filter[i] == filtros[3]:
                cor = colormap[3]
                lab = '$P_R$'
            elif Filter[i] == filtros[4]:
                cor = colormap[4]
                lab = '$P_I$'

    MJD[:] = [x for x in MJD]
    Pol[:] = [x for x in Pol]
    error[:] = [x for x in error]

# ------------------------------------------------------------------------------
# HPOL Plot (MAST files)
# Polarized Spectrum

# Johnson-Cousins UBVRI filter curves
lbd = np.array([3656., 4353., 5477., 6349., 8797.])

# Plot P vs lambda
Pol = []
error = []
wave = []
nu = 0.
nb = 0.
nv = 0.
nr = 0.
ni = 0.
wu = 0.
wb = 0.
wv = 0.
wr = 0.
wi = 0.
eu = 0.
eb = 0.
ev = 0.
er = 0.
ei = 0.

for i in range(len(Date)):
    for filt in filtros:

        if Filter[i].decode('UTF-8') == filt.decode('UTF-8'):
                filt = filt.decode('UTF-8')
                if filt == 'U':
                        wave.append(lbd[0])
                        Pol.append(P[i])
                        error.append(SIGMA[i])
                        wu = wu + P[i]
                        nu = nu + 1.
                        eu = eu + (SIGMA[i])**2.
                if filt == 'B':
                        wave.append(lbd[1])
                        Pol.append(P[i])
                        error.append(SIGMA[i])
                        wb = wb + P[i]
                        nb = nb + 1.
                        eb = eb + (SIGMA[i])**2.
                if filt == 'V':
                        wave.append(lbd[2])
                        Pol.append(P[i])
                        error.append(SIGMA[i])
                        wv = wv + P[i]
                        nv = nv + 1.
                        ev = ev + (SIGMA[i])**2.
                if filt == 'R':
                        wave.append(lbd[3])
                        Pol.append(P[i])
                        error.append(SIGMA[i])
                        wr = wr + P[i]
                        nr = nr + 1.
                        er = er + (SIGMA[i])**2.
                if filt == 'I':
                        wave.append(lbd[4])
                        Pol.append(P[i])
                        error.append(SIGMA[i])
                        wi = wi + P[i]
                        ni = ni + 1.
                        ei = ei + (SIGMA[i])**2.

eu = np.sqrt(eu)
eb = np.sqrt(eb)
ev = np.sqrt(ev)
er = np.sqrt(er)
ei = np.sqrt(ei)

ef = np.array([eu, eb, ev, er, ei])
mean = np.array([wu / nu, wb / nb, wv / nv, wr / nr, wi / ni])

# ------------------------------------------------------------------------------
# HPOL Plot (MAST files)
if period == 'all':
    files = ["hpolccd_51-oph_20000716r_hw.fits",
             "hpolccd_51-oph_19950325r_hw.fits",
             "hpolccd_51-oph_20020705b_hw.fits",
             "hpolccd_51-oph_19950520b_hw.fits",
             "hpolccd_51-oph_20020705r_hw.fits",
             "hpolccd_51-oph_19950520r_hw.fits",
             "hpolret_51-oph_19910527_hw.fits",
             "hpolccd_51-oph_19950615b_hw.fits",
             "hpolret_51-oph_19920801_hw.fits",
             "hpolccd_51-oph_19950615r_hw.fits",
             "hpolret_51-oph_19940723_hw.fits",
             "hpolccd_51-oph_20000716b_hw.fits"]


else:
    files1 = ['hpolccd_bet-cmi_19950208b_hw.fits',
              'hpolccd_bet-cmi_19950208r_hw.fits']

path = len(files) * ['hpol/']
files = [x + y for x, y in zip(path, files)]

fmt = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # 2 para arquivo .s do
shift = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# DO FIGURE
# fig, axes = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True)
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, axisbg='white')
# ax.set_xlim(3000., 9000.)
# ax.set_ylabel('P $[\%]$')
# ax.set_xlabel('$\lambda$ [$\AA$]')

hot = plt.get_cmap('hot')
cNorm = pcolors.Normalize(vmin=0, vmax=(len(files) + 1))
jet = cm = plt.get_cmap('gist_rainbow')
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

splot(files=files, fig_folder=fig_folder, fmt=fmt, shift=shift,
      wave_t=wave, Pol_t=Pol, error_t=error, lbd_t=lbd, mean_t=mean,
      ef_t=ef, scalarMap=scalarMap)


# ------------------------------------------------------------------------------
# plot HDUST model
linestyle = ['-', '--', ':', '-.', '-']
label = ['${:.0f}^{{\circ}}$'.format(incl[0]),
         '${:.0f}^{{\circ}}$'.format(incl[1]),
         '${:.0f}^{{\circ}}$'.format(incl[2]),
         '${:.0f}^{{\circ}}$'.format(incl[3]),
         '${:.0f}^{{\circ}}$'.format(incl[4])]
for iobs in range(nobs):
    plt.plot(lbd_hdust[iobs] * 1e4, 1e2 * pol_hdust[iobs],
             color='red', ls=linestyle[iobs], label=label[iobs])

plt.xlabel('$\lambda\, \mathrm{[\AA]}$', fontsize=16)
plt.ylabel('$P\, \mathrm{[\%]}$', fontsize=16)
plt.xlim(3000, 9000)
plt.tight_layout()
# plt.legend(loc=0)

# ==============================================================================
# Halpha

lbd0 = 656.28
c = 29979245800. / 1e5

plt.subplot(224)


# ------------------------------------------------------------------------------
# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A * np.exp(-(x-mu)**2 / (2. * sigma**2))


# ------------------------------------------------------------------------------
def delta_v(vel, flux):
    '''
    Calcula o shift na velocidade, baseado no centro ajustado p gaussiana
    '''

    # from plot_corr import gauss
    #coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
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


# ------------------------------------------------------------------------------
# plot ESPaDOnS
lines = glob('../data/51 Oph 20150917.fits')

for n in range(len(lines)):
    fname = lines[n]
    # read fits
    hdr_list = pyfits.open(fname)
    fits_data = hdr_list[0].data
    fits_header = hdr_list[0].header
    lbd = fits_data[0, :]
    ordem = lbd.argsort()
    lbd = lbd[ordem]
    flux_norm = fits_data[1, ordem]   
    vel, flux = spt.lineProf(lbd, flux_norm, lbc=lbd0)
    deltav = delta_v(vel, flux) 
    plt.plot(vel - deltav, flux, 'k-')


# ==============================================================================
# plot HDUST model
for iobs in range(nobs):
    vel_hdust, flux_halpha = spt.lineProf(lbd_hdust[iobs], flux_hdust[iobs],
                                          lbc=lbd0 * 1e-3)
    deltav = delta_v(vel_hdust, flux_halpha)
    plt.plot(vel_hdust - deltav, flux_halpha,
             color='red', ls=linestyle[iobs])  # label=label[iobs]

plt.xlabel('$velocity\, \mathrm{[km\, s^{-1}]}$', fontsize=16)
plt.ylabel('$\mathrm{Intensity}$', fontsize=16)
plt.xlim(-650, 650)
plt.tight_layout()
plt.savefig('plot_all.png', dpi=100)
plt.savefig('plot_all.eps', format='eps', dpi=100)


plt.figure(2)
if show_figure is True:
    plt.show()

plt.close()

# ------------------------------------------------------------------------------

