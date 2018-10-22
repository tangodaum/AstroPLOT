# ==============================================================================
# Plot HPOL data
# Version 2.0 - Bruno - 2015/03/17
# Para rodar no ipython: %run plot_pol_vs_lambda_mjd_bcmi.py

import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import pyfits
import jdcal
import jdutil as jd
import datetime as dt
from PyAstronomy.pyasl import binningx0dt
import pyhdust as phd
from pyhdust import phc
import os
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
import csv
from pyhdust import stats

# ==============================================================================
# Plot configurations
show_figure = True
ymin = 0.0
ymax = 1.2
xmin = 3170
xmax = 9200
autoscale = False
legend_position = (0.25, 0.75)
star = 'aara'  # 'aara'
suffix = 'hpol'
period = 'all'
individual = True
computer = 'tangodown'
folder_commom = '/Users/tangodown/Dropbox/2_Artigos/tex_aara/polarization/'
common_folder_2 = '/Users/' + computer +\
    '/Dropbox/2_Artigos/tex_aara/photometry/'
tab_folder = folder_commom + 'data/'
fig_folder = folder_commom + 'figures/'
model_disk = common_folder_2 + 'hdust_models/'
folder_tab = common_folder_2 + 'tables/'
folder_tab_pol = folder_commom + 'tables/'
# table = 'aara_pol_v1.txt'  # PLOT GOOD DATA
# table_csv = "28cma_iscor.csv"  # 'aara_iscor.csv'
table_csv = 'aara_iscor.csv'
extension = '.pdf'
plot_mast = False  # There are MAST data
plot_all_models = False
read_my_reduction = False
plot_adjust = False
ref_angle = 45.
plot_gregorian_date = True
fontsize = 10


# ==============================================================================
def create_txt_file(x, y, z, w, file_name):
    '''
    Create a txt file.

    :param x: array with n elements (array)
    :param y: array with n elements (array)
    :param file_name: file's name (string)
    :return: txt file
    '''

    writer = open(file_name, 'w')
    writer.write('#MJD, Pol, Error, Filter\n')
    writer.close()

    with open(file_name, 'a') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(x, y, z, w))

    return


# ==============================================================================
def find_nearest(array, value):
    '''
    Find the nearest value inside an array
    '''

    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


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
            if filename.endswith(".sed2"):
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
def plot_hdust_model(fullsedfile):
    input_file = str(fullsedfile)
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
               lbda[i] <= 9000. and err[i] < 1 * avg_err:
                Pol.append(FF * np.sqrt((Q[i] / FF)**2 + (U[i] / FF)**2.))
                PA.append(np.fabs(0.5 * np.arctan2(U[i],
                          Q[i]) / Raddeg + 180.))
                nlbda.append(lbda[i])
                error.append(err[i])
        else:
            # print('estou aqui')
            if FF != 0. and err[i] != 0. and lbda[i] >= 3180. and\
               lbda[i] <= 7700. and err[i] < 1 * avg_err:
                Pol.append(FF * np.sqrt((Q[i] / FF)**2 + (U[i] / FF)**2.))
                PA.append(np.fabs(0.5 * np.arctan2(U[i], Q[i]) /
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
def splot(files, fullsedfile, fig_folder, fmt, shift, wave_t, Pol_t, error_t,
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

    fig = plt.figure(figsize=(8, 6), dpi=120)
    ax = fig.add_subplot(111, axisbg='white')
    ax = plt.subplot2grid((4, 4), (0, 0), colspan=3, rowspan=4)
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
            ax.step(bin_lbda, bin_pol, color=colors[i], where='mid',
                    label=('HPOL: %s-0%s-0%s' % (year[i], month[i], day[i])),
                    linewidth=1.5)
            ax.errorbar(bin_lbda, bin_pol, yerr=error,
                        ecolor=colors[i], fmt=' ')
        if month[i] / 10. < 1. and day[i] / 10. > 1.:
            ax.step(bin_lbda, bin_pol, color=colors[i],
                    label=('HPOL: %s-0%s-%s' % (year[i], month[i], day[i])),
                    linewidth=1.5)
            ax.errorbar(bin_lbda, bin_pol, yerr=error,
                        ecolor=colors[i], fmt=' ')
        if month[i] / 10. > 1. and day[i] / 10. > 1.:
            ax.step(bin_lbda, bin_pol, color=colors[i],
                    label=('HPOL: %s-%s-%s' % (year[i], month[i], day[i])),
                    linewidth=1.5)
            ax.errorbar(bin_lbda, bin_pol, yerr=error,
                        ecolor=colors[i], fmt=' ')

# ------------------------------------------------------------------------------
# Plot OPD data
    ax.errorbar(wave_t, Pol_t, error_t, marker='o', linestyle=' ', alpha=0.1,
                label='OPD-Data:\n2009/09/20 - 2014/03/09', color='blue')
    ax.errorbar(lbd_t, mean_t, ef_t, marker='o', color='red', alpha=0.5,
                linestyle=' ', label='OPD-Average Data')

# ------------------------------------------------------------------------------
# Plot hdust model
    array = plot_hdust_model(fullsedfile=fullsedfile)
    lbdarr = array[:, :, 2] * 10**4  # Angstrom
    pol = array[:, :, 7]
    nobs = len(array[:, :, 7])
    label = ['70.0', '70.0', '70.0', '70.0', '70.0']
    for i in range(nobs):
        plt.plot(lbdarr[i], pol[i] * 100, '-',
                 label=('hdust-model: {}'.format(label[i])))

    plt.ylim(ymin, ymax)
    plt.xlim(xmin, xmax)
    plt.xlabel(r'$\lambda$ [$\AA$]')
    plt.ylabel('P $[\%]$')
    ax.legend(ncol=1, shadow=False, title="", fancybox=True,
              prop={'size': 7.8}, numpoints=1,
              bbox_to_anchor=(1.01, 1), loc=2)

    plt.minorticks_on()
    plt.rc('font', family='serif', size=13)

    if period == 'all':
        figure_name = star + '_all' + extension
    else:
        figure_name = star +\
            '_{}_{}-{}-{}.eps'.format(suffix, year[0], month[0], day[0])
    plt.tight_layout()
    plt.savefig(fig_folder + figure_name, dpi=600)
    # plt.show()
    plt.close()
    print('\nArquivo criado:  %s\n' % figure_name)

# ------------------------------------------------------------------------------
# Plotting  Pa vs lambda
    print(72 * '-')
    print('\nPlotting third plot\n')
    print(72 * '-')

    for i in range(len(files)):
        lbd = np.array(lbdas[i][:])
        pa = np.array(pas[i][:])
        plt.plot(lbd, pa, color=colors[i], marker='o',
                 linestyle=' ', alpha=0.5)

    plt.xlabel(r'$\lambda$ [$\AA$]')
    plt.ylabel('PA $[\deg]$')
    figure_name = star + '_pa_all' + extension
    # plt.ylim(ymin_pa, ymax_pa)
    plt.autoscale(enable=True)
    plt.xlim(xmin, xmax)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(fig_folder + figure_name, dpi=600)
    print('\nArquivo criado:  %s\n' % figure_name)
    plt.close()

    print(72 * '-')
    return lbd_t, mean_t, ef_t


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
# First plot
print(72 * '-')
print('\nPlotting first plot\n')
print(72 * '-')

# Reading opd data'
if read_my_reduction is True:
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

if read_my_reduction is False:
    # Reading opd data from csv (beacon site)
    csv_file = tab_folder + table_csv

    df = pd.read_csv(csv_file)
    JD = df['#MJD'] + 2400000
    Filter = df['filt']
    Q = df['Q']
    U = df['U']
    SIGMA = df['sigP']
    P = df['P']
    Theta = df['th']
    SIGMAtheor = df['sigth']
    flag = df['flag']

colormap = np.array(['violet', 'blue', 'green', 'red', 'orange'])

# fig = plt.figure()
fig = plt.figure(figsize=(10, 6), dpi=120)

ax = fig.add_subplot(111)
file_name = 'aara'

Pol_U, Pol_B, Pol_V, Pol_R, Pol_I = [], [], [], [], []
MJD_U, MJD_B, MJD_V, MJD_R, MJD_I = [], [], [], [], []
err_U, err_B, err_V, err_R, err_I = [], [], [], [], []
filtros = np.unique(Filter)
for filt in filtros:
    print('Filter: %s' % filt)
    Pol = []
    MJD = []
    error = []
    for i in range(len(JD)):
        if Filter[i] == filt and flag[i] is not 'W':
            band = filt
            Pol.append(P[i])
            MJD.append(JD[i] - 2400000)
            error.append(SIGMA[i])
            if Filter[i] == 'u':
                cor = colormap[0]
                lab = '$P_U$'
                Pol_U.append(P[i])
                MJD_U.append(JD[i] - 2400000)
                err_U.append(SIGMA[i])
            elif Filter[i] == 'b':
                cor = colormap[1]
                lab = '$P_B$'
                Pol_B.append(P[i])
                MJD_B.append(JD[i] - 2400000)
                err_B.append(SIGMA[i])
            elif Filter[i] == 'v' or Filter[i] == 'v2':
                cor = colormap[2]
                lab = '$P_V$'
                Pol_V.append(P[i])
                MJD_V.append(JD[i] - 2400000)
                err_V.append(SIGMA[i])
            elif Filter[i] == 'r':
                cor = colormap[3]
                lab = '$P_R$'
                Pol_R.append(P[i])
                MJD_R.append(JD[i] - 2400000)
                err_R.append(SIGMA[i])
            elif Filter[i] == 'i':
                cor = colormap[4]
                lab = '$P_I$'
                Pol_I.append(P[i])
                MJD_I.append(JD[i] - 2400000)
                err_I.append(SIGMA[i])

    MJD[:] = [x for x in MJD]
    Pol[:] = [x for x in Pol]
    error[:] = [x for x in error]
    error = np.array(error)
    # Pol = np.array(Pol) * 100  # in percent
    # ax.errorbar(MJD, Pol, error, marker='o', linestyle=' ', label=lab,
    #             color=cor, alpha=0.6)

# ------------------------------------------------------------------------------
# upper subplot
# define two subplots with different sizes
# gs1 = gridspec.GridSpec(2, 2, width_ratios=[100, 1], height_ratios=[5, 1])
MJD_adj = np.linspace(min(MJD) - 400, max(MJD) + 100, 200)
# MJD_adj = np.linspace(min(MJD) - 100, 43000, 200)
# MJD_adj_2 = np.linspace(54000, max(MJD) + 100, 200)
ymin, ymax = 0.1, 0.8

cmap = 'rainbow'
cor = phc.gradColor(np.arange(len(filtros)), cmapn=cmap)
f = plt.figure(figsize=(10, 5))
gs1 = gridspec.GridSpec(3, 2)
# gs1.update(hspace=0.00, wspace=0.025)
gs1.update(hspace=0.10, wspace=0.025, top=0.85, bottom=0.44)
# upper subplot
upperplot1 = plt.subplot(gs1[0, 0])
plt.errorbar(MJD_U, Pol_U, err_U, marker='o', linestyle=' ',
             color=cor[0], alpha=0.6)
plt.setp(upperplot1.get_xticklabels(), visible=False)
# plt.ylim(min(Pol_U) - 0.15, max(Pol_U) + 0.05)
plt.ylim(ymin, ymax)
# plt.yticks(np.linspace(min(Pol_U) - 0.05, max(Pol_U) + 0.10, 4),
#                        fontsize=fontsize)
linspace = np.array([0.20, 0.40, 0.60, 0.80])
plt.yticks(linspace, fontsize=fontsize)

linspace = [40000, 41000, 42000, 43000]
plt.xticks(linspace, fontsize=fontsize)

# plt.xticks(np.linspace(min(MJD_V) - 100, 43000, 4), fontsize=fontsize)

upperplot1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
avg_pol = np.median(Pol_U)

sigma_pol_avg = stats.mad(Pol_U)
plt.hlines(avg_pol - sigma_pol_avg, xmin=min(MJD) - 400,
           xmax=max(MJD) + 200, colors='k', linestyles='--')
plt.hlines(avg_pol + sigma_pol_avg, xmin=min(MJD) - 400,
           xmax=max(MJD) + 200, colors='k', linestyles='--',
           label='$P \pm \sigma = {0:.2f} \pm {1:.2f} \%$'.
           format(avg_pol, sigma_pol_avg))
plt.fill_between(MJD_adj, (avg_pol - sigma_pol_avg) * np.ones(len(MJD_adj)),
                 (avg_pol + sigma_pol_avg) * np.ones(len(MJD_adj)),
                 facecolor='orange', alpha=0.3, linestyles='--', color='gray')
plt.ylabel(r'$P_{\rm U}$ $[\%]$', fontsize=fontsize)
# plt.legend(fontsize=fontsize, loc=1, frameon=False)
plt.minorticks_on()

leg = plt.legend(fontsize=fontsize, loc=1, frameon=False,
                 handlelength=0, handletextpad=0)

for item in leg.legendHandles:
    item.set_visible(False)
plt.xlim(min(MJD) - 400, 43200)

# ------------------------------------------------------------------------------
# Cria datas gregorianas no eixo superior
if plot_gregorian_date is True:
    MJD_adj_2 = np.linspace(min(MJD) - 100, 43000 + 100, 200)
    # plt.xlim(min(MJD_adj_2), max(MJD_adj_2))
    year_0, month_0, day_0, frac_0 = jdcal.jd2gcal(2400000.5, min(MJD_adj_2))
    ticks = gentkdates(min(MJD) + 200, 43000, 12, 'm',
                       dtstart=dt.datetime(year_0, month_0, day_0).date())
    mjdticks = [jdcal.gcal2jd(date.year, date.month, date.day)[1]
                for date in ticks]
    ax2 = plt.twiny()
    ax2.set_xlim(min(MJD), 43000)

    ax2.set_xticks(mjdticks)
    # ax2.set_xticklabels([date.strftime("%d %b %y") for date in ticks],
    #                     fontsize=fontsize)
    ax2.set_xticklabels([date.strftime("%Y") for date in ticks],
                        fontsize=fontsize)
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)
    plt.subplots_adjust(left=0.1, right=0.95, top=0.84, bottom=0.1)

# ------------------------------------------------------------------------------
# upper subplot
upperplot1 = plt.subplot(gs1[0, 1])
plt.errorbar(MJD_U, Pol_U, err_U, marker='o', linestyle=' ',
             color=cor[0], alpha=0.6)
plt.setp(upperplot1.get_xticklabels(), visible=False)
plt.setp(upperplot1.get_yticklabels(), visible=False)

plt.xlim(54000, max(MJD) + 100)
# plt.ylim(min(Pol_U) - 0.15, max(Pol_U) + 0.05)
plt.ylim(ymin, ymax)
# plt.yticks(np.linspace(min(Pol_U) - 0.05, max(Pol_U) + 0.10, 4),
#            fontsize=fontsize)
linspace = np.array([0.20, 0.40, 0.60, 0.80])
plt.yticks(linspace, fontsize=fontsize)

# plt.xticks(np.linspace(54000, max(MJD_V) + 10, 4), fontsize=fontsize)
linspace = [54500, 55500, 56500, 57500]
plt.xticks(linspace, fontsize=fontsize)

upperplot1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
avg_pol = np.median(Pol_U)

sigma_pol_avg = stats.mad(Pol_U)
plt.hlines(avg_pol - sigma_pol_avg, xmin=min(MJD) - 100,
           xmax=max(MJD) + 100, colors='k', linestyles='--')
plt.hlines(avg_pol + sigma_pol_avg, xmin=min(MJD) - 100,
           xmax=max(MJD) + 100, colors='k', linestyles='--',
           label='$P \pm \sigma = {0:.2f} \pm {1:.2f} \%$'.
           format(avg_pol, sigma_pol_avg))
plt.fill_between(MJD_adj, (avg_pol - sigma_pol_avg) * np.ones(len(MJD_adj)),
                 (avg_pol + sigma_pol_avg) * np.ones(len(MJD_adj)),
                 facecolor='orange', alpha=0.3, linestyles='--', color='gray')
# plt.ylabel(r'$P_{\rm U}$ $[\%]$')
leg = plt.legend(fontsize=fontsize, loc=1, frameon=False,
                 handlelength=0, handletextpad=0)

for item in leg.legendHandles:
    item.set_visible(False)

plt.minorticks_on()

# ------------------------------------------------------------------------------
# Cria datas gregorianas no eixo superior
if plot_gregorian_date is True:
    MJD_adj_2 = np.linspace(54000, max(MJD) + 100, 200)
    # plt.xlim(min(MJD_adj_2), max(MJD_adj_2))
    year_0, month_0, day_0, frac_0 = jdcal.jd2gcal(2400000.5, min(MJD_adj_2))
    ticks = gentkdates(54000 + 200, max(MJD), 12, 'm',
                       dtstart=dt.datetime(year_0, month_0, day_0).date())
    mjdticks = [jdcal.gcal2jd(date.year, date.month, date.day)[1]
                for date in ticks]
    ax2 = plt.twiny()
    ax2.set_xlim(54000, max(MJD))

    ax2.set_xticks(mjdticks)
    # ax2.set_xticklabels([date.strftime("%d %b %y") for date in ticks],
    #                     fontsize=fontsize)
    ax2.set_xticklabels([date.strftime("%Y") for date in ticks],
                        fontsize=fontsize)
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)
    plt.subplots_adjust(left=0.1, right=0.95, top=0.84, bottom=0.1)

# ------------------------------------------------------------------------------
# second subplot
second = plt.subplot(gs1[1, 0])
plt.errorbar(MJD_B, Pol_B, err_B, marker='o', linestyle=' ',
             color=cor[1], alpha=0.6)
plt.setp(second.get_xticklabels(), visible=False)
plt.xlim(min(MJD) - 400, 43200)
# plt.ylim(min(Pol_B) - 0.10, max(Pol_B) + 0.10)
plt.ylim(ymin, ymax)
plt.yticks(np.linspace(min(Pol_B) - 0.05, max(Pol_B) + 0.05, 4),
           fontsize=fontsize)
linspace = np.array([0.20, 0.40, 0.60, 0.80])
plt.yticks(linspace, fontsize=fontsize)

# plt.xticks(np.linspace(min(MJD_V) - 100, 43000, 4), fontsize=fontsize)
linspace = [40000, 41000, 42000, 43000]
plt.xticks(linspace, fontsize=fontsize)

second.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
avg_pol = np.median(Pol_B)
sigma_pol_avg = stats.mad(Pol_B)
plt.hlines(avg_pol - sigma_pol_avg, xmin=min(MJD) - 400,
           xmax=max(MJD) + 100, colors='k', linestyles='--')
plt.hlines(avg_pol + sigma_pol_avg, xmin=min(MJD) - 400,
           xmax=max(MJD) + 100, colors='k', linestyles='--',
           label='$P \pm \sigma = {0:.2f} \pm {1:.2f} \%$'.
           format(avg_pol, sigma_pol_avg))
plt.fill_between(MJD_adj, (avg_pol - sigma_pol_avg) * np.ones(len(MJD_adj)),
                 (avg_pol + sigma_pol_avg) * np.ones(len(MJD_adj)),
                 facecolor='orange', alpha=0.3, linestyles='--', color='gray')
plt.ylabel(r'$P_{\rm B}$ $[\%]$', fontsize=fontsize)
leg = plt.legend(fontsize=fontsize, loc=4, frameon=False,
                 handlelength=0, handletextpad=0)

for item in leg.legendHandles:
    item.set_visible(False)

plt.minorticks_on()

second = plt.subplot(gs1[1, 1])
plt.errorbar(MJD_B, Pol_B, err_B, marker='o', linestyle=' ',
             color=cor[1], alpha=0.6)
plt.setp(second.get_xticklabels(), visible=False)
plt.setp(second.get_yticklabels(), visible=False)

plt.xlim(54000, max(MJD) + 100)
# plt.ylim(min(Pol_B) - 0.10, max(Pol_B) + 0.10)
plt.ylim(ymin, ymax)
# plt.yticks(np.linspace(min(Pol_B) - 0.05, max(Pol_B) + 0.05, 4),
#            fontsize=fontsize)
linspace = np.array([0.20, 0.40, 0.60, 0.80])
plt.yticks(linspace, fontsize=fontsize)

# plt.xticks(np.linspace(54000 - 100, max(MJD_V) + 100, 4), fontsize=fontsize)
linspace = [54500, 55500, 56500, 57500]
plt.xticks(linspace, fontsize=fontsize)

second.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
avg_pol = np.median(Pol_B)
sigma_pol_avg = stats.mad(Pol_B)
plt.hlines(avg_pol - sigma_pol_avg, xmin=min(MJD) - 200,
           xmax=max(MJD) + 100, colors='k', linestyles='--')
plt.hlines(avg_pol + sigma_pol_avg, xmin=min(MJD) - 200,
           xmax=max(MJD) + 100, colors='k', linestyles='--',
           label='$P \pm \sigma = {0:.2f} \pm {1:.2f} \%$'.
           format(avg_pol, sigma_pol_avg))
plt.fill_between(MJD_adj, (avg_pol - sigma_pol_avg) * np.ones(len(MJD_adj)),
                 (avg_pol + sigma_pol_avg) * np.ones(len(MJD_adj)),
                 facecolor='orange', alpha=0.3, linestyles='--', color='gray')
# plt.ylabel(r'$P_{\rm B}$ $[\%]$')
leg = plt.legend(fontsize=fontsize, loc=4, frameon=False,
                 handlelength=0, handletextpad=0)

for item in leg.legendHandles:
    item.set_visible(False)

plt.minorticks_on()

# ------------------------------------------------------------------------------
# third subplot
third = plt.subplot(gs1[2, 0])
plt.errorbar(MJD_V, Pol_V, err_V, marker='o', linestyle=' ',
             color=cor[2], alpha=0.6)
# plt.setp(third.get_xticklabels(), visible=False)
plt.xlim(min(MJD) - 400, 43200)
# plt.ylim(min(Pol_V) - 0.10, max(Pol_V) + 0.10)
plt.ylim(ymin, ymax)
# plt.yticks(np.linspace(min(Pol_V) - 0.05, max(Pol_V) + 0.05, 4),
#            fontsize=fontsize)
linspace = np.array([0.20, 0.40, 0.60, 0.80])
plt.yticks(linspace, fontsize=fontsize)

# plt.xticks(np.linspace(min(MJD_V) + 200, 43000 - 200, 4), fontsize=fontsize)
linspace = [40000, 41000, 42000, 43000]
plt.xticks(linspace, fontsize=fontsize)

third.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
avg_pol = np.median(Pol_V)
sigma_pol_avg = stats.mad(Pol_V)
plt.hlines(avg_pol - sigma_pol_avg, xmin=min(MJD) - 400,
           xmax=max(MJD) + 100, colors='k', linestyles='--')
plt.hlines(avg_pol + sigma_pol_avg, xmin=min(MJD) - 400,
           xmax=max(MJD) + 100, colors='k', linestyles='--',
           label='$P \pm \sigma = {0:.2f} \pm {1:.2f} \%$'.
           format(avg_pol, sigma_pol_avg))
plt.fill_between(MJD_adj, (avg_pol - sigma_pol_avg) * np.ones(len(MJD_adj)),
                 (avg_pol + sigma_pol_avg) * np.ones(len(MJD_adj)),
                 facecolor='orange', alpha=0.3, linestyles='--', color='gray')
plt.ylabel(r'$P_{\rm V}$ $[\%]$', fontsize=fontsize)
# plt.legend(fontsize=fontsize, loc=4, frameon=False)
leg = plt.legend(fontsize=fontsize, loc=4, frameon=False,
                 handlelength=0, handletextpad=0)

for item in leg.legendHandles:
    item.set_visible(False)

plt.minorticks_on()
plt.xlabel('MJD', fontsize=fontsize)

third = plt.subplot(gs1[2, 1])
plt.errorbar(MJD_V, Pol_V, err_V, marker='o', linestyle=' ',
             color=cor[2], alpha=0.6)
# plt.setp(third.get_xticklabels(), visible=False)
plt.xlim(54000, max(MJD) + 100)
# plt.ylim(min(Pol_V) - 0.10, max(Pol_V) + 0.10)
plt.ylim(ymin, ymax)
# plt.yticks(np.linspace(min(Pol_V) - 0.05, max(Pol_V) + 0.05, 4),
#            fontsize=fontsize)
linspace = np.array([0.20, 0.40, 0.60, 0.80])
plt.yticks(linspace, fontsize=fontsize)

linspace = [54500, 55500, 56500, 57500]
plt.xticks(linspace, fontsize=fontsize)
# plt.xticks(np.linspace(54000 + 300, max(MJD_V) - 200, 4), fontsize=fontsize)

plt.setp(third.get_yticklabels(), visible=False)

third.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
avg_pol = np.median(Pol_V)
sigma_pol_avg = stats.mad(Pol_V)
plt.hlines(avg_pol - sigma_pol_avg, xmin=min(MJD) - 100,
           xmax=max(MJD) + 100, colors='k', linestyles='--')
plt.hlines(avg_pol + sigma_pol_avg, xmin=min(MJD) - 100,
           xmax=max(MJD) + 100, colors='k', linestyles='--',
           label='$P \pm \sigma = {0:.2f} \pm {1:.2f} \%$'.
           format(avg_pol, sigma_pol_avg))
plt.fill_between(MJD_adj, (avg_pol - sigma_pol_avg) * np.ones(len(MJD_adj)),
                 (avg_pol + sigma_pol_avg) * np.ones(len(MJD_adj)),
                 facecolor='orange', alpha=0.3, linestyles='--', color='gray')
# plt.ylabel(r'$P_{\rm V}$ $[\%]$')
leg = plt.legend(fontsize=fontsize, loc=4, frameon=False,
                 handlelength=0, handletextpad=0)

for item in leg.legendHandles:
    item.set_visible(False)

plt.minorticks_on()
plt.xlabel('MJD', fontsize=fontsize)

# ------------------------------------------------------------------------------
# fourth subplot
# gs1.update(hspace=0.1)
gs2 = gridspec.GridSpec(2, 2)
gs2.update(hspace=0.10, top=0.35, bottom=0.10)

fourth = plt.subplot(gs2[0, :])
# gs1.update(hspace=0.80)

plt.errorbar(MJD_R, Pol_R, err_R, marker='o', linestyle=' ',
             color=cor[3], alpha=0.6)
plt.setp(fourth.get_xticklabels(), visible=False)
# plt.xlim(55000, max(MJD) + 100)
plt.xlim(54900, max(MJD) + 100)
# plt.ylim(min(Pol_R) - 0.10, max(Pol_R) + 0.10)
plt.ylim(ymin, ymax)
# plt.yticks(np.linspace(min(Pol_R) - 0.05, max(Pol_R) + 0.05, 4),
#            fontsize=fontsize)
linspace = np.array([0.20, 0.40, 0.60, 0.80])
plt.yticks(linspace, fontsize=fontsize)

fourth.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
avg_pol = np.median(Pol_R)
sigma_pol_avg = stats.mad(Pol_R)
plt.hlines(avg_pol - sigma_pol_avg, xmin=min(MJD) - 100,
           xmax=max(MJD) + 100, colors='k', linestyles='--')
plt.hlines(avg_pol + sigma_pol_avg, xmin=min(MJD) - 100,
           xmax=max(MJD) + 100, colors='k', linestyles='--',
           label='$P \pm \sigma = {0:.2f} \pm {1:.2f} \%$'.
           format(avg_pol, sigma_pol_avg))
plt.fill_between(MJD_adj, (avg_pol - sigma_pol_avg) * np.ones(len(MJD_adj)),
                 (avg_pol + sigma_pol_avg) * np.ones(len(MJD_adj)),
                 facecolor='orange', alpha=0.3, linestyles='--', color='gray')
plt.ylabel(r'$P_{\rm R}$ $[\%]$', fontsize=fontsize)
leg = plt.legend(fontsize=fontsize, loc=4, frameon=False,
                 handlelength=0, handletextpad=0)

for item in leg.legendHandles:
    item.set_visible(False)

plt.minorticks_on()

# ------------------------------------------------------------------------------
# fifth subplot
fifth = plt.subplot(gs2[1, 0:])
plt.errorbar(MJD_I, Pol_I, err_I, marker='o', linestyle=' ',
             color=cor[4], alpha=0.6)
plt.xlim(54900, max(MJD) + 100)
# plt.ylim(min(Pol_I) - 0.10, max(Pol_I) + 0.10)
plt.ylim(ymin, ymax)
# plt.yticks(np.linspace(min(Pol_I) - 0.05, max(Pol_I) + 0.05, 4),
#            fontsize=fontsize)
linspace = np.array([0.20, 0.40, 0.60, 0.80])
plt.yticks(linspace, fontsize=fontsize)

# plt.xticks(np.linspace(55000 + 100, max(MJD_V) - 100, 6), fontsize=fontsize)
linspace = [55000, 55500, 56000, 56500, 57000, 57500]
plt.xticks(linspace, fontsize=fontsize)

fifth.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
avg_pol = np.median(Pol_I)
sigma_pol_avg = stats.mad(Pol_I)
plt.hlines(avg_pol - sigma_pol_avg, xmin=min(MJD) - 100,
           xmax=max(MJD) + 100, colors='k', linestyles='--')
plt.hlines(avg_pol + sigma_pol_avg, xmin=min(MJD) - 100,
           xmax=max(MJD) + 100, colors='k', linestyles='--',
           label='$P \pm \sigma = {0:.2f} \pm {1:.2f} \%$'.
           format(avg_pol, sigma_pol_avg))
plt.fill_between(MJD_adj, (avg_pol - sigma_pol_avg) * np.ones(len(MJD_adj)),
                 (avg_pol + sigma_pol_avg) * np.ones(len(MJD_adj)),
                 facecolor='orange', alpha=0.3, linestyles='--', color='gray')
plt.ylabel(r'$P_{\rm I}$ $[\%]$', fontsize=fontsize)
plt.xlabel('MJD', fontsize=fontsize)
leg = plt.legend(fontsize=fontsize, loc=4, frameon=False,
                 handlelength=0, handletextpad=0)

for item in leg.legendHandles:
    item.set_visible(False)

plt.minorticks_on()

# ------------------------------------------------------------------------------
# Saving tables and computing averages

avg_U = np.mean(Pol_U)
avg_B = np.mean(Pol_B)
avg_V = np.mean(Pol_V)
avg_R = np.mean(Pol_R)
avg_I = np.mean(Pol_I)

avg_U = Pol_U - avg_U
avg_B = Pol_B - avg_B
avg_V = Pol_V - avg_V
avg_R = Pol_R - avg_R
avg_I = Pol_I - avg_I

avg_err_U = np.mean(err_U)
avg_err_B = np.mean(err_B)
avg_err_V = np.mean(err_V)
avg_err_R = np.mean(err_R)
avg_err_I = np.mean(err_I)

avg_err_U = err_U - avg_err_U
avg_err_B = err_B - avg_err_B
avg_err_V = err_V - avg_err_V
avg_err_R = err_R - avg_err_R
avg_err_I = err_I - avg_err_I

avg_pol = np.concatenate((avg_U, avg_B, avg_V, avg_R, avg_I))
avg_mjd = np.concatenate((MJD_U, MJD_B, MJD_V, MJD_R, MJD_I))
avg_err = np.concatenate((avg_err_U, avg_err_B, avg_err_V,
                          avg_err_R, avg_err_I))

all_filt = np.concatenate((['u'] * np.size(Pol_U), ['b'] * np.size(Pol_B),
                           ['v'] * np.size(Pol_V), ['r'] * np.size(Pol_R),
                           ['i'] * np.size(Pol_I)))

create_txt_file(x=MJD_U, y=Pol_U, z=err_U, w=['u'] * np.size(Pol_U),
                file_name=folder_tab_pol + star + '_pol_U.txt')
create_txt_file(x=MJD_B, y=Pol_B, z=err_B, w=['b'] * np.size(Pol_B),
                file_name=folder_tab_pol + star + '_pol_B.txt')
create_txt_file(x=MJD_V, y=Pol_V, z=err_V, w=['v'] * np.size(Pol_V),
                file_name=folder_tab_pol + star + '_pol_V.txt')
create_txt_file(x=MJD_R, y=Pol_R, z=err_R, w=['r'] * np.size(Pol_R),
                file_name=folder_tab_pol + star + '_pol_R.txt')
create_txt_file(x=MJD_I, y=Pol_I, z=err_I, w=['i'] * np.size(Pol_I),
                file_name=folder_tab_pol + star + '_pol_I.txt')
print(len(avg_mjd), len(avg_pol), len(avg_err), len(all_filt))
create_txt_file(x=avg_mjd, y=avg_pol, z=avg_err, w=all_filt,
                file_name=folder_tab_pol + star + '_pol_all.txt')

# ------------------------------------------------------------------------------
# Plotting fits
if plot_adjust is True:

    # Calculate polynomial
    z = np.polyfit(MJD, Pol, 1)
    f = np.poly1d(z)

    # calculate new x's and y's
    x_new = np.linspace(MJD[0], MJD[-1], 500)
    # x_new = np.linspace(55000., 56800., 5000)
    y_new = f(x_new)
    plt.plot(x_new, y_new, color='black', label="Linear fit")

# ------------------------------------------------------------------------------
ax.legend(ncol=2, shadow=False, title="", fancybox=False, loc=4, fontsize=12)
if plot_adjust is True:
    ax.text(55220, 0.7, 'P$_{adj}$ = %.1e MJD$^1$ + %.1e ' %
            (z[0], z[1]), fontsize=11.)


figure_name = star + '_pol_x_JD' + extension

plt.savefig(fig_folder + figure_name)

print('\nArquivo criado:  %s\n' % figure_name)
print(72 * '-')

# ==============================================================================
# PLOT TWO
show_figure = True
legend_position = (0.25, 0.75)

# ------------------------------------------------------------------------------
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

wu_arr, wb_arr, wv_arr, wr_arr, wi_arr = [], [], [], [], []
for i in range(len(JD)):
    for filt in filtros:

        if Filter[i] == filt and flag[i] is not 'W' and\
           P[i] > 0.38 and P[i] < 0.6:
                if filt == 'u':
                        wave.append(lbd[0])
                        Pol.append(P[i])
                        error.append(SIGMA[i])
                        wu = wu + P[i]  # * 100.
                        wu_arr.append(P[i])
                        nu = nu + 1.
                        eu = eu + (SIGMA[i])**2.
                if filt == 'b':
                        wave.append(lbd[1])
                        Pol.append(P[i])
                        error.append(SIGMA[i])
                        wb = wb + P[i]  # * 100.
                        wb_arr.append(P[i])
                        nb = nb + 1.
                        eb = eb + (SIGMA[i])**2.
                if filt == 'v':
                        wave.append(lbd[2])
                        Pol.append(P[i])
                        error.append(SIGMA[i])
                        wv = wv + P[i]  # * 100.
                        wv_arr.append(P[i])
                        nv = nv + 1.
                        ev = ev + (SIGMA[i])**2.
                if filt == 'r':
                        wave.append(lbd[3])
                        Pol.append(P[i])   # 100. * P[i])
                        error.append(SIGMA[i])
                        wr = wr + P[i]  # * 100.
                        wr_arr.append(P[i])
                        nr = nr + 1.
                        er = er + (SIGMA[i])**2.
                if filt == 'i':
                        wave.append(lbd[4])
                        Pol.append(P[i])
                        error.append(SIGMA[i])
                        wi = wi + P[i]  # * 100.
                        wi_arr.append(P[i])
                        ni = ni + 1.
                        ei = ei + (SIGMA[i])**2.

eu = np.sqrt(eu)
eb = np.sqrt(eb)
ev = np.sqrt(ev)
er = np.sqrt(er)
ei = np.sqrt(ei)

ef = np.array([eu, eb, ev, er, ei])  # * 100

try:
    mean_u = wu / nu
    # mean_u = np.median(wu_arr)
except:
    mean_u = 0.
    print('u sem dado')

try:
    mean_b = wb / nb
    # mean_b = np.median(wb_arr)
except:
    mean_b = 0.
    print('b sem dado')

try:
    mean_v = wv / nv
    # mean_v = np.median(wv_arr)
except:
    mean_v = 0.
    print('v sem dado')

try:
    mean_r = wr / nr
    # mean_r = np.median(wr_arr)
except:
    mean_r = 0.
    print('r sem dado')

try:
    mean_i = wi / ni
    # mean_i = np.median(wi_arr)
except:
    mean_i = 0.
    print('i sem dado')

mean = np.array([mean_u, mean_b, mean_v, mean_r, mean_i])

plt.clf()
plt.close()

fig1 = plt.figure(figsize=(8, 6), dpi=120)
# fig1.subplots_adjust(hspace=0.01, wspace=0.0001, bottom=0.1,
#                      top=0.96, left=0.2, right=1.5)

# define two subplots with different sizes
gs = gridspec.GridSpec(2, 2, width_ratios=[100, 1], height_ratios=[3, 1])
# upper plot
upperplot = plt.subplot(gs[0])

# ------------------------------------------------------------------------------
# Plot OPD data
plt.errorbar(wave, Pol, error, marker='o', linestyle=' ', alpha=0.1,
             label='OPD-Data:\n2009/09/20 - 2014/03/09', color='blue')
# lbd = lbd[1:]
# mean = mean[1:]
# ef = ef[1:]
# plt.errorbar(lbd, mean, ef, marker='o', color='red', alpha=0.5,
#              linestyle=' ', label='OPD-Average Data')

# ------------------------------------------------------------------------------
# Plot hdust model
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
        pol = np.sqrt(((array[j, :, 7]) ** 2)) * 100
        for k in range(len(lbd)):
            tmp = find_nearest(array=lbdarr, value=lbd[k])
            index = tmp[1]
            if angle >= 30. and angle <= 70.:
                peso = 0.1
            else:
                peso = 1.0
            chi2red[j][i] = (chi2red[j][i] + (mean[k] - pol[index])**2 /
                             (ef[k]) ** 2) * peso
            if plot_all_models is True:
                plt.plot(lbdarr, pol, '-', alpha=np.random.random())

chi2min = np.min(chi2red)

# encontrando o indice do minimo chi2
index = np.argwhere(chi2red == chi2min)
index = index[0]

# calculando o chi2 reduzido
m = 2  # numero de graus de liberdade
n = len(lbdarr)
chi2_red_min = chi2min / (n - m - 1)
chi2_red = chi2red / (n - m - 1)

# imprimindo os valores de T1 e T2, correspondentes ao melhor ajuste
bestmodel = file_disk[index[1]]
# nobs = np.arange(10)
best_angle = index[0]
# print(index)
# best_angle = 6

print('Best model: %s' % bestmodel)
n = 'n = ' + str(float(bestmodel[85:88]) + 1.5)
n0 = r'$n_0 = $' + bestmodel[92:100]
Rd = r'$R_D = $' + str(float(bestmodel[108:113])) + r'$R_\odot$'
M = r'$M_\odot = $' + str(float(bestmodel[118:123])) + r'$M_\odot$'
W = r'$W = $' + str(float(bestmodel[125:129]))
beta = r'$\beta = $' + str(float(bestmodel[131:135]))
Rpole = r'$R_{pole} = $' + str(float(bestmodel[138:143])) + r'$R_\odot$'
Lum = r'$L = $' + str(float(bestmodel[145:150])) + r'$ L_\odot$'
label_best_param = [n, n0, Rd, M, W, beta, Rpole, Lum]

# Obs 1 = 90.0, 0.
# Obs 2 = 83.6, 0.
# Obs 3 = 77.2, 0.
# Obs 4 = 70.5, 0.
# Obs 5 = 63.6, 0.
# Obs 6 = 56.3, 0.
# Obs 7 = 48.2, 0.
# Obs 8 = 38.9, 0.
# Obs 9 = 27.3, 0.
# Obs 10 = 00.0, 0.

i_angle = r'$i = $' + str(i_angles[best_angle]) + r' $^o$'

xc = 2200
xc_2 = 4000
delt_y = 0.2
fontsize = 11
plt.annotate(s=n, xy=(xc, 0.94 + delt_y), xycoords='data',
             xytext=(xc, 0.94 + delt_y), fontsize=fontsize)
plt.annotate(s=n0, xy=(xc, 0.88 + delt_y), xycoords='data',
             xytext=(xc, 0.88 + delt_y), fontsize=fontsize)
plt.annotate(s=Rd, xy=(xc, 0.83 + delt_y), xycoords='data',
             xytext=(xc, 0.83 + delt_y), fontsize=fontsize)
plt.annotate(s=i_angle, xy=(xc, 0.78 + delt_y), xycoords='data',
             xytext=(xc, 0.78 + delt_y), fontsize=fontsize)

plt.annotate(s=M, xy=(xc_2, 0.94 + delt_y), xycoords='data',
             xytext=(xc_2, 0.94 + delt_y), fontsize=fontsize)
plt.annotate(s=W, xy=(xc_2, 0.88 + delt_y), xycoords='data',
             xytext=(xc_2, 0.88 + delt_y), fontsize=fontsize)
plt.annotate(s=beta, xy=(xc_2, 0.83 + delt_y), xycoords='data',
             xytext=(xc_2, 0.83 + delt_y), fontsize=fontsize)
plt.annotate(s=Rpole, xy=(xc_2, 0.78 + delt_y), xycoords='data',
             xytext=(xc_2, 0.78 + delt_y), fontsize=fontsize)
plt.annotate(s=Lum, xy=(xc_2, 0.73 + delt_y), xycoords='data',
             xytext=(xc_2, 0.73 + delt_y), fontsize=fontsize)

# array = plot_hdust_model(fullsedfile=file_disk[index])
array = plot_hdust_model(fullsedfile=bestmodel)
lbdarr = array[:, :, 2] * 1.e4  # Angstrom
pol = np.sqrt(array[best_angle, :, 7] ** 2) * 100.

cmap = 'gist_rainbow'
cor = phc.gradColor(np.arange(nobs), cmapn=cmap)
for j in range(nobs):
    if i_angles[j] <= 1.4 * ref_angle and i_angles[j] >= 0.6 * ref_angle:
        lbdarr = array[j, :, 2] * 1.e4  # Angstrom
        pol = np.sqrt(array[j, :, 7] ** 2) * 100.
        plt.plot(lbdarr, pol, linestyle='--', color=cor[j],
                 alpha=0.7, label='i = {} $^o$'.format(i_angles[j]),
                 linewidth=2)

lbdarr = array[best_angle, :, 2] * 1.e4  # Angstrom
pol = np.sqrt(array[best_angle, :, 7] ** 2) * 100.
plt.plot(lbdarr, pol, '-', color=cor[best_angle], linewidth=3,
         label='Best model')
plt.errorbar(lbd, mean, ef, marker='o', color='red', alpha=1.,
             linestyle=' ', label='OPD-Average Data', markersize=10,
             markeredgecolor='black', markeredgewidth=2)

plt.legend(loc='best', fontsize=10, numpoints=1, ncol=2, framealpha=0.)

plt.xlim(2000, 12000)
plt.ylim(0, 1.2)

# plt.xlabel(r'$\lambda$ [$\AA$]')
plt.ylabel('P $[\%]$')
ax.legend(ncol=1, shadow=False, title="", fancybox=True,
          prop={'size': 7.8}, numpoints=1,
          bbox_to_anchor=(1.01, 1), loc=2)

plt.minorticks_on()
plt.rc('font', family='serif', size=13)

plt.tight_layout()
plt.setp(upperplot.get_xticklabels(), visible=False)

# plt.autoscale(enable=True)

# lower plot - residuals
res = np.zeros(len(mean))
err_res = []
for i in range(len(mean)):
    tmp = find_nearest(array=lbdarr, value=lbd[i])
    index = tmp[1]
    lbd_tmp = lbdarr[index]
    tmp_2 = find_nearest(array=lbd, value=lbd_tmp)
    lbd_tmp = tmp_2[0]

    # lbd = np.array([3656., 4353., 5477., 6349., 8797.])
    if lbd_tmp == 3656.:
        filt = 'u'
    elif lbd_tmp == 4353.:
        filt = 'b'
    elif lbd_tmp == 5477.:
        filt = 'v'
    elif lbd_tmp == 6349.:
        filt = 'r'
    elif lbd_tmp == 8797.:
        filt = 'i'

    pol_conv = phd.doFilterConv(lbdarr, pol,
                                filt=filt.upper(), zeropt=True)
    print(lbd_tmp, filt, pol_conv, mean[i])
    res[i] = (mean[i] - pol_conv) / ef[i]
    err_tmp = np.sqrt(ef[i])
    err_res.append(err_tmp)

lowerplot = plt.subplot(gs[2])
plt.errorbar(lbd, res, err_res, marker=r'o', markersize=8,
             linestyle='', markeredgecolor='black', markeredgewidth=2,
             markerfacecolor='white')
plt.xlim(2000, 12000)
plt.hlines(y=0, xmin=2000, xmax=12000, linestyle='--')
plt.minorticks_on()
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.ylabel('Residuals')

# plt.colorbar()
plt.tight_layout()
figure_name = star + '_specpol' + extension
plt.savefig(fig_folder + figure_name, dpi=600)
# plt.show()
plt.close()
print('\nArquivo criado:  %s\n' % figure_name)


# ==============================================================================
# HPOL Plot (MAST  files)
if plot_mast is True:
    print('\nPlotting second plot\n')
    print(72 * '-')
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

    path = len(files) * [tab_folder + 'hpol/']
    files = [x + y for x, y in zip(path, files)]

    fmt = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # 2 para arquivo .s do
    shift = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    # DO FIGURE
    if individual is True:
        fig, axes = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True)
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, axisbg='white')
        ax.set_xlim(3000., 9000.)
        ax.set_ylabel('P $[\%]$')
        ax.set_xlabel('$\lambda$ [$\AA$]')

        cmap = 'gist_rainbow'
        colors = phc.gradColor(np.arange(len(files)), cmapn=cmap)

        splot(files=files, fig_folder=fig_folder, fmt=fmt, shift=shift,
              wave_t=wave, Pol_t=Pol, error_t=error, lbd_t=lbd, mean_t=mean,
              ef_t=ef, scalarMap=colors)

    else:
        fig, axes = plt.subplots(nrows=4, ncols=2, sharex=True, sharey=True)
        colors = phc.gradColor(np.arange(len(files)), cmapn='rainbow')
        for i in range(len(axes)):
            for j in range(len(axes[0][:])):
                for k in range(np.size(files)):
                    file = [files[k]]
                    ax = fig.add_subplot(axes[i, j], axisbg='white')

                    cmap = 'gist_rainbow'
                    colors = phc.gradColor(np.arange(len(files)), cmapn=cmap)

                    bin_lbda, bin_pol, error = plot_fits(file, fmt, shift,
                                                         colors)
                    plt.subplots_adjust(right=0.85)
                    plt.plot(bin_lbda, bin_pol, color=colors[k])
                    plt.errorbar(bin_lbda, bin_pol, yerr=error,
                                 ecolor=colors, fmt='kp--')
                    plt.ylim(ymin, ymax)
                    plt.xlim(xmin, xmax)
                    plt.xlabel(r'$\lambda$ [$\AA$]')
                    plt.ylabel('P $[\%]$')
                    plt.minorticks_on()
                    plt.autoscale(enable=autoscale, axis=u'both', tight=False)
                    plt.rc('font', family='serif', size=10)

# ------------------------------------------------------------------------------

# for i in range(len(file_disk)):
#     array = plot_hdust_model(fullsedfile=file_disk[i])
#     lbdarr = array[:, :, 2] * 1.e4
#     if np.shape(lbdarr)[1] != 452:
#         pol = np.sqrt(array[:, :, 7] ** 2) * 100.
#         print(np.shape(lbdarr)[1], np.shape(pol)[1], file_disk[i][76:])

