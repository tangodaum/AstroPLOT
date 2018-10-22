# PLOT VOSA AND VOSPEC XML FILES

import numpy as np
import matplotlib.pyplot as plt
import astropy
import atpy
from matplotlib.widgets import Cursor
import pyhdust.spectools as spt
import pyhdust
from PyAstronomy import pyasl
import pyfits
import os
import pyhdust.phc as phc
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib.ticker import MaxNLocator
import pyhdust.beatlas as bat
import pyhdust as hdt
from scipy.interpolate import griddata

# ==============================================================================
# Basic data
star_name = 'alf_ara'  # 'alf_ara'  'acol'
show_figure = True
autoscale = False
plot_kurucz = False
loop = False  # Plot more than one star
computer = 'tangodown'
fits = True
plot_hdust_model = True
plot_vospec = True
plot_vosa = True
markersize = 5
cmap = 'jet'  # 'gist_rainbow', 'jet'
plot_all_models = False
show_integral = False

# Define the data folder
folder_data = '/Users/' + computer + '/Dropbox/2_Artigos/tex_aara/emcee/iue/'
folder_data_bump = '/Users/' + computer +\
    '/Dropbox/2_Artigos/1_Projeto/emcee/iue/'
commom_folder = '/Users/tangodown/Dropbox/2_Artigos/tex_aara/Bemcee/'
# commom_folder = '/home/tangodown/Dropbox/2_artigos/tex_aara/Bemcee/'
folder_models = commom_folder + 'models/'

# ==============================================================================
# Stellar Parameters and constants
# Read the grid models, with the interval of parameters.
xdrPL = folder_models + 'aara_final.xdr'  # 'PL.xdr'
# xdrPL = folder_models + 'aara_acs.xdr'  # 'PL.xdr'
# xdrPL = folder_models + 'disk_flx.xdr'  # 'PL.xdr'

listpar, lbdarr, minfo, models = bat.readBAsed(xdrPL, quiet=False)
lbdarr = lbdarr * 1e4  # AA
# F(lbd)] = 10^-4 erg/s/cm2/Ang

for i in range(np.shape(minfo)[0]):
    for j in range(np.shape(minfo)[1]):
        if minfo[i][j] < 0:
            minfo[i][j] = 0.

for i in range(np.shape(models)[0]):
    for j in range(np.shape(models)[1]):
        if models[i][j] < 0. or models[i][j] == 0.:
            models[i][j] = (models[i][j + 1] + models[i][j - 1]) / 2.

# n0 to logn0
listpar[4] = np.log10(listpar[4])
listpar[3].sort()
# minfo[:, 4] = np.log10(minfo[:, 4])
for i in range(len(minfo)):
    minfo[i][4] = np.log10(minfo[i][4])

# delete columns of fixed par
cols2keep = [0, 1, 3, 4, 5, 7, 8]
cols2delete = [2, 6]
listpar = [listpar[i] for i in cols2keep]
minfo = np.delete(minfo, cols2delete, axis=1)
listpar[3].sort()

# print(params_tmp)
dims = ['M', 'ob', 'Hfrac', 'sig0', 'Rd', 'mr', 'cosi']
dims = dict(zip(dims, range(len(dims))))
isig = dims["sig0"]

# ==============================================================================
# Stellar Parameters and constants
G = astropy.constants.G.cgs
M_sun = astropy.constants.M_sun.cgs
R_sun = astropy.constants.R_sun.cgs
pc = astropy.constants.pc.cgs
L_sun = astropy.constants.L_sun.cgs  # erg/s

# Stellar Parameters from BeFaVOr
# logg = 4.15  # 3.95
# R_star = 3.34 * R_sun  # 5.85 * R_sun
# L_star = 1500 * L_sun
L_star = 1940.54 * L_sun

# Put here the best parameters obtained with the emcee
params = [6.85, 1.30, 0.44, 12.87, 15.96, 1.35, 0.68594536278939489]
distance = 86.12  # * pc
distance_nodisk = distance * pc

# Define here the ebmv inferred by some method
ebmv = 0.03
ebmv_ps = 0.01
ebmv_ms = 0.01

angle_index = 7
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

# file_star = 'fullsed_mod01_noCS_Be_M06.57_W0.76_b0.18_rp04.57_L02180_Ell.sed2'
file_star = 'fullsed_mod01_noCS_Be_M06.85_W0.96_b0.18_rp03.81_L01941_Ell.sed2'

common_folder = '/Users/' + computer +\
    '/Dropbox/2_Artigos/tex_aara/photometry/'
# model_nodisk = common_folder + 'hdust_models/' + 'fullsed_mod01a-phot.sed2'
# model_disk = common_folder + 'hdust_models/' + 'fullsed_mod01a.sed2'

model_nodisk = common_folder + 'hdust_models/' + file_star
# model_disk = common_folder + 'hdust_models/' + file_disk
model_disk = common_folder + 'hdust_models/'
# model_disk = '/Users/tangodown/Dropbox/2_Artigos/tex_aara/tobe_run/' +\
   # 'aara_acs/fullsed/'
folder_tab = common_folder + 'tables/'



# ==============================================================================
def griddataBA(minfo, models, params, listpar, dims):
    '''
    Moser's routine to interpolate BeAtlas models
    obs: last argument ('listpar') had to be included here
    '''
    idx = np.arange(len(minfo))
    lim_vals = len(params) * [[], ]
    for i in range(len(params)):
        # print(i, listpar[i], params[i], minfo[:, i])
        lim_vals[i] = [
            phc.find_nearest(listpar[i], params[i], bigger=False),
            phc.find_nearest(listpar[i], params[i], bigger=True)]
        tmp = np.where((minfo[:, i] == lim_vals[i][0]) |
                       (minfo[:, i] == lim_vals[i][1]))
        idx = np.intersect1d(idx, tmp[0])

    out_interp = griddata(minfo[idx], models[idx], params)[0]

    if (np.sum(out_interp) == 0 or np.sum(np.isnan(out_interp)) > 0):

        mdist = np.zeros(np.shape(minfo))
        ichk = range(len(params))
        for i in ichk:
            mdist[:, i] = np.abs(minfo[:, i] - params[i])/(np.max(listpar[i]) - 
                np.min(listpar[i]))
        idx = np.where(np.sum(mdist, axis=1) == np.min(np.sum(mdist, axis=1)))
        if len(idx[0]) != 1:
            out_interp = griddata(minfo[idx], models[idx], params)[0]
        else:
            out_interp = models[idx][0]

    # if (np.sum(out_interp) == 0 or np.sum(np.isnan(out_interp)) > 0) or\
    #    bool(np.isnan(np.sum(out_interp))) is True:
    #     print("# Houve um problema na grade e eu nao consegui arrumar...")

    return out_interp


# ==============================================================================
def mass_density(n_0):
    NA = 6.022140858e23  # mol−1
    M_H = 1.00798  # g/mol
    p_m = (M_H * n_0) / NA
    return p_m


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
def find_nearest(array, value):
    '''
    Find the nearest value inside an array
    '''

    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


# ==============================================================================
# Use a vo format xml
os.chdir('/Users/' + computer + '/Dropbox/2_Artigos/tex_aara/photometry/')
if plot_vosa is True:
    # VOSA VO file
    fits_file1 = common_folder + 'vo_data/' + 'vosa_aara.xml'

# VOSPEC VO File
# 'vospec_aara.xml'
fits_file2 = common_folder + 'vo_data/' + 'vospec_european_hst_ssap.xml'
fits_file3 = common_folder + 'vo_data/' + 'avg_iue_aara.xml'
# ==============================================================================
if plot_vosa is True:
    t1 = atpy.Table(fits_file1)  # read VOSA

t2 = atpy.Table(fits_file2, tid=1)  # read VOSpec
t3 = atpy.Table(fits_file3, tid=1)

# ==============================================================================
# Para acessar os nomes das colunas:
if plot_vosa is True:
    t1.columns

t2.columns

# ==============================================================================
# Define os arrays:
if plot_vosa is True:
    # VOSA
    flux1 = t1['Flux'][:]  # erg/cm2/s/A
    wave1 = t1['Wavelength'][:]  # Angstrom
    Error1 = t1['Error'][:]  # erg/cm2/s/A
    catalogue = t1['Service'][:]
    flux1_der = pyasl.unred(wave1, flux1, ebv=ebmv, R_V=3.1)
    flux1_der_psig = pyasl.unred(wave1, flux1, ebv=ebmv_ps, R_V=3.1)
    flux1_der_msig = pyasl.unred(wave1, flux1, ebv=ebmv_ms, R_V=3.1)


# VOSpec
flux2 = t2['Flux0'][:]  # erg/cm2/s/A
wave2 = t2['SpectralAxis0'][:]  # Angstrom
flux3 = t3['Flux0'][:]  # erg/cm2/s/A
wave3 = t3['SpectralAxis0'][:]  # Angstrom


# To calculate the Reddening and extinction, we have used this link
# Fitzpatrick law
flux2_der = pyasl.unred(wave2, flux2, ebv=ebmv, R_V=3.1)
flux2_der_psig = pyasl.unred(wave2, flux2, ebv=ebmv_ps, R_V=3.1)
flux2_der_msig = pyasl.unred(wave2, flux2, ebv=ebmv_ms, R_V=3.1)
flux3_der = pyasl.unred(wave3, flux3, ebv=ebmv, R_V=3.1)
flux3_der_psig = pyasl.unred(wave3, flux3, ebv=ebmv_ps, R_V=3.1)
flux3_der_msig = pyasl.unred(wave3, flux3, ebv=ebmv_ms, R_V=3.1)

# Units
erg = astropy.units.erg
second = astropy.units.s
cm2 = astropy.units.cm**2

# ==============================================================================
# Plotting Observational Data
fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, axisbg='white')
ax = plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=3)  # projection='3d')

if plot_vosa is True:
    missions = np.unique(catalogue)
    symbols = np.unique(catalogue)
    colors = phc.gradColor(np.arange(len(missions)), cmapn=cmap)
    colorsdict = dict(zip(missions, colors))
    list_mission = []
    Error1 = Error1 * wave1  # / flux1
    # Error1 = 0.434 * Error1 / (wave1 * flux1)
    for mission in catalogue:
        idx = np.where(catalogue == mission)
        colors = colorsdict[mission]

        if (mission in list_mission) is False:
            list_mission.append(mission)
            plt.errorbar(wave1[idx], flux1[idx] * wave1[idx],
                         yerr=Error1[idx],
                         fmt='.', label=mission.decode('UTF-8'), color=colors,
                         markersize=2 * markersize, elinewidth=1,
                         markeredgewidth=1)
        else:
            plt.errorbar(wave1[idx], flux1[idx] * wave1[idx],
                         yerr=Error1[idx],
                         fmt='.', color=colors, markersize=2 * markersize,
                         elinewidth=1, markeredgewidth=1)

            for i in range(len(wave1)):
                if wave1[i] > 1e9:
                    plt.errorbar(wave1[i], flux1[i] * wave1[i],
                                 yerr=np.full(1, 0.9 * flux1[i] * wave1[i]),
                                 uplims=True, color=colors,
                                 markeredgecolor="r", markerfacecolor="r",
                                 linestyle="None", alpha=0.5)

# Plotting VOSpec
if plot_vospec is True:
    # ax.plot(wave2, flux2 * wave2, '-', label='VOSpec',
    #         color='red', markersize=markersize)
    ax.plot(wave2, flux2_der * wave2, '-', label='HST',
            color='gray', markersize=markersize, alpha=0.2)
    ax.plot(wave3, flux3_der * wave3, '-', label='IUE',
            color='gray', markersize=markersize, alpha=0.5)
    # ax.plot(wave2, flux2_der_psig * wave2, 'o', label='E(B-V) = 0.18',
    #        color='green', markersize=2)
    # ax.plot(wave2, flux2_der_msig * wave2, 'o', label='E(B-V) = 0.12',
    #        color='orange', markersize=2)


# ==============================================================================
# Convert to boolean
def to_bool(value):
    valid = {'true': True, 't': True, '1': True,
             'false': False, 'f': False, '0': False}

    if isinstance(value, bool):
        return value

    if not isinstance(value, basestring):
        raise ValueError('invalid literal for boolean. Not a string.')

    lower_value = value.lower()
    if lower_value in valid:
        return valid[lower_value]
    else:
        raise ValueError('invalid literal for boolean: "%s"' % value)


# ==============================================================================
# Creating the list of stars
if loop is True:
    folder_table = '/home/' + computer +\
        '/MEGA/1 Tese/1 Projeto/48lib/emcee/tables/'
    file_name = 'emcee_bestars.txt'
    typ = (0, 1, 2, 3)
    file_data = folder_table + file_name
    a = np.loadtxt(file_data, usecols=typ, unpack=True, delimiter='\t',
                   dtype={'names': ('star', 'plx', 'sig_plx', 'bump'),
                          'formats': ('S9', 'f2', 'f2', 'S5')})
    stars, plx, sig_plx, bump0 = a[0], a[1], a[2], a[3]
    bump = []
    for i in range(len(bump0)):
        bump.append(to_bool(bump0[i]))
    bump = np.array(bump)

# ==============================================================================
# Opening FIT files from INES
    for i in range(len(stars)):
        print(stars[i])
        if fits is True:
            table = 'list.txt'
            os.chdir(folder_data + 'ines-' + stars[i] + '/')
        if os.path.isfile(table) is False or os.path.isfile(table) is True:
            os.system('ls *.FITS | xargs -n1 basename >list.txt')
            table = folder_data + 'ines-' + stars[i] + '/' + table
            iue_list = np.loadtxt(table, comments='#',
                                  dtype=[('file_name', '|S41')])
            file_name = iue_list['file_name']
        fluxes = []
        waves = []
        errors = []
        if np.size(iue_list) == 1:
            file_iue = str(folder_data) + 'ines-' + stars[i] + '/' +\
                str(file_name)
            hdulist = pyfits.open(file_iue)
            tbdata = hdulist[1].data
            wave = tbdata.field('WAVELENGTH')  # Angs
            flux = tbdata.field('FLUX')  # erg/cm2/s/A
            sigma = tbdata.field('SIGMA')
            plt.errorbar(wave, flux * wave, yerr=sigma, fmt=".k",
                         markersize=markersize)
            fluxes.append(flux)
            waves.append(wave)
            errors.append(sigma)
        else:
            for j in range(len(file_name)):
                file_iue = folder_data + 'ines-' + stars[i] + '/' +\
                    file_name[j]
                hdulist = pyfits.open(file_iue)
                tbdata = hdulist[1].data
                wave = tbdata.field('WAVELENGTH')  # Angs
                flux = tbdata.field('FLUX')  # erg/cm2/s/A
                sigma = tbdata.field('SIGMA')
                quality = tbdata.field('QUALITY')
                plt.errorbar(wave, flux * wave, yerr=sigma, fmt=".k")
                fluxes.append(flux)
                waves.append(wave)
                errors.append(sigma)
        # bump = False #********
        if bump[i] is True:
            table_bump = 'list.txt'
            os.chdir(folder_data_bump + 'ines-' + stars[i] + '/bump/')
            if os.path.isfile(table_bump) is False or \
               os.path.isfile(table_bump) is True:
                os.system('ls *.FITS | xargs -n1 basename >list.txt')
                table_bump = folder_data_bump + 'ines-' + stars[i] +\
                    '/bump/' + table_bump
                iue_bump_list = np.loadtxt(table_bump, comments='#',
                                           dtype=[('file_name_iue', '|S41')])
                file_name_bump = iue_bump_list['file_name_iue']
            fluxes_bump = []
            waves_bump = []
            errors_bump = []
            if np.size(iue_bump_list) == 1:
                file_iue_bump = str(folder_data_bump) + 'ines-' + stars[i] +\
                    '/bump/' + str(file_name_bump)
                hdulist = pyfits.open(file_iue_bump)
                tbdata = hdulist[1].data
                wave_bump = tbdata.field('WAVELENGTH')  # Angs
                flux_bump = tbdata.field('FLUX')  # erg/cm2/s/A
                sigma_bump = tbdata.field('SIGMA')
                quality_bump = tbdata.field('QUALITY')
                plt.errorbar(wave_bump, flux_bump * wave_bump,
                             yerr=sigma_bump, fmt=".k")
                fluxes_bump.append(flux_bump)
                waves_bump.append(wave_bump)
                errors_bump.append(sigma_bump)
            else:
                for j in range(len(file_name_bump)):
                    file_iue_bump = folder_data_bump + 'ines-' + stars[i] +\
                        '/bump/' + file_name_bump[j]
                    hdulist = pyfits.open(file_iue_bump)
                    tbdata = hdulist[1].data
                    wave_bump = tbdata.field('WAVELENGTH')  # Angs
                    flux_bump = tbdata.field('FLUX')  # erg/cm2/s/A
                    sigma_bump = tbdata.field('SIGMA')
                    quality = tbdata.field('QUALITY')
                    plt.errorbar(wave_bump, flux_bump * wave_bump,
                                 yerr=sigma_bump, fmt=".k")
                    fluxes_bump.append(flux_bump)
                    waves_bump.append(wave_bump)
                    errors_bump.append(sigma_bump)

# ==============================================================================
# Kurucz Models - Plotting some Kurucz's models
if plot_kurucz is True:
    Ts = np.arange(16000., 24000., 4000.)
    factor = len(Ts)

    wav = []
    flx = []
    inf = []
    for i in range(len(Ts)):
        # [flx] = ergs/cm**2/s/Hz/ster   [wav] = nm
        wav, flx, inf = spt.kuruczflux(Ts[i], logg, [100, 1E5])
        wave = 10. * wav[:]  # Angstrom
        flux = ((3.336 * 1E-19 * (wave**2) *
                ((4. * np.pi)**(-1.)))**(-1)) * flx
        flux = flux * ((R_star / distance)**2)
        ax.plot(wave, flux * wave, '-',
                label='T={0} g={1}'.format(inf[0], inf[1]))

# ==============================================================================
if plot_hdust_model is True:
    # Plotting a hdust model
    global label_best_param
    norm = L_star / (4. * np.pi * (distance_nodisk)**2)

    model_hdt = pyhdust.readfullsed2(model_nodisk)
    wave_model = 1E4 * model_hdt[angle_index, :, 2]  # Angstrom
    wave_model_nodisk = np.copy(wave_model)
    flux_model = 1E-4 * model_hdt[angle_index, :, 3] * norm  # erg/cm2/s
    flux_model_nodisk = np.copy(flux_model)
    ax.plot(wave_model_nodisk, flux_model_nodisk * wave_model_nodisk,
            '-', label='$\star$ model', color='blue')

    # model_hdt = pyhdust.readfullsed2(model_disk)
    # wave_model = 1E4 * model_hdt[angle_index, :, 2]  # Angstrom
    # wave_model_disk = np.copy(wave_model)
    # flux_model = 1E-4 * model_hdt[angle_index, :, 3] * norm  # erg/cm2/s
    # flux_model_disk = np.copy(flux_model)
    # ax.plot(wave_model_disk, flux_model_disk * wave_model_disk, '--',
    #         label='$\star$+disk model', color='red')

    os.system('rm ' + folder_tab + 'files_disks.txt')
    create_list_files(list_name='files_disks', folder=model_disk,
                      folder_table=folder_tab)

    file_disk = read_list_files_all(table_name='files_disks.txt',
                                    folder_table=folder_tab)

    chi2_red = []
    file_disk = file_disk[1:]
    for i in range(len(file_disk)):
        # model_disk = common_folder + 'hdust_models/' + str(file_disk[i])
        model_disk = str(file_disk[i])
        model_hdt = pyhdust.readfullsed2(model_disk)
        wave_model = 1E4 * model_hdt[angle_index, :, 2]  # Angstrom
        # print(file_disk[i])
        # print(len(wave_model))
        wave_model_disk = np.copy(wave_model)
        flux_model = 1E-4 * model_hdt[angle_index, :, 3] * norm  # erg/cm2/s
        flux_model_disk = np.copy(flux_model)
        if plot_all_models is True:
            ax.plot(wave_model_disk, flux_model_disk * wave_model_disk, '--',
                    label='$\star$+disk model', color='red',
                    alpha=np.random.random())
        chi2red = 0.
        for j in range(len(wave1)):
            tmp = find_nearest(array=wave_model_disk, value=wave1[j])
            index = tmp[1]
            if wave1[j] < 1.e9 and wave1[j] > 1.e6:
                chi2red = chi2red + (flux1_der[j] -
                                     flux_model_disk[index])**2 /\
                    (flux_model_disk[index]) ** 2
            else:
                chi2red = np.inf

        chi2_red.append(chi2red)

    # Plotting best model
    min_chi2red = np.min(chi2_red)
    tmp = find_nearest(array=chi2_red, value=min_chi2red)
    index = tmp[1]

    # model_disk = common_folder + 'hdust_models/' + str(file_disk[index])
    model_disk = str(file_disk[index])

    model_hdt = pyhdust.readfullsed2(model_disk)
    best_wave_model_disk = 1E4 * model_hdt[angle_index, :, 2]  # Angstrom
    # Flux in erg/cm2/s
    best_flux_model_disk = 1E-4 * model_hdt[angle_index, :, 3] * norm

    # norm = L_star / (4. * np.pi * (distance)**2)
    norma = (10. / distance)**2

    # select lbdarr to coincide with lbd
    models = np.log10(models)
    best_flux_model_disk = griddataBA(minfo, models, params, listpar, dims)
    best_flux_model_disk += np.log10(norma)
    best_flux_model_disk = 10**best_flux_model_disk

    # best_flux_model_disk = best_flux_model_disk * norm
    ax.plot(lbdarr, best_flux_model_disk * lbdarr * 1e-4,
            '--', label='best model', color='green')
    # ax.plot(best_wave_model_disk, best_flux_model_disk * best_wave_model_disk,
    #         '--', label='best model', color='green')

    # print('Best model: %s' % file_disk[index])
    # n = 'n = ' + str(float(model_disk[85:88]) + 1.5)
    # n0 = r'$n_0 = $' + model_disk[92:100]
    # Rd = r'$R_D = $' + str(float(model_disk[108:113])) + r'$R_\odot$'
    # M = r'$M_\odot = $' + str(float(model_disk[118:123])) + r'$M_\odot$'
    # W = r'$W = $' + str(float(model_disk[125:129]))
    # beta = r'$\beta = $' + str(float(model_disk[131:135]))
    # Rpole = r'$R_{pole} = $' + str(float(model_disk[138:143])) + r'$R_\odot$'
    # Lum = r'$L = $' + str(float(model_disk[145:150])) + r'$ L_\odot$'
    # label_best_param = [n, n0, Rd, M, W, beta, Rpole, Lum]

    # remove the errorbars
    # use them in the legend
    # legend1 = ax.legend(label_best_param, ncol=1,
    #                     shadow=False, fancybox=True,
    #                     prop={'size': 7.9}, numpoints=1,
    #                     scatterpoints=0,
    #                     markerscale=0., markerfirst=False,
    #                     bbox_to_anchor=(1.25, -0.01), loc=4,
    #                     fontsize='x-large')
    # plt.gca().add_artist(legend1)

# ==============================================================================
    # Fill between regions
    flux_model_nodisk = np.interp(lbdarr, wave_model_nodisk,
                                  flux_model_nodisk)
    wave_model_nodisk = np.copy(lbdarr)
    # wave_model_disk = np.copy(best_wave_model_disk)
    wave_model_disk = np.copy(lbdarr)

    ax.fill_between(lbdarr,
                    best_flux_model_disk * lbdarr * 1e-4,
                    flux_model_nodisk * lbdarr, facecolor='orange',
                    alpha=0.5, linestyles='--', color='orange')
    # ax.fill_between(lbdarr,
    #                 best_flux_model_disk * lbdarr * 1e-4,
    #                 flux_model_nodisk * lbdarr, facecolor='orange',
    #                 alpha=0.5, linestyles='--', color='orange')


    # print(len(wave_model_disk), len(best_flux_model_disk),
    #       len(best_wave_model_disk), len(flux_model_nodisk),
    #       len(wave_model_nodisk))
    # int_flux_disk = np.trapz(y=best_flux_model_disk * best_wave_model_disk,
    #                          x=wave_model_disk, dx=10)

    # int_flux_star = np.trapz(y=flux_model_nodisk * wave_model_nodisk,
    #                          x=wave_model_nodisk, dx=10)

    # flux_ratio = int_flux_disk / int_flux_star
    # flux_ratio = flux_ratio.value

    # Annotate
    # if show_integral is True:
    #     ax.annotate(s=r'$\int_{-\infty}^{+\infty}$' +
    #                 r'$\frac{F_{(\star + csd)}}{F_\star} d\lambda = $' + '{:.2f}'.
    #                 format(flux_ratio), xy=(1.3e7, 1e-14),
    #                 xycoords='data', xytext=(2.e5, 1.3e-7),
    #                 arrowprops=dict(arrowstyle="<-",
    #                                 connectionstyle="angle, angleA=180," +
    #                                                 "angleB=270, rad=5"),)

# ==============================================================================
    # Zooming some region
    zoom = True
    zoom_ay = 2.5
    axylim = [0.005, 0.05]
    ayylim = [0.5e-11, 1e-8]
    wave_lim_min = 1.2e3
    wave_lim_max = 5.e4
    wave_lim_min2 = 1.e-8
    wave_lim_max2 = 1.e-5

    if zoom is True:

        axins = zoomed_inset_axes(ax, zoom_ay, loc=3)
        axins.axis([wave_lim_min, wave_lim_max, wave_lim_min2, wave_lim_max2])
        axins.set_xlim(wave_lim_min, wave_lim_max)
        axins.xaxis.tick_top()
        axins.set_ylim(wave_lim_min2, wave_lim_max2)
        mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="0.5", alpha=0.3)
        axins.xaxis.set_major_locator(MaxNLocator(nbins=1, prune='lower'))

        # UV region
        axins.annotate(s='', xy=(1.2e3, 5.e-8),
                       xycoords='data', xytext=(4.e3, 5.e-8),
                       annotation_clip=False,
                       arrowprops=dict(arrowstyle="<->",
                                       connectionstyle="angle"))
        axins.text(x=2.e3, y=2.e-8, s='UV', fontsize=10)
        axins.vlines(x=4.e3, ymin=1.e-8, ymax=1.e-5,
                     linestyles='--', alpha=0.3)

        # Visual region
        axins.annotate(s='', xy=(4.e3, 5.e-8),
                       xycoords='data', xytext=(6.5e3, 5.e-8),
                       annotation_clip=False,
                       arrowprops=dict(arrowstyle="<->",
                                       connectionstyle="angle"))
        axins.text(x=4.2e3, y=2.e-8, s='VIS', fontsize=10)
        axins.vlines(x=6.5e3, ymin=1.e-8, ymax=1.e-5,
                     linestyles='--', alpha=0.3)

        # NIR region
        axins.annotate(s='', xy=(6.5e3, 5.e-8),
                       xycoords='data', xytext=(5.1e4, 5.e-8),
                       annotation_clip=False,
                       arrowprops=dict(arrowstyle="<->",
                                       connectionstyle="angle"))
        axins.text(x=1.3e4, y=2.e-8, s='NIR', fontsize=10)
        axins.vlines(x=5.e4, ymin=1.e-8, ymax=1.e-5,
                     linestyles='--', alpha=0.3)

        xc = 7.5e3
        xc_2 = 2e4
        fontsize = 7.5
        # axins.annotate(s=n, xy=(xc, 6e-6), xycoords='data',
        #                xytext=(xc, 6e-6), fontsize=fontsize)
        # axins.annotate(s=n0, xy=(xc, 3.5e-6), xycoords='data',
        #                xytext=(xc, 3.5e-6), fontsize=fontsize)
        # axins.annotate(s=Rd, xy=(xc, 2.0e-6), xycoords='data',
        #                xytext=(xc, 2.0e-6), fontsize=fontsize)
        # axins.annotate(s=M, xy=(xc, 1.1e-6), xycoords='data',
        #                xytext=(xc, 1.1e-6), fontsize=fontsize)

        # axins.annotate(s=W, xy=(xc_2, 6e-6), xycoords='data',
        #                xytext=(xc_2, 6e-6), fontsize=fontsize)
        # axins.annotate(s=beta, xy=(xc_2, 3.5e-6), xycoords='data',
        #                xytext=(xc_2, 3.5e-6), fontsize=fontsize)
        # axins.annotate(s=Rpole, xy=(xc_2, 2e-6), xycoords='data',
        #                xytext=(xc_2, 2e-6), fontsize=fontsize)
        # axins.annotate(s=Lum, xy=(xc_2, 1.1e-6), xycoords='data',
        #                xytext=(xc_2, 1.1e-6), fontsize=fontsize)

        # In the normal plot (NIR region)
        ymax_vline = 1.e-4
        ymin_vline = 2.e-5
        alpha = 0.2
        arrow_y_position = 0.4e-5
        text_y_position = 0.8e-5

        # MIR region  7.5 to 25.
        ax.annotate(s='', xy=(7.5e4, arrow_y_position),
                    xycoords='data', xytext=(25.e4, arrow_y_position),
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle="<->",
                                    connectionstyle="angle"))
        ax.text(x=0.9e5, y=text_y_position, s='MIR', fontsize=10)
        ax.vlines(x=7.5e4, ymin=ymin_vline, ymax=ymax_vline,
                  linestyles='--', alpha=alpha)
        ax.vlines(x=25.e4, ymin=ymin_vline, ymax=ymax_vline,
                  linestyles='--', alpha=alpha)

        # FIR region  28 to 450.
        ax.annotate(s='', xy=(28.e4, arrow_y_position),
                    xycoords='data', xytext=(450.e4, arrow_y_position),
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle="<->",
                                    connectionstyle="angle"))
        ax.text(x=8.e5, y=text_y_position, s='FIR', fontsize=10)
        ax.vlines(x=450.e4, ymin=ymin_vline, ymax=ymax_vline,
                  linestyles='--', alpha=alpha)

        # MM region  1e7 to 1e10.
        # RADIO region  1e7 to 1e6.
        ax.annotate(s='', xy=(450.e4, arrow_y_position),
                    xycoords='data', xytext=(1.e10, arrow_y_position),
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle="<->",
                                    connectionstyle="angle"))
        ax.text(x=5e7, y=text_y_position, s='MICROWAVE', fontsize=10)

        #
        axins.set_yscale("log")
        axins.set_xscale("log")
        axins.minorticks_on()
        axins.set_yticks([])
        axins.set_xticks([])

        # Plotting
        # axins.plot(best_wave_model_disk, best_flux_model_disk *
        #            best_wave_model_disk, '--', label='best model',
        #            color='green')
        axins.plot(lbdarr, best_flux_model_disk * lbdarr * 1e-4,
                   '--', label='best model', color='green')

        axins.plot(wave_model_nodisk, flux_model_nodisk * wave_model_nodisk,
                   '-', label='$\star$ model', color='blue')

    if plot_vospec is True:
        axins.plot(wave2, flux2_der * wave2, '-',
                   color='gray', markersize=markersize, alpha=0.2)
        axins.plot(wave3, flux3_der * wave3, '-', label='IUE',
                   color='gray', markersize=markersize, alpha=0.5)
    if plot_vosa is True:

        missions = np.unique(catalogue)
        symbols = np.unique(catalogue)
        colors = phc.gradColor(np.arange(len(missions)), cmapn=cmap)
        colorsdict = dict(zip(missions, colors))
        list_mission = []
        for mission in catalogue:
            idx = np.where(catalogue == mission)
            colors = colorsdict[mission]

            if (mission in list_mission) is False:
                list_mission.append(mission)
                axins.errorbar(wave1[idx], flux1[idx] * wave1[idx],
                               yerr=Error1[idx],
                               fmt='.', label=mission.decode('UTF-8'),
                               color=colors, markersize=2 * markersize,
                               elinewidth=2)
            else:
                plt.errorbar(wave1[idx], flux1[idx] * wave1[idx],
                             yerr=Error1[idx],
                             fmt='.', color=colors, markersize=2 * markersize,
                             elinewidth=2)


# ==============================================================================
# Plot's limits
cursor = Cursor(ax, useblit=True, color='red', linewidth=2)

# Axis scales
ax.set_yscale("log")
ax.set_xscale('log')
# ax.set_xlim(1e3, 2e7)
ax.set_xlim(1e3, 1e10)
# ax.set_ylim(1e-16, 1e-4)
ax.set_ylim(1e-25, 1e-4)

# Axis label
ax.set_ylabel('$\lambda \, F_\lambda$ [erg/cm$^2$/s]')
ax.set_xlabel('$\lambda$ [$\AA$]')

# Extras
ax.minorticks_on()
ax.autoscale(enable=autoscale, axis=u'both', tight=False)
ax.legend(ncol=1, shadow=False, title="", fancybox=True, prop={'size': 7.9},
          numpoints=1, bbox_to_anchor=(1.02, 1.01), loc=2, fontsize='x-large')

plt.rc('font', family='serif', size=13)
plt.subplots_adjust(right=1.2)

# Saving
figure_name = star_name + '.pdf'
# plt.tight_layout()
plt.savefig(figure_name)

# Showing
plt.grid(True)
plt.show()
plt.close()


