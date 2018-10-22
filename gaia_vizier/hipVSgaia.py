# ==============================================================================
# !/usr/bin/env python
# -*- coding:utf-8 -*-

# Created by B. Mota 2016-02-16 to present...

# import packages

import matplotlib.pyplot as plt
import matplotlib as mpl
# import matplotlib.font_manager as fm
import numpy as np
import pyhdust.phc as phc
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import astropy.coordinates as coord
from gaia.tap import cone_search
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import pylab
from astropy.coordinates import Angle
# font = fm.FontProperties(size=16)
font = {'family': 'serif', 'size': 20}

# pylab.ticklabel_format(axis='y', style='sci', scilimits=(1, 4))
pylab.tick_params(direction='in', size=3, which='both')
mpl.rc('xtick', labelsize=18)
mpl.rc('ytick', labelsize=18)
mpl.rc('font', **font)
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
fontsize_label = 10  # 'x-large'


__version__ = "0.0.1"
__author__ = "Bruno Mota"


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
# Parameters that must be defined
num_spa = 75
PyGaia = False
log_scale = False
dist_in_pc = False
# tag = 'hd_be_beacon+bsc'
# tag = 'bstars_bsc'
# tag = 'simbad_be'
tag = 'article'

if tag == 'bstars_bsc':
    stars = read_txt(table_name='tables/hd_b_bsc.txt')
    angle = 0.0005
if tag == 'hd_be_beacon+bsc':
    stars = read_txt(table_name='tables/hd_be_beacon+bsc.txt')
    angle = 0.0005
if tag == 'simbad_be':
    stars = read_txt(table_name='tables/simbad_be.txt')
    stars = stars.readlines()
    angle = 0.0005
    for i in range(len(stars)):
        temp = stars[i]
        stars[i] = temp.replace('\n', "")
if tag == 'article':
    stars = open('tables/list_iue_article.txt')
    # stars = open('tables/list_beacon_comp.txt')

    stars = stars.readlines()
    angle = 0.0005
    for i in range(len(stars)):
        temp = stars[i]
        stars[i] = temp.replace('\n', "")

print(num_spa * '=')
print('\nComparing Hipparcos VS Gaia\n')
print(num_spa * '=')


# ==============================================================================
customSimbad = Simbad()
customSimbad.TIMEOUT = 2000  # sets the timeout to 2000s
customSimbad.add_votable_fields('plx', 'plx_error', 'sptype', 'flux(V)')
customSimbad.get_votable_fields()

cmap = 'rainbow'
count = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
cor = phc.gradColor(count, cmapn=cmap, vmin=0, vmax=10)


# ==============================================================================
def main():
    global plx_gaia, eplx_gaia, plx_hip, eplx_hip, vmag_arr,\
        ratio_arr, dupsource
    plx_gaia, eplx_gaia = [], []
    plx_hip, eplx_hip = [], []
    vmag_arr, ratio_arr = [], []
    dupsource = []

    list_size = len(stars)
    for i in range(len(stars)):
        if tag == 'simbad_be':
            hd = str(stars[i])
        else:
            hd = 'HD' + str(int(stars[i]))
        # print(i / len(stars) * 100, '\%', )

        if hd == 'HD158427':
            color = 'lightgreen'
            alpha = 0.8
        elif hd == 'HD37795':
            color = 'red'
            alpha = 0.8
        else:
            color = 'blue'
            alpha = 0.2

        if PyGaia is True:
            coordinates = coord.SkyCoord.from_name(hd)
            ra = coordinates.ra.value
            dec = coordinates.dec.value

            pos = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree),
                           frame='icrs')

            width = u.Quantity(0.03, u.deg)
            height = u.Quantity(0.03, u.deg)
            r = Gaia.query_object_async(coordinate=pos, width=width,
                                        height=height)
            gaia_plx = r[0]['parallax']
            gaia_splx = r[0]['parallax_error']
            comp = r[0]['duplicated_source']
        else:
            r = Vizier.query_object(hd, radius=Angle(angle, "deg"),
                                    catalog='I/345/gaia2')
            try:
                table = r.values()[0]
                gaia_plx = table[0]['Plx']
                gaia_splx = table[0]['e_Plx']
                comp = table[0]['Dup']
            except:
                gaia_plx = np.nan
                gaia_splx = np.nan
                comp = 0

        if comp == 1:
            fmt = 's'
        else:
            fmt = 'o'

# ------------------------------------------------------------------------------
        # HIpparcos
        r = customSimbad.query_object(hd)
        hip_plx = r['PLX_VALUE']
        hip_plx = hip_plx.item()
        hip_splx = r['PLX_ERROR']
        hip_splx = hip_splx.item()
        vmag = r['FLUX_V'].item()
        sptyp = r['SP_TYPE']
        sptyp = sptyp.item()

        if 'O' in str(sptyp):
            cor_temp = cor[0]
        elif 'B0' in str(sptyp):
            cor_temp = cor[1]
        elif 'B1' in str(sptyp):
            cor_temp = cor[2]
        elif 'B2' in str(sptyp):
            cor_temp = cor[3]
        elif 'B3' in str(sptyp):
            cor_temp = cor[4]
        elif 'B4' in str(sptyp):
            cor_temp = cor[5]
        elif 'B5' in str(sptyp):
            cor_temp = cor[6]
        elif 'B6' in str(sptyp):
            cor_temp = cor[7]
        elif 'B7' in str(sptyp):
            cor_temp = cor[8]
        elif 'B8' in str(sptyp):
            cor_temp = cor[9]
        elif 'B9' in str(sptyp):
            cor_temp = cor[10]
        elif 'A0' in str(sptyp):
            cor_temp = cor[11]

        percent = i / list_size * 100

        print('{0:0.2f}, {1}, {2}, {3:0.2f}'.format(percent, hd,
                                                    sptyp, gaia_plx))

        if dist_in_pc is True:
            x = 1e3 / hip_plx
            y = 1e3 / gaia_plx
            ex = 1e3 * hip_splx / hip_plx**2
            ey = 1e3 * gaia_splx / gaia_plx**2

            plt.errorbar(x=x, y=y, xerr=ex,
                         yerr=ey, color=cor_temp, markersize=10,
                         alpha=alpha, fmt=fmt, capsize=4, capthick=2,
                         ecolor='black', elinewidth=2, markeredgewidth=2)
        else:
            plt.errorbar(x=hip_plx, y=gaia_plx, xerr=hip_splx,
                         yerr=gaia_splx, color=cor_temp, markersize=10,
                         alpha=alpha, fmt=fmt, capsize=4, capthick=2,
                         ecolor='black', elinewidth=2, markeredgewidth=2)

        plx_gaia.append(gaia_plx)
        eplx_gaia.append(gaia_splx)
        plx_hip.append(hip_plx)
        eplx_hip.append(hip_splx)
        vmag_arr.append(vmag)
        ratio_arr.append(gaia_plx / hip_plx)
        dupsource.append(comp)

    if dist_in_pc is True:
        plt.plot([0.1, 10000], [0.1, 10000], linestyle='--')
    else:
        plt.plot([0, 1000], [0, 1000], linestyle='--')
    # for i in range(len(plx_gaia)):
    #     plt.errorbar(x=plx_hip[i], y=plx_gaia[i],
    #                  xerr=eplx_hip[i], yerr=eplx_gaia[i],
    #                  color='blue', markersize=5, alpha=0.2, fmt='o',
    #                  ecolor='black', elinewidth=3, capsize=4,
    #                  capthick=2, markeredgewidth=1)

    norm = plt.Normalize(vmin=0, vmax=11)
    s_m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    colourbar = plt.colorbar(s_m, ticks=list(count))
    # colourbar.set_label('ST')
    colourbar.set_ticklabels(['O', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                              'B7', 'B8', 'B9', 'A0'])

    if tag == 'hd_be_beacon+bsc':
        if dist_in_pc is True:
            plt.xlim(10, 1e4)
            plt.ylim(10, 1e4)
        else:
            plt.xlim(-0.1, 25)
            plt.ylim(-0.1, 25)
    if tag == 'bstars_bsc':
        if dist_in_pc is True:
            plt.xlim(10, 1e4)
            plt.ylim(10, 1e4)
        else:
            plt.xlim(-0.1, 25)
            plt.ylim(-0.1, 25)
    if tag == 'simbad_be':
        if dist_in_pc is True:
            plt.xlim(10, 1e4)
            plt.ylim(10, 1e4)
        else:
            plt.xlim(-0.1, 25)
            plt.ylim(-0.1, 25)

    if tag == 'article':
        if dist_in_pc is True:
            plt.xlim(10, 1e4)
            plt.ylim(10, 1e4)
        else:
            plt.xlim(-0.1, 30)
            plt.ylim(-0.1, 30)

    if dist_in_pc is True:
        plt.xlabel(r'$d_\mathrm{Hipparcos} \mathrm{ [pc]}$')
        plt.ylabel(r'$d_\mathrm{Gaia} \mathrm{ [pc]}$')
    else:
        plt.xlabel(r'$\pi_\mathrm{Hipparcos} \mathrm{ [mas]}$')
        plt.ylabel(r'$\pi_\mathrm{Gaia} \mathrm{ [mas]}$')

    plt.minorticks_on()
    plt.tight_layout()

    if log_scale is True:
        plt.xscale('log')
        plt.yscale('log')

    plt.savefig('figures/' + 'hipVSgaia' + '_' + tag + '.pdf')
    plt.show()

    # Vmag vs Ratio
    plt.clf()
    for i in range(len(ratio_arr)):
        try:
            if dupsource[i] == 1:
                mark = 's'
            else:
                mark = 'o'

            hd = 'HD' + str(int(stars[i]))
            if hd == 'HD158427':
                color = 'lightgreen'
                alpha = 0.8
            elif hd == 'HD37795':
                color = 'red'
                alpha = 0.8
            else:
                color = 'blue'
                alpha = 0.2

            plt.scatter(vmag_arr[i], ratio_arr[i], color=color,
                        s=80, alpha=alpha, edgecolors='black', linewidths=2,
                        marker=mark)

        except:
            pass

    plt.xlabel(r'$V \mathrm{ [mag]}$')
    plt.ylabel(r'$\pi_\mathrm{Gaia} / \pi_\mathrm{Hipparcos}$')
    plt.minorticks_on()
    plt.tight_layout()
    plt.yscale('log')
    plt.xlim(0, 10)
    plt.ylim(0.1, 10)
    plt.hlines(xmin=0, xmax=10, y=1, linestyle='--')
    plt.savefig('figures/' + 'vmagVSratio' + '_' + tag +'.pdf')
    plt.show()


# ==============================================================================
if __name__ == '__main__':
    main()
