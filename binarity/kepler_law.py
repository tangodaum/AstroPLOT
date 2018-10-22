#!/usr/bin/env python

from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import *
# import matplotlib.gridspec as gridspec
from astropy import constants
from astropy import units as u
import matplotlib.ticker as ticker

sfmt = ticker.ScalarFormatter(useMathText=True)
font = {'family': 'Times New Roman', 'size': 12}
plt.rcParams["font.family"] = "serif"

# ==============================================================================
# Constants
G = constants.G.si  # SI
M_sun = constants.M_sun.si  # SI
R_sun = constants.R_sun.si  # SI
L_sun = constants.L_sun.si  # SI


# ==============================================================================
def period_kepler(a, M_1, M_2):
    '''
    Calculates the period in days
    M_1: Primary mass
    M_2: Secondary mass
    '''

    P = np.sqrt(((4 * np.pi**2) /
                 (G * M_sun * (M_1 + M_2))) * (a * R_sun) ** 3)
    P = P.to(u.day)
    P = P.value

    return P


# ==============================================================================
# name the output file
pdfname = 'kepler_period.png'

# Parameters
M_1 = 7.28
M_2 = 2.0
M_3 = 0.5
rang = 250

R_pol = 3.70
despina_corr = 1.30
R_T = 19.29 * despina_corr * R_pol
err_R_T = (2.28 + 1.66) / 2
err_period = err_R_T * despina_corr * R_pol
R_T_mins_sig = (R_T - err_R_T) * despina_corr * R_pol

P_true = period_kepler(a=R_T, M_1=M_1, M_2=M_2)
P_true_minus_sigma = period_kepler(a=R_T - err_R_T, M_1=M_1, M_2=M_2)
P_true_plus_sigma = period_kepler(a=R_T + err_R_T, M_1=M_1, M_2=M_2)

# Period and semi-major axis arrays
y = np.arange(rang) * R_pol
x = period_kepler(a=y, M_1=M_1, M_2=M_2)  # days
y = y / R_pol

# Period and semi-major axis arrays
y2 = np.arange(rang) * R_pol
x2 = period_kepler(a=y2, M_1=M_1, M_2=M_3)  # days
y2 = y2 / R_pol

# create the main figure
fig = plt.figure(figsize=(8, 6), dpi=120)
# fig.subplots_adjust(hspace=0.05, wspace=0.0001, bottom=0.1,
#                     top=0.96, left=0.2, right=1.5)

# define two subplots with different sizes
# gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[3, 1])

# upper plot
# upperplot = plt.subplot(gs[0])
plot(x, y, '-', markeredgecolor="k", markerfacecolor="w",
     label=r'$m = 2.0 \, \mathrm{M}_\odot$', color='black',
     linewidth=2)
plot(x2, y2, '--', markeredgecolor="k", markerfacecolor="w",
     label=r'$m = 0.5 \, \mathrm{M}_\odot$', color='black',
     linewidth=2)

errorbar(P_true, R_T / R_pol, yerr=err_R_T, xerr=err_period,
         fmt='ko', label=r'$P = {0:0.2f} \pm {1:0.2f}$ days'.
         format(P_true, err_period))

errorbar(P_true, R_T / R_pol, yerr=err_R_T, xerr=err_period,
         fmt='ko', label=r'$a = {0:0.2f} \pm {1:0.2f} \, R_\star$'.
         format(R_T / R_pol, err_R_T))

plot([-0.05, 1000.], [-7.826, -7.826], 'k:')
ylabel(r"Semi-major axis $\, [R_\star]$", fontsize=14)
xlim([-0.05, 300.])
ylim([0, 120])
# plt.setp(upperplot.get_xticklabels(), visible=False)
plt.legend(loc='best', numpoints=1)
# plt.yscale('log')
# plt.xscale('log')

# # lower plot
# lowerplot = plt.subplot(gs[2])
# plot(x, y, 'ko', markeredgecolor="k", markerfacecolor="w")
# # errorbar(x, y, err, fmt='ko')
# plot([-0.05, 1.05], [0., 0.], 'k:')
# ylabel("Orbital Phase")
xlabel("Period [days]", fontsize=14)
# xlim([-0.05, 100.])
# ylim([0, 50])
minorticks_on()
# close and save file

savefig(pdfname)
plt.show()
# clf()
