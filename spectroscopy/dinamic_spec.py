#! /usr/bin/env python

# Copyright 2015 Malcolm Locke
#
# dynaspec is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dynaspec is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with dynaspec.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import pyfits
import matplotlib.pyplot as plt
import argparse
from onedspec import OneDSpec
from scipy.interpolate import griddata
from scipy.constants import c
from matplotlib.ticker import ScalarFormatter


def convert_wavelengths_to_velocity(wavelengths, frequency):
    return [((w / frequency) - 1) * (c / 1000) for w in wavelengths]


def get_phase(timestamp, epoch, period):
    phase = ((timestamp - epoch) % period) / period
    print("%f: phase = %f" % (timestamp, phase))
    return phase


parser = argparse.ArgumentParser(description='Make 2D grayscale image from multiple 1D spectra')
parser.add_argument('master', type=str)
parser.add_argument('file', type=str, nargs='+')
parser.add_argument('--cmap', '-c', type=str, default='gray')
parser.add_argument('--title', '-t', type=str, default='')
parser.add_argument('--vmin', type=float)
parser.add_argument('--vmax', type=float)
parser.add_argument('--velocity', type=float,
                    help='Specify rest frequency for velocity plot')
parser.add_argument('--phase', type=str,
                    help='Generate phase plot.  Argument must be in the form epoch:period.  Uses JD-MID header')

args = parser.parse_args()

x = np.array([])
y = np.array([])
z = np.array([])

ticks = []
print "Processing %s" % (args.master)
master = OneDSpec(pyfits.open(args.master))
header = master.hdulist[0].header
epoch = None
period = None

if args.phase:
    epoch, period = args.phase.split(':')
    epoch = float(epoch)
    period = float(period)

length = len(y)
wavelengths = np.copy(lbdas)
if args.velocity:
    wavelengths = convert_wavelengths_to_velocity(wavelengths, args.velocity)
start_jd, _ = divmod(header['JD-MID'], 1)
a = np.empty(length)

if args.phase:
    phase = get_phase(header['JD-MID'], epoch, period)
    ticks.append(phase)
    a.fill(phase)
else:
    ticks.append(header['JD-MID'] - start_jd)
    a.fill(header['JD-MID'] - start_jd)

x = np.append(x, wavelengths)
y = np.append(y, a)
z = np.append(z, master.data())
xlabel = ''
if 'CUNIT1' in header:
    xlabel = header['CUNIT1']
if args.velocity:
    xlabel = 'Velocity (km/s)'

for f in args.file:
    print "Processing %s" % (f)
    s = OneDSpec(pyfits.open(f))
    a = np.empty(length)

    if args.phase:
        phase = get_phase(s.hdulist[0].header['JD-MID'], epoch, period)
        ticks.append(phase)
        a.fill(phase)
    else:
        ticks.append(s.hdulist[0].header['JD-MID'] - start_jd)
        a.fill(s.hdulist[0].header['JD-MID'] - start_jd)

    x = np.append(x, wavelengths)
    y = np.append(y, a)
    interp = s.interpolate_to(master)
    if args.velocity:
        convert_wavelengths_to_velocity(interp, args.velocity)
    z = np.append(z, interp)


# define regular grid spatially covering input data
n = length
xg = np.linspace(x.min(), x.max(), len(wavelengths))
yg = np.linspace(y.min(), y.max(), length)

X, Y = np.meshgrid(xg, yg)

# interpolate Z values on defined grid
Z = griddata(np.vstack((x.flatten(), y.flatten())).T,
             np.vstack(z.flatten()), (X, Y), method='linear').reshape(X.shape)
# mask nan values, so they will not appear on plot
Zm = np.ma.masked_where(np.isnan(Z), Z)

# plot
fig, ax = plt.subplots()
colormesh_args = {'cmap': args.cmap, 'shading': 'gouraud'}

if args.vmax:
    colormesh_args['vmax'] = args.vmax

if args.vmin:
    colormesh_args['vmin'] = args.vmin

pcm = ax.pcolormesh(X, Y, Zm, **colormesh_args)
fig.colorbar(pcm, ax=ax, use_gridspec=True)
ax.set_title(args.title)
ax.set_xlabel(xlabel)
if args.phase:
    ax.set_ylabel("Phase")
    ax.invert_yaxis()
else:
    ax.set_ylabel("JD - %d" % start_jd)
plt.yticks(ticks)

plt.tight_layout()
plt.show()
