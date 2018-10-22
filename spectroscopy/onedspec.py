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

import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import argparse
import os

class SinglePointCalibration:
    def __init__(self, reference, angstrom_per_pixel):
        self.reference = reference
        self._angstrom_per_pixel = angstrom_per_pixel

    def angstrom(self, pixel):
        return (self.angstrom_per_pixel() * (pixel - self.reference.pixel)) + self.reference.angstrom

    def angstrom_per_pixel(self):
        return self._angstrom_per_pixel

    def intercept(self):
        return self.angstrom(0)

    def coefficient(self):
        return self.angstrom_per_pixel()



class CalibrationReference:

    def __init__(self, pixel, angstrom):
        self.pixel = pixel
        self.angstrom = angstrom

    @classmethod
    def from_string(cls, string):
        pixel, angstrom = string.split(':')

        # This could be a string, e.g. Ha, Hb or just a number.
        angstrom = cls.get_angstrom(angstrom)

        return cls(float(pixel), float(angstrom))

    @classmethod
    def get_angstrom(cls, angstrom):
        if angstrom in ElementLine.presets():
            return ElementLine.presets()[angstrom].angstrom
        else:
            return angstrom

    def __repr__(self):
        return "<Calibration reference: pixel: %d, angstrom: %f>" % (
            self.pixel, self.angstrom
        )

class Plotable:

    can_plot_image = False
    grayscale = False
    linestyle = '-'

    def plot_onto(self, axes, offset = 0):
        plot_args = {'label': self.label, 'linestyle': self.linestyle}

        if self.grayscale:
          plot_args['color'] = 'k'

        data = self.data() + offset

        axes.plot(self.wavelengths(), data, **plot_args)

    def max(self):
        return self.data().max()

    def divide_by(self, other_spectra):
        divided = np.divide(self.data(), other_spectra.interpolate_to(self))
        divided[divided==np.inf]=0
        return divided

    def interpolate_to(self, spectra):
        return np.interp(spectra.wavelengths(), self.wavelengths(), self.data())


class OneDSpec(Plotable):

    label = 'Raw spectra'
    label_header = 'DATE-OBS'

    def __init__(self, hdulist, calibration = False):
        self.hdulist = hdulist
        if calibration:
            self.calibration = calibration
        else:
            self.calibration = SinglePointCalibration(
                CalibrationReference(
                    self.get_header('CRPIX1'), self.get_header('CRVAL1')
                ),
                self.get_header('CDELT1')
            )
        self.set_label()

    def set_label(self):
        if self.label_header in self.header():
            self.label = self.get_header(self.label_header)

    def set_label_header(self, label_header):
        self.label_header = label_header
        self.set_label()

    def plot_image_onto(self, axes):
        return False

    def header(self):
        return self.hdulist[0].header

    def get_header(self, header):
        return self.header()[header]

    def wavelengths(self):
        intercept = self.calibration.intercept()
        coefficient = self.calibration.coefficient()
        return np.arange(len(self.data())) * coefficient + intercept

    def data(self):
        return self.hdulist[0].data

