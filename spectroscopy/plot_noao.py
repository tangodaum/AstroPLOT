import pyfits
import matplotlib.pylab as plt
import numpy as np

# file = 'ct2686868.fits.fz'
# file = 'c09i_080909_231129_ori.fits.fz'
# hdulist = pyfits.open(file)

# h1 = hdulist[1]
# data = h1.compressed_data

# image = []

# for i in range(len(data)):
# 	image.append(data[i])


# ==============================================================================
# VISIR

import pyfits
import numpy as np
from glob import glob
import matplotlib.pylab as plt
import time
files = glob('*.fits')


def read_visir(filename):
    try:
        hdulist = pyfits.open(filename)
        data = hdulist[0].data
        header = hdulist[0].header
        mjd = header['MJD-OBS']
        print(np.shape(data))
        plt.clf()
        plt.imshow(data[0])
        plt.show()
        time.sleep(10)
        return data, header, mjd
    except:
        print('ERRO:')
        print(filename)


for i in range(len(files)):
    try:
        data, header, mjd = read_visir(files[i])
    except:
        print('ERRO')



# =============================================================
# Exemplo de conversão de datas
import pyfits
from astropy.time import Time

file = 'IUE_SWP55941HL_lin.fits'

hdulist = pyfits.open(file)

data = hdulist[0].data
header = hdulist[0].header
date = header['DATE']
# substituir 't' por ' '

t = Time(date, format='iso', scale='utc')

# Para converter
mjd = t.mjd




