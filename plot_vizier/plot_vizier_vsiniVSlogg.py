#######
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import numpy as np
import pyhdust.phc as phc
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import csv
import os
import re

Vizier.ROW_LIMIT = -1
# cat = ['J/MNRAS/361/1055/table1', 'J/A+A/393/887/SMC-Be']
cat = 'J/ApJ/722/605/table4'  # Huang 2010
catalogs = Vizier.get_catalogs(cat)
# catalog1 = catalogs[0]
# catalog2 = catalogs[1]
catalog1 = catalogs[0]

# Operating with the data
data = catalog1.as_array()

# Print available data
print(data.dtype)

# Taking data
vsini = data['vsini']
evsini = data['e_vsini']
teff = data['Teff']
eteff = data['e_Teff']
logg = data['log_g_']
elogg = data['e_log_g_']
hd = data['HD']

# Plotting
plt.errorbar(x=vsini, xerr=evsini, y=logg, yerr=elogg, linestyle='')
plt.show()

