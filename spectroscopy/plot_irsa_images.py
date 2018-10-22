from astropy.io import fits
import pyfits
from matplotlib.colors import LogNorm
import os
import matplotlib.pyplot as plt
import pyhdust.phc as phc
import numpy as np


# ==============================================================================

common_folder = '/Users/tangodown/Dropbox/2_Artigos/tex_aara/spectroscopy/'
folder_results = common_folder + 'results/'
folder_table = common_folder + 'tables/'
folder_data = common_folder + 'dados_aara/irsa'
folder_data = common_folder + 'dados_aara/hst'
cmap = 'inferno'


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
create_list_files('list_files_irsa', folder=folder_data,
                  folder_table=folder_table)

fits_list = read_list_files_all(table_name='list_files_irsa.txt',
                                folder_table=folder_table)

dsstore = common_folder + 'dados_aara/irsa/.DS_Store'
# plt.clf()
for i in range(len(fits_list)):
    image_file = fits_list[i]
    print(image_file)
    if image_file != dsstore:

        hdu_list = pyfits.open(image_file)

        # Reading the header
        header = hdu_list[0].header
        image_data = hdu_list[0].data[0]
        hdu_list = fits.open(image_file)
        hdu_list.info()
        plt.imshow(image_data, cmap=cmap, norm=LogNorm(), alpha=0.2)
# plt.colorbar()
plt.savefig(folder_results + 'irsa_image.png')
plt.show()

