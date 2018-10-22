import pyfits
import matplotlib.pylab as plt
import numpy as np
from astroquery.simbad import Simbad
import aplpy
from astropy import units as u
import math


# ------------------------------------------------------------------------------

obj = 'HR6510'
angle_deg = 0.004  # degrees
file = 'VISIR.2004-05-06T08_43_58.562.fits'
file = 'VISIR.2004-05-06T08:39:11.938.fits'
file = 'VISIR.2004-05-06T08:39:11.938.fits'
common_folder = '/Users/tangodown/Dropbox/2_Artigos/tex_aara/spectroscopy/'
folder_results = common_folder + 'results/'
hdulist = pyfits.open(file)
data = hdulist[0]

image = data.data

header = data.header
mjd = header['MJD-OBS']

ra = header['RA']
dec = header['DEC']

# ------------------------------------------------------------------------------
# image = image[0] + image[1] + image[2] + image[3] + image[4] + image[5] + image[6]
# fig = aplpy.FITSFigure(file)
fig = aplpy.FITSFigure(file, slices=[4])
fig.show_grayscale()

median_fig = np.median(fig.image.get_array())
median_fig = median_fig.data

fig.show_grayscale(vmin=median_fig, stretch=scale)

# ------------------------------------------------------------------------------
# Colorbar
# fig.add_colorbar()
# fig.colorbar.show()
# fig.colorbar.set_location('right')
# fig.colorbar.set_width(0.2)  # arbitrary units, default is 0.2
# fig.colorbar.set_pad(0.05)   # arbitrary units, default is 0.05
# fig.colorbar.set_font(size='medium', weight='medium', stretch='normal',
#                       family='sans-serif', style='normal',
#                       variant='normal')
# fig.colorbar.set_axis_label_text('Counts')
# fig.colorbar.set_axis_label_font(size=12, weight='bold')
# fig.colorbar.set_axis_label_pad(10)

# ------------------------------------------------------------------------------
# Coordinates
fig.set_yaxis_coord_type('longitude')

# Add a scalebar to indicate the size of the image (for 0.01 degrees)
fig.add_scalebar(angle_deg, str(angle_deg) + " degrees", color='black',
                 corner='top right')

fig.scalebar

# SIMBAD query
# result_table = Simbad.query_object(obj)
customSimbad = Simbad()
customSimbad.get_votable_fields()
customSimbad.add_votable_fields('plx', 'rot')
customSimbad.get_votable_fields()
result_table = customSimbad.query_object(obj)
parallax = result_table['PLX_VALUE']

distance_pc = 1000. / parallax
pc_to_au = u.pc.to(u.au)
distance_au = distance_pc * pc_to_au
angle_rad = math.radians(angle_deg)
scale = float(distance_au * angle_rad)

fig.scalebar.set_label("{:.2f}".format(scale) + " AU")
# ------------------------------------------------------------------------------
# Change the grid spacing (to be 0.1 degree)
# fig.ticks.set_xspacing("auto")
# fig.ticks.set_yspacing("auto")

# Change the formating of the tick labels
fig.tick_labels.set_xformat('hh:mm:ss')
fig.tick_labels.set_yformat('dd:mm:ss')

# Font size and ticks
# fig.axis_labels.set_font(size='large')
# fig.tick_labels.set_font(size='large')

# ------------------------------------------------------------------------------
# Add a grid over the image
fig.add_grid()
fig.grid.set_color('black')
fig.grid.set_alpha(0.3)
fig.grid.set_linewidth(0.4)
fig.grid.x_auto_spacing = True
fig.grid.y_auto_spacing = True


# ------------------------------------------------------------------------------
# Add a marker to indicate the position
fig.show_markers(ra, dec, layer='markers', edgecolor='white',
                 facecolor='none', marker='o', s=10, alpha=0.5)

# We can plot an array of ra,dec = [...],[...]
# All arguments of the method scatter() from matplotlib can be used

# Add a label to indicate the location
fig.add_label(ra + 0.0055, dec + 0.0055, obj, layer='source',
              color='black')

# ------------------------------------------------------------------------------

fig.show_contour(cmap='Blues', smooth=3)
fig.show_arrows(x=ra, y=dec, dx=0.005, dy=0.005)

# ------------------------------------------------------------------------------
# Use a present theme for publication
fig.set_theme('publication')  # Or 'pretty'  - for screen visualisation

# Saving...
plt.tight_layout()
fig.save(folder_results + str(obj) + '_' + str(mjd) + '_img.png')

# ------------------------------------------------------------------------------
zoom = True
if zoom is True:
    plt.tight_layout()
    fig.recenter(ra, dec, width=0.02, height=0.02)
    fig.save(folder_results + '/' + str(obj) + '_' + str(mjd) +
             '_img_zoom.png')

fig.close()
# ------------------------------------------------------------------------------
