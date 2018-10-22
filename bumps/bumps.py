#Este programa serve para realizar a leitura de uma tabela VO e para plotar a SED na
#regiao do bump de 2200Angs

#===============================================================================
import numpy as np
#import matplotlib.animation as animation #verificar 
#import astropy
import atpy
#import pylab
#from matplotlib.widgets import Cursor
from PyAstronomy import pyasl
from astropysics import obstools
from scipy.optimize import curve_fit
#import sympy as sym
import pyhdust.phc as phc
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
#from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset
#from matplotlib.ticker import MaxNLocator
#import seaborn

#=======================================================================
# import packages
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import *
import pyhdust.phc as phc
#import seaborn

#===============================================================================
#Basic data

star_name = 'rhooph'
# star_name = 'rhooph'
file_name = 'bump_'+star_name
xml_files = ['iue/bump2200_'+star_name+'.xml']

labels      = ['aara','dsco','acol','48lib','ztau','rhooph'] 
labels      = ['']*10
autoscale   = True
option      = 2   #Option 1: delta = f_adjust/f_spectra / option 2:  delta = f_adjust/f_adjust2
zoom        = False
colors      = ['black','red','green','blue','orange','darkmagenta','salmon',\
'darkturquoise','seagreen', 'MediumTurquoise','Aquamarine','Khaki',\
'SlateBlue','Coral','YellowGreen','LawnGreen','Peru','SeaGreen','Chocolate',\
'Aqua','Orchid']

#===============================================================================
def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax


#===============================================================================
#Use a vo format: xml.
fluxes      = []
waves       = []
inv_waves   = []
fluxesUnred = []
ff0s        = []
ff0s2       = []
rv          = 3.1
if star_name == 'acol':
    ebmv     = np.arange(0.00,0.1,0.02)
    zoom_ay  = 1.5
    axylim   = [0.005,0.1]
    ayylim   = [0,0]
if star_name == 'ztau':
    ebmv     = np.arange(0.00,0.1,0.02)
    zoom_ay  = 1.5
    axylim   = [0.005,0.1]
    ayylim   = [0,0]
if star_name == '48lib':
    ebmv     = np.arange(0.00,0.3,0.05)
    zoom_ay  = 1.5
    axylim   = [0.005,0.1]
    ayylim   = [0,0]
if star_name == 'dsco':
    ebmv     = np.arange(0.06,0.3,0.05)
    zoom_ay  = 1.5
    axylim   = [0.005,0.05]
    ayylim   = [3e-10,0.5e-7]
if star_name == 'aara':    
    ebmv     = np.arange(0.00,0.2,0.04)
    zoom_ay  = 1.5
    axylim = [0.005,0.05]
    ayylim = [0.5e-11,1e-8]
if star_name == 'rhooph':
    ebmv     = np.arange(0.0,1.2,0.2)
    zoom_ay  = 0.8
    axylim   = [0.005,0.5]
    ayylim   = [8e-12,1e-5]

#===============================================================================
for j in range(len(ebmv)):
    for i in range(len(xml_files)):
        t = atpy.Table(xml_files[i],tid=1)
        #Para acessar os nomes das colunas:
        t.columns
        #Define os arrays:
        flux = t['Flux0'][:] #erg/cm2/s/A
        wave = t['SpectralAxis0'][:] #Angstrom
        inv_wave   = (1.e4)/wave
        #Fitzpatrick
        fluxUnred  = pyasl.unred(wave, flux, ebv=ebmv[j], R_V=rv)    
        ff0  = flux/fluxUnred
        ff0 = ff0/np.max(ff0)
        #Cardelli
        extinction_law = obstools.CardelliExtinction(EBmV=ebmv[j], Rv=rv)
        Av = extinction_law(wave)
        m  = -2.5*np.log10(flux)
        m  = m-Av
        flux_der = 10.**(-m/2.5)
        ff02 = flux/flux_der
        ff02 = ff02/np.max(ff02)
        fluxes.append(flux)
        waves.append(wave)
        inv_waves.append(inv_wave)
        fluxesUnred.append(fluxUnred)
        ff0s.append(ff0)
        ff0s2.append(ff02)

#===============================================================================
#Model
p     = 1.
C0    = 1.
fs    = 0.1
#~ alfp  = 2.
Alamb = 1.
lambp = (waves[0]*1e-4)**p
#~ y2 = C0*np.e**((-1.*alfp*ebmv)/lambp)  #fs<<1
#~ y3 = C0*(1.+fs)*np.e**((-1.*alfp*ebmv)/lambp) #Scattered starlight

#Ploting the figure
fig = plt.figure(figsize=(12, 6))

#===============================================================================
#Left figure
for i in range(len(fluxes)):
    ax = fig.add_subplot(121, axisbg='white')
    #~ ax.plot(inv_waves[i], fluxes[i], '-',label=labels[i],color=colors[i])
    #~ ax.plot(inv_waves[i],fluxesUnred[i], '--',label=labels[i]+'dered',color=colors[i+1])
    ax.plot(inv_waves[i],ff0s[i], '-', label=('IUE Spectra'),color=colors[i])
    ax.plot(inv_waves[i],ff0s2[i], '--',label=('E(B-V) = %0.2f' % ebmv[i]),\
    color=colors[i])
    ax.annotate(('E(B-V) = %0.2f' % (ebmv[i])),xy=(inv_waves[i][0],\
    ff0s2[i][0]),xycoords='data',xytext=(inv_waves[i][0]+0.5,\
    ff0s2[i][0]),arrowprops=dict(arrowstyle="<-",\
    connectionstyle="angle,angleA=-90,angleB=180,rad=5"),size='small')
ax.axvline(x=1e4/2175,color='red',alpha=0.3)

#===============================================================================
#Right Figure
def func(x, a, b):
    return a*x**-4 + b

def func2(x, a, b, c, d, e):
    return a*x**-4 + b*x**-3 +c*x**-2 + d*x**-1+e

wave_lim_min  = 1900.
wave_lim_max  = 2020.
wave_lim_min2 = 2400.
wave_lim_max2 = 2550.
wave_lim_min3 = 1900.
wave_lim_max3 = 2550.


delta     = [] 
zoom_y    = [] 
y_bump    = []
x_bump    = []
z_bump    = []
y_adjust  = []
y_adjust2 = []

for i in range(len(fluxes)):
    ay = fig.add_subplot(122, axisbg='white')
    ay.plot(waves[i],fluxesUnred[i], '-',label=('E(B-V) = %0.2f' % ebmv[i]),\
    color=colors[i])
    new_lbdarr = []
    new_fluxes = []
    new_lbdarr2 = []
    new_fluxes2 = []
    for k in range(len(waves[i])):
        if waves[i][k]>wave_lim_min and waves[i][k]<=wave_lim_max or\
        waves[i][k]>=wave_lim_min2 and waves[i][k]<=wave_lim_max2:
            new_lbdarr.append(waves[i][k])
            new_fluxes.append(fluxesUnred[i][k])
    new_fluxes = np.array(new_fluxes)
    new_lbdarr = np.array(new_lbdarr)
    popt, pcov = curve_fit(func, new_lbdarr, new_fluxes)

    for k in range(len(waves[i])):
        if waves[i][k]>wave_lim_min3 and waves[i][k]<=wave_lim_max3:
            new_lbdarr2.append(waves[i][k])
            new_fluxes2.append(fluxesUnred[i][k])
    new_fluxes2 = np.array(new_fluxes2)
    new_lbdarr2 = np.array(new_lbdarr2)
    popt2, pcov2 = curve_fit(func2, new_lbdarr2, new_fluxes2)

    if option == 1:
        wave_bump = phc.find_nearest(waves[i], 2175)
        index  = np.where(waves[i]==wave_bump)
        if func(wave_bump,*popt)>=fluxesUnred[i][index]:
            Delta = func(wave_bump,*popt)/fluxesUnred[i][index]
        else:
            Delta = fluxesUnred[i][index]/func(wave_bump,*popt)

    if option == 2:    
        wave_bump = phc.find_nearest(waves[i], 2175)
        index = np.where(waves[i] == wave_bump)
        if func(wave_bump, *popt) >= func2(wave_bump, *popt2):
            Delta = func(wave_bump, *popt) / func2(wave_bump, *popt2)
        else:
            Delta = func2(wave_bump, *popt2) / func(wave_bump, *popt)

    #~ print(Delta, func(wave_bump,*popt),fluxesUnred[i][index])
    delta.append(Delta)
    if option == 1:
        ay.plot(wave_bump,func(wave_bump,*popt),'o',wave_bump,fluxesUnred[i][index],\
        's',color=colors[i])
        ay.plot(waves[i],func(waves[i],*popt), '--',color=colors[i])
        y_bump.append((func(wave_bump,*popt)))
        z_bump.append(fluxesUnred[i][index])
        x_bump.append(wave_bump)
    if option == 2:
        ay.plot(new_lbdarr2,func2(new_lbdarr2,*popt2), '--',color=colors[i])
        ay.plot(wave_bump,func(wave_bump,*popt),'o',wave_bump,func2(wave_bump,*popt2),\
        's',color=colors[i])
        ay.plot(waves[i],func(waves[i],*popt), '--',color=colors[i])
        y_bump.append((func(wave_bump,*popt)))
        z_bump.append((func2(wave_bump,*popt2)))
        x_bump.append(wave_bump)    
    zoom_y.append(fluxesUnred[i][index])
    y_adjust.append(func(waves[i],*popt))
    y_adjust2.append(func2(waves[i],*popt2))

ay.axvline(x=2175,color='red',alpha=0.3)    

for i in range(len(ebmv)):
    #Window: Delta vs EBmV
    rect = [0.65,0.65,0.3,0.3]
    ax1 = add_subplot_axes(ay,rect)
    ax1.plot(ebmv[i],delta[i],'o',color=colors[i])

#Zoom
if zoom == True:
    axins = zoomed_inset_axes(ay,zoom_ay, loc=3) # zoom = 2
    axins.set_xlim(wave_lim_min, wave_lim_max2)
    axins.xaxis.tick_top()
    #~ axins.set_ylim(np.min(zoom_y)-1e-9,  np.max(zoom_y)+1e-9)
    axins.set_ylim(np.min(zoom_y),  np.max(zoom_y)*1.3)
    mark_inset(ay, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    axins.xaxis.set_major_locator(MaxNLocator(nbins=1, prune='lower'))
    #~ axins.invert_yaxis()
    #~ axins.set_axis_bgcolor('none')

    
    for i in range(len(fluxesUnred)):    
        axins.plot(x_bump[i],y_bump[i], '--',color=colors[i])
        if option == 1:
            axins.plot(waves[i],fluxesUnred[i],color=colors[i])
            axins.plot(waves[i],y_adjust[i],'--',color=colors[i])
            axins.plot(x_bump[i],y_bump[i],'o',wave_bump,z_bump[i],\
            's',color=colors[i])
        if option == 2:
            axins.plot(waves[i],y_adjust[i],'--',color=colors[i])
            axins.plot(waves[i],y_adjust2[i],'-',color=colors[i])
            axins.plot(x_bump[i],y_bump[i],'o',wave_bump,z_bump[i],\
            's',color=colors[i])
    axins.set_yscale("log")

#===============================================================================
#Plot cursor
#cursor = Cursor(ax, useblit=True, color='red', linewidth=2)

#Log scale?
#ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(3,11.5)
ax.set_ylim(np.min(ff0s2)+axylim[0],np.max(ff0s2)+axylim[1])
ax1.set_xlim(min(ebmv)-0.02,max(ebmv)+0.02)
ay.set_ylim(np.min(fluxesUnred)-ayylim[0],np.max(fluxesUnred)+ayylim[1])
#~ ax1.set_yscale("log")
ay.set_yscale("log")



#Axis label
ax.set_ylabel('F/F$_0$')
ax.set_xlabel('$\mu$ m$^{-1}$')
ax1.set_ylabel('$\Delta$')
ax1.set_xlabel('E(B-V)')
ay.set_ylabel('$F_\lambda$ [erg/cm$^2$/s/$\AA$]')
ay.set_xlabel('Wavelength [$\AA$]')

#Plot's limits
#~ pylab.ylim([min(flux),max(flux)])
#~ pylab.xlim([min(inv_wave),max(inv_wave)])
#ax.autoscale(enable=autoscale, axis=u'both', tight=False)
#~ ay.autoscale(enable=autoscale, axis=u'both', tight=False)

#frescuras
ax.minorticks_on()
ax1.minorticks_on()
ay.minorticks_on()

if zoom == True:
    axins.minorticks_on()
    #~ ax1.set_yscale("log")
    axins.set_yticks([])
    axins.set_xticks([])

#ax.legend(loc=3, ncol=1, shadow=False, title="", fancybox=True,fontsize='small')
ay.legend(loc=4, ncol=2, shadow=False, title="", fancybox=False,fontsize='x-small',frameon=False)

#Saving
figure_name = file_name+'.pdf'
plt.savefig(figure_name)

#Showing
plt.grid(True)
plt.show()
