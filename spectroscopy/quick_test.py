# ==============================================================================
# Heros
plt.clf()
cor = phc.gradColor(np.arange(len(fitslist_red)), cmapn=cmap)
for i in range(len(fitslist_red)):
    try:
        hdr = pyfits.getheader(fitslist_red[i])
        data0 = pyfits.getdata(fitslist_red[i])
        x0 = hdr['CDELT1'] * (np.arange(len(data0)) -
                              hdr['CRPIX1']) + hdr['CRVAL1']
        x0 = np.array(x0)
        x0 = x0 * 1e-4
        mjd = getMJD(hdr)
        ind = np.where((x0 < cont1) & (x0 > cont0))
        data = data0[ind]
        x = x0[ind]

        f = spt.linfit(x, data)
        # f = doNorm(l, f)
        plt.plot(x, data, color=cor[i])
    except:
         pass

plt.show()


# ==============================================================================
plt.clf()
cor = phc.gradColor(np.arange(len(uves_red)), cmapn=cmap)
for i in range(len(uves_red)):
    try:
         # print(uves[i])
         l, f, m = read_uves(uves_red[i])
         ind = np.where((l < 0.6590) & (l > 0.6530))
         
         f = f[ind]
         l = l[ind]
         # f = spt.linfit(l, f)
         # print(m, max(f))
         # print(f, l)
         if np.max(f) > 208700:
             f = doNorm(l, f)
             # print(m)
             plt.plot(l, f, color=cor[i], label=m)
    except:
         pass
plt.legend()
plt.show()


# ==============================================================================
plt.clf()
cor = phc.gradColor(np.arange(len(bess_data)), cmapn=cmap)
for i in range(len(bess_data)):
    try:
         l, f, m = read_fits(bess_data[i])
         ind = np.where((l < 0.664) & (l > 0.650))
         f = f[ind]
         l = l[ind]
         f = spt.linfit(l, f)
         # f = doNorm(l, f)
         print(m)
         plt.plot(l, f, color=cor[i], label=bess_data[i])
    except:
         pass
plt.legend()
plt.show()


# ==============================================================================
plt.clf()
cor = phc.gradColor(np.arange(len(feros)), cmapn=cmap)
for i in range(len(feros)):
    try:
         # print(feros[i])
         l, f, m = read_feros(feros[i])
         # ind = np.where((l < cont1) & (l > cont0))
         ind = np.where((l < 0.660) & (l > 0.653))
         f = f[ind]
         l = l[ind]
         # print(f, l)
         f = spt.linfit(l, f)
         # f = doNorm(l, f)
         plt.plot(l, f, color=cor[i])
    except:
         pass

plt.show()


# ==============================================================================
plt.clf()
for i in range(len(harps)):
    try:
         l, f, m = read_harps(harps[i])
         ind = np.where((l < cont1) & (l > cont0))
         ind = np.where((l < 0.660) & (l > 0.652))
         f = f[ind]
         l = l[ind]
         f = spt.linfit(l, f)
         plt.plot(l, f, color=cor[i])
    except:
         pass

plt.show()


# ==============================================================================
plt.clf()
cor = phc.gradColor(np.arange(len(opd)), cmapn=cmap)
for i in range(len(opd)):
    try:
         l, f, m = read_opd(opd[i])
         # ind = np.where((l < 0.6575) & (l > 0.6545))
         ind = np.where((l < 0.6575) & (l > 0.6552))
         ind = np.where((l < cont1) & (l > cont0))
         # print(l, f)
         f = f[ind]
         l = l[ind]
         # print(f, l)
         f = spt.linfit(l, f)
         # f = doNorm(l, f)
         # if i == 7:
         #     print(opd[i])
         #     plt.plot(l, f, color=cor[i], label=opd[i][-30:-1])
         plt.plot(l, f, color=cor[i], label=opd[i][-30:-1])
    except:
         pass
plt.hlines(y=1, xmin=cont0, xmax=cont1, linestyle='--', alpha=0.5)
plt.legend()
plt.show()



# ==============================================================================
plt.clf()
cor = phc.gradColor(np.arange(len(hanuschik)), cmapn=cmap)
for i in range(len(hanuschik)):
    try:
         l, f, m = read_hanushick(hanuschik[i])
         ind = np.where((l < cont1) & (l > cont0))
         f = f[ind]
         l = l[ind]
         f = spt.linfit(l, f)
         plt.plot(l, f, color=cor[i], label=m)
         print(m, hanuschik[i])
    except:
         pass
plt.legend()
plt.show()


# ==============================================================================
plt.clf()
cor = phc.gradColor(np.arange(len(dachs)), cmapn=cmap)
for i in range(len(dachs)):
    try:
         l, f, m = read_dachs(dachs[i])
         ind = np.where((l < cont1) & (l > cont0))
         f = f[ind]
         l = l[ind]
         f = spt.linfit(l, f)
         plt.plot(l, f, color=cor[i], label=m)
         print(m, dachs[i])
    except:
         pass
plt.legend()
plt.show()


# ==============================================================================
plt.clf()
for i in range(len(banerjee)):
    try:
         l, f, m = read_banerjee(banerjee[i])
         ind = np.where((l < cont1) & (l > cont0))
         f = f[ind]
         l = l[ind]
         f = spt.linfit(l, f)
         plt.plot(l, f, color=cor[i])
         print(m)
    except:
         pass

plt.show()

# ==============================================================================

plt.clf()
for i in range(len(pucheros_data)):
    try:
         l, f, m = read_pucheros(pucheros_data[i])
         ind = np.where((l < cont1) & (l > cont0))
         f = f[ind]
         l = l[ind]
         f = spt.linfit(l, f)
         plt.plot(l, f, color=cor[i])
         print(m)
    except:
         pass

plt.show()


# ==============================================================================

plt.clf()
cor = phc.gradColor(np.arange(len(nelson_feros)), cmapn=cmap)
for i in range(len(nelson_feros)):
    try:
        print(nelson_feros[i])
        l, f, m = read_feros_nelson(nelson_feros[i])
        ind = np.where((l < cont1) & (l > cont0))
        f = f[ind]
        l = l[ind]
        f = spt.linfit(l, f)
        # print(len(l), len(f), len(m))
        # print(l, f)
        plt.plot(l, f, color=cor[i], label=m)
        # plt.show()
    except:
         continue
# plt.legend()
# plt.autoscale()
# plt.vlines(x=0.65628, ymin=0.94, ymax=1.02, linestyle='--')
plt.show()


# ==============================================================================
def read_opd(fname):
    # file = 'aara.halpha.fits'
    # print(fname)
    spec = fitsio.read(fname)
    header = fitsio.read_header(fname)
    lbda = np.arange(len(spec)) * header['CDELT1'] + header['CRVAL1']
    if ('JD' in header) is False:
        date = header['DATE-OBS']
        date = date[0:10]
        fmt = '%Y-%m-%d'
        mjd = dt.datetime.strptime(date, fmt)
        mjd = date_to_jd(mjd.year, mjd.month, mjd.day)
        mjd = jd_to_mjd(mjd)
    else:
        jd = header['JD']
        mjd = jd - 2400000.5
    to_u = 1e-4
    lbda = to_u * lbda
    OBJ = header['OBJECT']
    # print(OBJ, header['INSTRUME'], mjd, )

    return lbda, spec, mjd


plt.clf()
cor = phc.gradColor(np.arange(len(nelson_lna)), cmapn=cmap)
for i in range(len(nelson_lna)):
    try:
        print(nelson_lna[i]) 
        l, f, m = read_opd(nelson_lna[i])
        # ind = np.where((l < cont1) & (l > cont0))
        # f = f[ind]
        # l = l[ind]
        # f = spt.linfit(l, f)
        plt.show()
        plt.plot(l, f, color=cor[i], label=m)
# print(m, nelson_lna[i])
    except:
         pass
# plt.legend()
# plt.vlines(x=0.65628, ymin=0.94, ymax=1.02, linestyle='--')
# plt.show()
