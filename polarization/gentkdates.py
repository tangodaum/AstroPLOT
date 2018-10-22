def gentkdates(mjd0, mjd1, fact, step, dtstart=None):
    """ Generates round dates between > mjd0 and < mjd1 in a given step.
    Valid steps are:

        'd/D/dd/DD' for days;
        'm/M/mm/MM' for months;
        'y/Y/yy/YY/yyyy/YYYY' for years.

    dtstart (optional) is expected to be in datetime.datetime.date() format
    [i.e., datetime.date(yyyy, m, d)]

    fact must be an integer
    """
    #check sanity of dtstart
    if dtstart is None:
        dtstart = dt.datetime(*jdcal.jd2gcal(jdcal.MJD_0,mjd0)[:3]).date()
        mjdst = jdcal.gcal2jd(dtstart.year,dtstart.month,dtstart.day)[1]
    else:
        mjdst = jdcal.gcal2jd(dtstart.year,dtstart.month,dtstart.day)[1]
        if mjdst < mjd0-1 or mjdst > mjd1:
            print('# Warning! Invalid "dtstart". Using mjd0.')
            dtstart = dt.datetime(*jdcal.jd2gcal(jdcal.MJD_0,mjd0)[:3]).date()
    #define step 'position' and vector:
    basedata = [dtstart.year, dtstart.month, dtstart.day]
    dates =  []
    mjd = mjdst
    if step.upper() in ['Y','YY','YYYY']:
        i = 0
        while mjd < mjd1+1:
            dates +=  [dt.datetime(*basedata).date()]
            basedata[i] += fact
            mjd = jdcal.gcal2jd(*basedata)[1]
    elif step.upper() in ['M','MM']:
        i = 1
        while mjd < mjd1+1:
            dates += [dt.datetime(*basedata).date()]
            basedata[i] += fact
            while basedata[i] > 12:
                basedata[0] += 1
                basedata[1] -= 12
            mjd = jdcal.gcal2jd(*basedata)[1]
    elif step.upper() in ['D','DD']:
        i = 2
        daysvec = np.arange(1,29,fact)
        if basedata[i] not in daysvec:
            j = 0
            while daysvec[j+1] < basedata[i]:
                j += 1
            daysvec += basedata[i]-daysvec[j]
            idx = np.where(daysvec < 29)
            daysvec = daysvec[idx]
        else: 
            j = np.where(daysvec == basedata[i])[0]
        while mjd < mjd1+1:
            dates += [dt.datetime(*basedata).date()]
            j += 1
            if j == len(daysvec):
                j = 0
                basedata[1] += 1
                if basedata[1] == 13:
                    basedata[1] = 1
                    basedata[0] += 1
            basedata[i] = daysvec[j]
            mjd = jdcal.gcal2jd(*basedata)[1]
    else:
        print('# ERROR! Invalid step')
        raise SystemExit(1)
    return dates
