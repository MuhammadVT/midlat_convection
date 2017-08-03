def cos_fit(ftype="fitacf", azbin_nvel_min=10, naz_min=3, az_span_min=30, baseLocation="./",
                   dbName=None, sqrt_weighting=True):
    
    """ Do the cosine fitting to the LOS data, and store the results in a different table
    named "master_cosfit". This table only have qualifying latc-lonc grid points
    azbin_nvel_min : int
        minimum number of measurements each azimuthal bin should have to be
        qualified for cosfitting. 
    naz_min : int
        the minimum number of azimuthal bins within each grid cell.
        cosine fitting is done if a grid cell has at least
        naz_min number of qualifying azimuthal bins
    az_span_min : int
        minimum azimuhtal span a grid cell should have to be qualified for cosfitting.
    sqrt_weighting : bool
        if set to True, the fitting is weighted by the number of points at each azimuthal bin
        if set to False, then all azimuthal bins are equla regardless of the nubmer of points
        they contain.
    
    """
    import sqlite3
    import datetime as dt
    import numpy as np 

    # make db connection for dopsearch
    if dbName is None:
        dbName = "master_" + ftype + ".sqlite"
    conn = sqlite3.connect(baseLocation + dbName)
    cur = conn.cursor()

    if sqrt_weighting:
        T1 = "master_cosfit"
    else:
        T1 = "master_cosfit_equal_weighting"
    
    #T1 = "master_cosfit"
    T2 = "master_summary"

#    # build a new master table where we store the cosfitted velocities at each lat-lon grid points
#    cur.execute("CREATE TABLE IF NOT EXISTS {tb1}\
#                (vel_mag REAL, vel_mag_err REAL, vel_dir REAL, vel_dir_err REAL,\
#                vel_count INTEGER, gazmc_count INTEGER,\
#                glatc REAL, glonc REAL,\
#                PRIMARY KEY (glatc, glonc))".format(tb1=T1))
#
#    # select the velocity data grouping by latc-lonc bins for the nigh side
#    command = "SELECT count(*), glatc, glonc FROM {tb2}\
#               WHERE (glonc BETWEEN 0 AND 135) OR (glonc BETWEEN 225 AND 360)\
#               GROUP BY glatc, glonc\
#               ".format(tb2=T2)

    cur.execute("CREATE TABLE IF NOT EXISTS {tb1}\
                (vel_mag REAL, vel_mag_err REAL, vel_dir REAL, vel_dir_err REAL,\
                vel_count INTEGER, gazmc_count INTEGER,\
                glatc REAL, glonc REAL,\
                PRIMARY KEY (glatc, glonc))".format(tb1=T1))

    # select the velocity data grouping by latc-lonc bins. Also each latc-lonc cell should
    #contain gazmc bins that have at least azbin_nvel_min amount of measurements
    command = "SELECT count(*), glatc, glonc, group_concat(gazmc) FROM\
               (SELECT * FROM {tb2} WHERE vel_count >= {azbin_nvel_min})\
               GROUP BY glatc, glonc\
               ".format(tb2=T2, azbin_nvel_min=azbin_nvel_min)
    #WHERE (glonc BETWEEN 0 AND 135) OR (glonc BETWEEN 225 AND 360)

    cur.execute(command)
    rws = cur.fetchall()

    # filter out lat-lon grid points that have less than 3 qualifying amimuthal bins 
    rws = [x for x in rws if x[0] >= naz_min]

    # filter out lat-lon grid points that have less than 30 degrees azimuthal span
    for rwi in rws:
        az_rwi = np.sort(np.array([int(x) for x in rwi[3].split(",")]))
        if len(az_rwi) == 3:
            if az_rwi.tolist()==[5, 345, 355] or az_rwi.tolist()==[5, 15, 355]:
                #print az_rwi
                rws.remove(rwi)
            elif az_rwi.tolist()==[az_rwi[0], az_rwi[0]+10, az_rwi[0]+20]:
                #print az_rwi
                rws.remove(rwi)
            else:
                continue
        else:
            continue

    azm_count = [x[0] for x in rws]
    lat = [x[1] for x in rws]
    lon = [x[2] for x in rws]

    for ii in xrange(len(lat)):
        command = "SELECT median_vel, vel_count, gazmc FROM {tb2}\
                   WHERE glatc={lat}\
                    AND glonc={lon}\
                    ORDER BY gazmc"\
                    .format(tb2=T2, lat=lat[ii], lon=lon[ii])
        cur.execute(command)
        rows = cur.fetchall()
        median_vel = np.array([x[0] for x in rows])
        vel_count = np.array([x[1] for x in rows])
        if sqrt_weighting:
            sigma =  1./np.sqrt(vel_count)
        else:
            sigma =  np.array([1.0 for x in rows])
        azm = np.array([x[2] for x in rows])

        # do cosine fitting with weight
        fitpars, perrs = cos_curve_fit(azm, median_vel, sigma)
        vel_mag = round(fitpars[0],2)
        vel_dir = round(np.rad2deg(fitpars[1]) % 360,1)
        vel_mag_err = round(perrs[0],2)
        vel_dir_err = round(np.rad2deg(perrs[1]) % 360, 1)

        # populate the table 
        command = "INSERT OR IGNORE INTO {tb1} (vel_mag,\
                    vel_mag_err, vel_dir, vel_dir_err, vel_count,\
                    gazmc_count, glatc, glonc) VALUES ({vel_mag},\
                    {vel_mag_err}, {vel_dir}, {vel_dir_err}, {vel_count},\
                    {gazmc_count}, {glatc}, {glonc})".format(tb1=T1, vel_mag=vel_mag,\
                    vel_mag_err=vel_mag_err, vel_dir=vel_dir,\
                    vel_dir_err=vel_dir_err, vel_count=np.sum(vel_count),\
                    gazmc_count =azm_count[ii], glatc=lat[ii], glonc=lon[ii])
        cur.execute(command)
        print "finish inserting cosfit result at " + str((lat[ii], lon[ii]))

    # commit the change
    conn.commit()

    # close db connection
    conn.close()

    return

def cosfunc(x, Amp, phi):
    import numpy as np
    return Amp * np.cos(1 * x - phi)

def cos_curve_fit(azms, vels, sigma):
    import numpy as np
    from scipy.optimize import curve_fit
    fitpars, covmat = curve_fit(cosfunc, np.deg2rad(azms), vels, sigma=sigma)
    perrs = np.sqrt(np.diag(covmat)) 

    return fitpars, perrs

def worker(baseLocation, dbName):

    import datetime as dt

    # input parameters
    azbin_nvel_min=10
    naz_min=3
    az_span_min=30

    ftype = "fitacf"
    #ftype = "fitex"
    sqrt_weighting=True
    # do the cosine fitting to the grid velocities
    print "doing cosine fitting to each of the grid cell velocities"
    cos_fit(ftype=ftype, azbin_nvel_min=azbin_nvel_min, naz_min=naz_min,
            az_span_min=az_span_min, baseLocation=baseLocation,
            dbName=dbName, sqrt_weighting=sqrt_weighting)
    print "finished cosine fitting"


    return
if __name__ == "__main__":

    import multiprocessing
    
    ftype = "fitacf"
    seasons = ["winter", "summer", "equinox"]
    #seasons = ["winter"]
    binned_season = False
    binned_F107 = False
    binned_imf = True

    # list of dbNames to be cosine fitted
    if binned_season:
        dbName_list = [None]


    if binned_F107:
        dbName_list = []
        #f107_bins = [[0, 100], [100, 175], [175, 500]]
        #f107_bins = [[0, 100], [100, 150], [150, 500]]
        #f107_bins = [[0, 105], [105, 125], [125, 500]]
        #f107_bins = [[0, 110], [110, 500]]
        #f107_bins = [[0, 120], [120, 500]]
        #f107_bins = [[140, 500]]
        f107_bins = [[0, 95], [105, 130], [140, 500]]
        db_prefix = "f107_"
        for f107_bin in f107_bins:
            dbName = db_prefix + str(f107_bin[0]) + "_to_" + str(f107_bin[1]) +\
                                        "_" + ftype + ".sqlite"
            dbName_list.append(dbName)


    if binned_imf:
        dbName_list = []
        #imf_clock_angle_bins = [[65, 115], [245, 295]] 
        #imf_clock_angle_bins = [[335, 25], [155, 205]]
        #imf_clock_angle_bins = [[330, 30], [150, 210]]
        #imf_clock_angle_bins = [[60, 120], [240, 300]] 
        #imf_clock_angle_bins = [[315, 45], [135, 225]]
        imf_clock_angle_bins = [[300, 60], [120, 240]]

        bins = imf_clock_angle_bins
        db_prefix = "imf_clock_angle_"
        for bn in bins:
            dbName = db_prefix + str(bn[0]) + "_to_" + str(bn[1]) +\
                                        "_" + ftype + ".sqlite"
            dbName_list.append(dbName)

    jobs = []
    for season in seasons:
        if binned_F107:
            baseLocation="../data/sqlite3/" + season + "/kp_l_3/" + "data_in_mlt" + "/binned_by_f107/"
        if binned_imf:
            baseLocation="../data/sqlite3/" + season + "/kp_l_3/" + "data_in_mlt" + "/binned_by_imf_clock_angle/"
        if binned_season:
            baseLocation="../data/sqlite3/" + season + "/kp_l_3/" + "data_in_mlt" + "/"
        #baseLocation="../data/sqlite3/" + season + "/kp_l_3/" + "data_in_mlt" + "/"
        #baseLocation="../data/sqlite3/" + season + "/kp_l_3/" + "data_in_geo" + "/"
        #baseLocation="../data/sqlite3/" + season + "/kp_l_2/" + "data_in_geo" + "/"
        #baseLocation="../data/sqlite3/" + season + "/kp_l_1/" + "data_in_mlt" + "/"
        #baseLocation="../data/sqlite3/" + season + "/kp_l_1/" + "data_in_geo" + "/"
        for dbName in dbName_list:
            p = multiprocessing.Process(target=worker, args=(baseLocation, dbName))
            jobs.append(p)
            p.start()

        for p in jobs:
            p.join()

