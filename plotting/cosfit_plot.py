def plot_cosfit(latc, lonc, ftype="fitacf", season="winter",
                baseLocation="../data/sqlite3/", sqrt_weighting=True):
    """ plots a the cosfit results for a give latc-lonc grid for a giving season"""

    import sqlite3
    import numpy as np
    import matplotlib.pyplot as plt

    dbName = "master_" + ftype + ".sqlite"

    # make db connection
    conn = sqlite3.connect(baseLocation + dbName)
    cur = conn.cursor()
    
    T1 = "master_summary"
    if sqrt_weighting:
        T2 = "master_cosfit"
    else:
        T2 = "master_cosfit_equal_weighting"

    command = "SELECT median_vel, vel_count, gazmc FROM {tb}\
               WHERE glatc={lat}\
                AND glonc={lon}\
                ORDER BY gazmc"\
                .format(tb=T1, lat=latc, lon=lonc)
    cur.execute(command)
    rows = cur.fetchall()
    median_vel = -np.array([x[0] for x in rows])
    gazmc_vel_count = np.array([x[1] for x in rows])
    weight =  np.sqrt(gazmc_vel_count)
    azm = np.array([x[2] for x in rows])
    azm = [x if x <= 180 else x-360 for x in azm]

    # select the cosine fitting results from db
    command = "SELECT vel_count, vel_mag, vel_mag_err, vel_dir, vel_dir_err FROM {tb}\
               WHERE glatc={lat}\
                AND glonc={lon}"\
               .format(tb=T2, lat=latc, lon=lonc)
    cur.execute(command)
    row = cur.fetchall()[0]
    vel_count = row[0]
    vel_mag = -row[1]
    vel_mag_err = row[2]
    vel_dir = row[3]
    vel_dir_err = row[4]

    # close db connection
    conn.close()

    # plot the LOS data
    fig, ax = plt.subplots()
    ax.scatter(azm, median_vel, marker='o',c='k', s=5*weight, edgecolors="face", label="binned velocities")

    # plot the cosfit curve
    #x_fit = np.arange(0, 360, 1)
    x_fit = np.arange(-180, 180, 0.01)
    y_fit = vel_mag * np.cos(np.deg2rad(x_fit) - np.deg2rad(vel_dir))
    ax.plot(x_fit, y_fit, 'y', linewidth=3, label="fit line")

    # mark the peak position
    ind_max = np.argmax(y_fit)
    y_max = y_fit[ind_max]
    x_max = x_fit[ind_max]
    fsz = 12
    ax.scatter(x_max, y_max, c='r', edgecolors="face", marker = '*', s = 150, label="2D vector", zorder=5)
    ax.annotate('vel=' + '{0:.01f}'.format(y_max) , xy = (0.02, 0.9), xycoords='axes fraction',\
       horizontalalignment='left', verticalalignment='bottom', fontsize=fsz) 
    ax.annotate('azm=' + '{0:.01f}'.format(x_max) +'$^\circ$' , xy = (0.015, 0.85), xycoords='axes fraction',\
       horizontalalignment='left', verticalalignment='bottom', fontsize=fsz) 

    ax.annotate('vel_std=' + '{0:.01f}'.format(vel_mag_err) , xy = (0.02, 0.80), xycoords='axes fraction',\
       horizontalalignment='left', verticalalignment='bottom', fontsize=fsz) 
    ax.annotate('azm_std=' + '{0:.01f}'.format(vel_dir_err) +'$^\circ$' , xy = (0.02, 0.75), xycoords='axes fraction',\
            horizontalalignment='left', verticalalignment='bottom', fontsize=fsz) 
    
    ax.set_xlim([-180, 180])
    ax.set_ylim([-150, 150])

    # put labels
    ax.set_title("Velocity Fitting Results, " + season[0].upper()+season[1:] +\
                 ", MLat = " + str(latc) + ", MLT = " + str(round(lonc/15., 2)))
    ax.set_xlabel("Azimuth [$^\circ$]")
    ax.set_ylabel("Velocity [m/s]")
    ax.legend()

    return fig
    
if __name__ == "__main__":

    import matplotlib.pyplot as plt
    import sys
    sys.path.append("../move_to_db/")
    from bin_data import grids 
    import numpy as np
    import sqlite3

    # points of interest
    latc, lonc = 56.5, 0

    # find points that best matches the user input latc, lonc.
#    grds = grids(lat_min=35, lat_max=90, dlat=1, half_dlat_offset=True)
#    lat_idx = np.digitize(latc, grds.lat_bins) - 1
#    latc = grds.center_lats[lat_idx]
#    lon_idx = np.digitize(lonc, grds.lon_bins[lat_idx]) - 1
#    lonc = grds.center_lons[lat_idx][lon_idx]


    #sqrt_weighting=False
    sqrt_weighting=True
    ftype = "fitacf"
    season = "winter"
    #season = "summer"
    #season = "equinox"
    #fig_dir = "./plots/cosfit_plot/kp_l_3/data_in_geo/" + season + "/"
    fig_dir = "./plots/cosfit_plot/kp_l_3/data_in_mlt/" + season + "/"
    #fig_name = "seasonal_cosfit_mlat"+str(latc) + "_mlt" + str(round(lonc/15., 2))
    fig_name = season + "_cosfit_geolat"+str(latc) + "_ltm" + str(round(lonc/15., 2))

    # select points of interest
    dbName = "master_" + ftype + ".sqlite"
    #baseLocation="../data/sqlite3/" + season + '/'
   #baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_geo/'
    baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_mlt/'

    # make db connection
    conn = sqlite3.connect(baseLocation + dbName)
    cur = conn.cursor()
    if sqrt_weighting:
        T2 = "master_cosfit"
    else:
        T2 = "master_cosfit_equal_weighting"

    command = "SELECT glatc, glonc FROM {tb}".format(tb=T2)
    cur.execute(command)
    rows = cur.fetchall()
    all_lats = np.array([x[0] for x in rows])
    matching_lats_idx = np.where(all_lats==latc)
    all_lons = np.array([x[1] for x in rows])
    possible_lons = all_lons[matching_lats_idx]
    lonc_idx = (np.abs(possible_lons - lonc)).argmin()
    lonc = round(possible_lons[lonc_idx],2)
    
#    lat_idx = np.digitize(latc, grds.lat_bins) - 1
    
    # close db connection
    conn.close()

    # plotting

    fig = plot_cosfit(latc, lonc, ftype=ftype, season=season, sqrt_weighting=sqrt_weighting,
                      baseLocation=baseLocation)
    #fig.savefig(fig_dir + fig_name + ".png", dpi=300)
    plt.show()
