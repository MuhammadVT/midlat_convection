def fetch_data(ftype="fitacf", nvel_min=250, lat_range=[52, 59],
               baseLocation="./", dbName=None):

    """ fetch the data from the master db into a dict
    nvel_min : int
        minimum number of velocity measurements in a lat-lon grid cell

    Return
    ------
    data_dict : dict
    """
    import sqlite3
    import datetime as dt
    import numpy as np 

    # make db connection for dopsearch
    if dbName is None:
        dbName = "master_" + ftype + ".sqlite"
    conn = sqlite3.connect(baseLocation + dbName)
    cur = conn.cursor()
    
    T1 = "master_cosfit"

    # select data from the master_cosfit table for the night side
    command = "SELECT vel_count, vel_mag, vel_dir, glatc, glonc, vel_mag_err, vel_dir_err FROM {tb1}\
               WHERE glatc BETWEEN {lat_min} AND {lat_max}\
               ".format(tb1=T1, lat_min=lat_range[0], lat_max=lat_range[1])
    cur.execute(command)
    rws = cur.fetchall()

    data_dict = {}
    # filter the data based on lattitude range. 
    data_dict['vel_count'] = np.array([x[0] for x in rws if x[0] >= nvel_min])
    data_dict['vel_mag'] = np.array([x[1] for x in rws if x[0] >= nvel_min])
    data_dict['vel_dir'] = np.array([x[2] for x in rws if x[0] >= nvel_min])
    data_dict['glatc'] = np.array([x[3] for x in rws if x[0] >= nvel_min])
    data_dict['glonc'] = np.array([x[4] for x in rws if x[0] >= nvel_min])
    data_dict['vel_mag_err'] = np.array([x[5] for x in rws if x[0] >= nvel_min])
    data_dict['vel_dir_err'] = np.array([x[6] for x in rws if x[0] >= nvel_min])


    # close db connection
    conn.close()

    return data_dict

def cart2pol(x, y):
    import numpy as np
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (phi, rho)

def pol2cart(phi, rho):
    import numpy as np
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)


def vector_plot(ax, data_dict, cmap, bounds, velscl=1, lat_min=50, title="xxx",
                for_HWM=False, hemi="north", fake_pole=False):
    
    """ plots the flow vectors in MLT coords """

    import matplotlib.pyplot as plt
    import matplotlib as mpl

    from matplotlib.collections import PolyCollection,LineCollection
    import numpy as np

    # build a custom color map
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # set axis limits
    if fake_pole:
        fake_pole_lat = 70
    else:
        fake_pole_lat = 90

    rmax = fake_pole_lat - lat_min
    ax.set_xlim([-rmax, rmax])
    ax.set_ylim([-rmax, 0])
    ax.set_aspect("equal")

    # remove tikcs
    ax.tick_params(axis='both', which='both', bottom='off', top='off',
                   left="off", right="off", labelbottom='off', labelleft='off')

    # plot the latitudinal circles
    for r in [10, 30, 50]:
        c = plt.Circle((0, 0), radius=r, fill=False)
        ax.add_patch(c)

    # plot the longitudinal lines
    for l in np.deg2rad(np.array([210, 240, 270, 300, 330])):
        x1, y1 = pol2cart(l, 10) 
        x2, y2 = pol2cart(l, 50) 
        ax.plot([x1, x2], [y1, y2], 'k')

    x1, y1 = pol2cart(np.deg2rad(data_dict['glonc']-90), fake_pole_lat-np.abs(data_dict['glatc']))

    # add the vector lines
    lines = []
    intensities = []
    vel_mag = data_dict['vel_mag']
    # calculate the angle of the vectors in a tipical x-y axis.
    theta = np.deg2rad(data_dict['glonc'] + 90 - data_dict['vel_dir']) 

    if for_HWM:
        x2 = x1+vel_mag/velscl*(1.0)*np.cos(theta)
        y2 = y1+vel_mag/velscl*(1.0)*np.sin(theta)
    else:
        x2 = x1+vel_mag/velscl*(-1.0)*np.cos(theta)
        y2 = y1+vel_mag/velscl*(-1.0)*np.sin(theta)
    
    lines.extend(zip(zip(x1,y1),zip(x2,y2)))
    #save the param to use as a color scale
    intensities.extend(np.abs(vel_mag))

    # plot the velocity locations
    ccoll = ax.scatter(x1, y1,
                        s=1.0,zorder=10,marker='o', c=np.abs(np.array(intensities)),
                        linewidths=.5, edgecolors='face'
                        ,cmap=cmap,norm=norm)
    lcoll = LineCollection(np.array(lines),linewidths=0.5,zorder=12
                        ,cmap=cmap,norm=norm)
    lcoll.set_array(np.abs(np.array(intensities)))
    ccoll.set_array(np.abs(np.array(intensities)))
    ax.add_collection(ccoll)
    ax.add_collection(lcoll)

    # add text
    ax.set_title(title, fontsize='small')
    # add latitudinal labels
    fnts = 'x-small'
    if hemi=="north":
        ax.annotate("80", xy=(0, -10), ha="left", va="bottom", fontsize=fnts)
        ax.annotate("60", xy=(0, -30), ha="left", va="bottom", fontsize=fnts)
    elif hemi=="south":
        ax.annotate("-80", xy=(0, -10), ha="left", va="bottom", fontsize=fnts)
        ax.annotate("-60", xy=(0, -30), ha="left", va="bottom", fontsize=fnts)
    # add mlt labels
    ax.annotate("0", xy=(0, -rmax), ha="center", va="top", fontsize=fnts)
    ax.annotate("6", xy=(rmax, 0), ha="left", va="center", fontsize=fnts)
    ax.annotate("18", xy=(-rmax, 0), ha="right", va="center", fontsize=fnts)

    return lcoll

def add_cbar(fig, coll, bounds, label="Velocity [m/s]", cax=None):

    # add color bar
    if cax:
        cbar=fig.colorbar(coll, cax=cax, orientation="vertical",
                          boundaries=bounds, drawedges=False) 
    else:
        cbar=fig.colorbar(coll, orientation="vertical", shrink=.65,
                          boundaries=bounds, drawedges=False) 


    #define the colorbar labels
    l = []
    for i in range(0,len(bounds)):
        if i == 0 or i == len(bounds)-1:
            l.append(' ')
            continue
        l.append(str(int(bounds[i])))
    cbar.ax.set_yticklabels(l)
    #cbar.ax.tick_params(axis='y',direction='out')
    cbar.set_label(label)

    return


def main(by_season=True, by_f107=False, by_imf_clock_angle=False):

    import datetime as dt
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    # input parameters
    #nvel_min=250
    #nvel_min=300
    nvel_min=50
    lat_range=[52, 59]
    #lat_range=[50, 90]
    #lat_range=[50, 70]
    lat_min = 50
    #lat_range=[39, 90]
    #lat_range=[40, 60]
    #lat_min = 39
    ftype = "fitacf"
    #ftype = "fitex"

    # cmap and bounds for color bar
    color_list = ['purple', 'b', 'c', 'g', 'y', 'r']
    cmap = mpl.colors.ListedColormap(color_list)
    bounds = [0., 8, 17, 25, 33, 42, 10000]

    seasons = ["winter", "summer", "equinox"]
    #seasons = ["winter"]

    if by_season:
        #fig_dir = "./plots/convection/kp_l_3/data_in_geo/"
        #fig_dir = "./plots/convection/kp_l_2/data_in_mlt/seasonal/"
        #fig_dir = "./plots/convection/kp_l_3/data_in_mlt/seasonal/"
        #fig_dir = "./plots/convection/kp_l_1/data_in_mlt/seasonal/"
        #fig_dir = "./plots/convection/kp_l_2/data_in_geo/seasonal/"
        fig_dir = "/home/muhammad/Dropbox/mypapers/paper_02/version_01/figures/"
        fig_name = "seasonal_convection_lat" + str(lat_range[0]) +"_to_lat" + str(lat_range[1])
       
        # create subplots
        fig, axes = plt.subplots(nrows=len(seasons), ncols=1, figsize=(6,8))
        fig.subplots_adjust(hspace=0.3)


        if len(seasons) == 1:
            axes = [axes]

        for i, season in enumerate(seasons):
            # fetches the data from db 
            #baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_geo/'
            #baseLocation="../data/sqlite3/" + season + '/kp_l_2/data_in_mlt/'
            baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_mlt/'
            #baseLocation="../data/sqlite3/" + season + '/kp_l_1/data_in_mlt/'
            #baseLocation="../data/sqlite3/" + season + '/kp_l_2/data_in_geo/'
            data_dict = fetch_data(ftype=ftype, nvel_min=nvel_min, 
                            lat_range=lat_range, baseLocation=baseLocation,
                            dbName=None)

            # plot the flow vectors
            title = "Velocities, " + season[0].upper()+season[1:] + ", Kp<3"
            coll = vector_plot(axes[i], data_dict, cmap, bounds, velscl=10,
                               lat_min=lat_min, title=title)

        # add colorbar
        fig.subplots_adjust(right=0.80)
        cbar_ax = fig.add_axes([0.85, 0.25, 0.02, 0.5])
        add_cbar(fig, coll, bounds, cax=cbar_ax, label="Velocity [m/s]")
        # save the fig
        fig.savefig(fig_dir + fig_name + ".png", dpi=500)
        #plt.show()

    if by_f107:

        #f107_bins = [[0, 100], [100, 175], [175, 500]]
        #f107_bins = [[0, 100], [100, 150], [150, 500]]
        #f107_bins = [[0, 100], [105, 130], [130, 500]]
        #f107_bins = [[0, 90], [105, 125], [130, 500]]
        f107_bins = [[0, 95], [105, 130], [140, 500]]
        #f107_bins = [[0, 105], [105, 125], [125, 500]]
        #f107_bins = [[0, 110], [110, 500]]
        #f107_bins = [[0, 120], [120, 500]]

        for season in seasons:

            # create subplots
            fig, axes = plt.subplots(nrows=len(f107_bins), ncols=1, figsize=(6,8))
            fig.subplots_adjust(hspace=0.3)
            if len(f107_bins) == 1:
                axes = [axes]

            fig_dir = "./plots/convection/kp_l_3/data_in_mlt/binned_by_f107/"
            fig_name = season + "_convection_by_f107_lat" + str(lat_range[0]) +"_to_lat" + str(lat_range[1])

            # fetches the data from db 
            #baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_mlt/'
            baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_mlt/binned_by_f107/'
            
            # loop through the f107 bins
            for i, f107_bin in enumerate(f107_bins):
                dbName = "f107_" + str(f107_bin[0]) + "_to_" + str(f107_bin[1]) +\
                            "_" + ftype + ".sqlite"

                data_dict = fetch_data(ftype=ftype, nvel_min=nvel_min, 
                                lat_range=lat_range, baseLocation=baseLocation,
                                dbName=dbName)

                # plot the flow vectors
                title = "Velocities, " + season[0].upper()+season[1:] + ", Kp<3" + ", F10.7=" + str(f107_bin)
                coll = vector_plot(axes[i], data_dict, cmap, bounds, velscl=10,
                                   lat_min=lat_min, title=title)

            # add colorbar
            fig.subplots_adjust(right=0.80)
            cbar_ax = fig.add_axes([0.85, 0.25, 0.02, 0.5])
            add_cbar(fig, coll, bounds, cax=cbar_ax, label="Velocity [m/s]")
            # save the fig
            fig.savefig(fig_dir + fig_name + ".png", dpi=500)
            #plt.show()

    if by_imf_clock_angle:
        
        # imf clock angle bins
        #bins = [[65, 115], [245, 295]] 
        #bins = [[335, 25], [155, 205]]

        #bins = [[330, 30], [150, 210]]
        #bins = [[315, 45], [135, 225]]
        bins = [[300, 60], [120, 240]]
        bins_txt = ["Bz+", "Bz-"]
        fig_txt =  "_span120_convection_by_imf_Bz_clock_angle_lat_"

#        bins = [[60, 120], [240, 300]] 
#        bins_txt = ["By+", "By-"]
#        fig_txt =  "_convection_by_imf_By_clock_angle_lat"

        for season in seasons:

            # create subplots
            fig, axes = plt.subplots(nrows=len(bins), ncols=1, figsize=(6,8))
            fig.subplots_adjust(hspace=0.3)
            if len(bins) == 1:
                axes = [axes]

            fig_dir = "./plots/convection/kp_l_3/data_in_mlt/binned_by_imf_clock_angle/"
            fig_name = season + fig_txt + str(lat_range[0]) +"_to_lat" + str(lat_range[1])

            # fetches the data from db 
            #baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_mlt/'
            baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_mlt/binned_by_imf_clock_angle/'
            
            # loop through the bins
            for i, bn in enumerate(bins):
                dbName = "imf_clock_angle_" + str(bn[0]) + "_to_" + str(bn[1]) +\
                            "_" + ftype + ".sqlite"

                data_dict = fetch_data(ftype=ftype, nvel_min=nvel_min, 
                                lat_range=lat_range, baseLocation=baseLocation,
                                dbName=dbName)

                # plot the flow vectors
                title = "Velocities, " + season[0].upper()+season[1:] + ", Kp<3" + ", IMF " + bins_txt[i]
                coll = vector_plot(axes[i], data_dict, cmap, bounds, velscl=10,
                                   lat_min=lat_min, title=title)

            # add colorbar
            fig.subplots_adjust(right=0.80)
            cbar_ax = fig.add_axes([0.85, 0.25, 0.02, 0.5])
            add_cbar(fig, coll, bounds, cax=cbar_ax, label="Velocity [m/s]")
            # save the fig
            fig.savefig(fig_dir + fig_name + ".png", dpi=500)
            #plt.show()


    return

if __name__ == "__main__":
    by_season=False
    by_f107=False
    by_imf_clock_angle=True
    main(by_season=by_season, by_f107=by_f107, by_imf_clock_angle=by_imf_clock_angle)
