def vel_vs_lt(ax, data_dict, veldir="zonal", center_at_zero_mlt=True,
               glatc_list=None, title="xxx", add_err_bar=False):
    
    """ plots the flow vectors in local time (MLT or SLT) coords

    parameters
    ----------
    veldir : str
        veocity component. if set to "all" then it means the velocity magnitude
    """

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.collections import PolyCollection,LineCollection
    import numpy as np

    # calculate velocity components
    vel_mag = data_dict['vel_mag']
    vel_dir = np.deg2rad(data_dict['vel_dir'])
    vel_mag_err = data_dict['vel_mag_err']

    if veldir == "zonal":
        vel_comp = vel_mag*(-1.0)*np.sin(vel_dir)
        vel_comp_err = vel_mag_err*(-1.0)*np.sin(vel_dir)
    elif veldir == "meridional":
        vel_comp = vel_mag*(-1.0)*np.cos(vel_dir)
        vel_comp_err = vel_mag_err*(-1.0)*np.cos(vel_dir)
    elif veldir == "all":
        vel_comp = np.abs(vel_mag)
        vel_comp_err = vel_mag_err
    vel_mlt = data_dict['glonc'] / 15.
    
    # colors of different lines
    color_list = ['darkblue', 'b', 'dodgerblue', 'c', 'g', 'orange', 'r']
    color_list.reverse()

    # MLATs
    if glatc_list is None:
        glatc_list = np.array([50.5])
    for jj, mlat in enumerate(glatc_list):
        vel_comp_jj = [vel_comp[i] for i in range(len(vel_comp)) if data_dict['glatc'][i] == mlat]
        vel_mlt_jj = [vel_mlt[i] for i in range(len(vel_comp)) if data_dict['glatc'][i] == mlat]
        vel_comp_err_jj = [vel_comp_err[i] for i in range(len(vel_comp_err)) if data_dict['glatc'][i] == mlat]
        if center_at_zero_mlt:
            # center at 0 MLT
            vel_mlt_jj = [x if x <=12 else x-24 for x in vel_mlt_jj]
        # plot the velocities for each MLAT
        ax.scatter(vel_mlt_jj, vel_comp_jj, c=color_list[jj],
                #marker='o', s=3, linewidths=.5, edgecolors='face', label=str(int(mlat)))
                marker='o', s=3, linewidths=.5, edgecolors='face', label=str(mlat))

        if add_err_bar:
            ax.errorbar(vel_mlt_jj, vel_comp_jj, yerr=vel_comp_err_jj, mfc=color_list[jj],
                    #marker='o', s=3, linewidths=.5, edgecolors='face', label=str(int(mlat)))
                    fmt='o', ms=2, elinewidth=.5, mec=color_list[jj], ecolor="k")

    # add text
    ax.set_title(title, fontsize="small")

    # add zero-line
    if veldir != "all":
        ax.axhline(y=0, color='k', linewidth=0.7)

    # set axis limits
    if center_at_zero_mlt:
        ax.set_xlim([-12, 12])
        # add legend
        ax.legend(bbox_to_anchor=(0.12, 0.96), fontsize=8)
    else:
        ax.set_xlim([0, 24])
        # add legend
        #ax.legend(loc='center right', fontsize=8)
        ax.legend(bbox_to_anchor=(0.65, 0.96), fontsize=8)

    if veldir == "all":
        ax.set_ylim([0, 60])
    else:
        ax.set_ylim([-60, 50])
    
    # axis labels
    ax.set_ylabel("Vel [m/s]")

    return

def main(by_season=True, by_f107=False, by_imf_clock_angle=False):

    import datetime as dt
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from convection import fetch_data
    from matplotlib.ticker import MultipleLocator
    import numpy as np

    # input parameters
    nvel_min=300
    #nvel_min=250
    #nvel_min=50
    #del_lat=3
    del_lat=1
    #lat_range=[58, 65]
    lat_range=[52, 59]

    #lat_range=[40, 54]
    #lat_range=[39, 57]
    #lat_range=[48, 55]
    #lat_range=[42, 49]
    #lat_range=[41, 48]
    glatc_list = np.arange(lat_range[1]-0.5, lat_range[0]-0.5, -del_lat)
    #lat_range=[40, 60]
    #glatc_list = np.array([41.5, 50.5])
    #print glatc_list
    if len(glatc_list) == 0:
        glatc_list = np.array([lat_range[0]]+0.5)
    
    #add_err_bar = True
    add_err_bar = False
    ftype = "fitacf"
    #ftype = "fitex"

    #veldir="all"
    #veldir="zonal"
    veldir="meridional"
    #center_at_zero_mlt=True
    center_at_zero_mlt=False

    seasons = ["winter", "summer", "equinox"]
    #seasons = ["winter"]

    if by_season:

        #fig_dir = "./plots/velcomp/kp_l_3/data_in_geo/seasonal/"
        #fig_dir = "./plots/velcomp/kp_l_2/data_in_geo/seasonal/"
        #fig_dir = "./plots/velcomp/kp_l_1/data_in_geo/seasonal/"
        #fig_dir = "./plots/velcomp/kp_l_1/data_in_mlt/seasonal/"
        #fig_dir = "./plots/velcomp/kp_l_2/data_in_mlt/seasonal/"
        #fig_dir = "./plots/velcomp/kp_l_3/data_in_mlt/seasonal/"
        fig_dir = "/home/muhammad/Dropbox/mypapers/paper_02/version_01/figures/"
        if center_at_zero_mlt:
            fig_name = "seasonal_" + veldir+ "_vel_vs_ltm_c0" +\
                       "_lat" + str(lat_range[0]) + "_to_lat" + str(lat_range[1])
        else:
            fig_name = "seasonal_" + veldir+ "_vel_vs_ltm" +\
                       "_lat" + str(lat_range[0]) + "_to_lat" + str(lat_range[1])

        # create subplots
        fig, axes = plt.subplots(nrows=len(seasons), ncols=1, figsize=None, sharex=True)
        fig.subplots_adjust(hspace=0.3)

        if len(seasons) == 1:
            axes = [axes]

        for i, season in enumerate(seasons):
            # fetches the data from db 
            #baseLocation="../data/sqlite3/" + season + '/'
            baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_mlt/'
            #baseLocation="../data/sqlite3/" + season + '/kp_l_2/data_in_mlt/'
            #baseLocation="../data/sqlite3/" + season + '/kp_l_1/data_in_mlt/'
            #baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_geo/'
            #baseLocation="../data/sqlite3/" + season + '/kp_l_2/data_in_geo/'
            #baseLocation="../data/sqlite3/" + season + '/kp_l_1/data_in_geo/'
            data_dict = fetch_data(ftype=ftype, nvel_min=nvel_min, 
                            lat_range=lat_range, baseLocation=baseLocation,
                            dbName=None)

            # plot the flow vector components
            if veldir == "all" :
                title = "Velocity Magnitude, " + season[0].upper()+season[1:] + ", Kp<3"
            else:
                title = veldir[0].upper()+veldir[1:] + " Velocities, " +\
                        season[0].upper()+season[1:] + ", Kp<3"
            vel_vs_lt(axes[i], data_dict, veldir=veldir, center_at_zero_mlt=center_at_zero_mlt,
                    glatc_list=glatc_list, title=title, add_err_bar=add_err_bar)

        # set axis label
        axes[-1].set_xlabel("MLT")
        #axes[-1].set_xlabel("Solar Local Time")
        axes[-1].xaxis.set_major_locator(MultipleLocator(base=3))
        if center_at_zero_mlt:
            xlabels = [item.get_text() for item in axes[-1].get_xticklabels()]
            xlabels = [str(x) for x in range(12, 24, 3) + range(0, 15, 3)]
            plt.xticks(range(-12, 15, 3), xlabels)

        # save the fig
        fig.savefig(fig_dir + fig_name + ".png", dpi=300)
        #plt.show()

    if by_f107:
        #f107_bins = [[0, 100], [100, 175], [175, 500]]
        #f107_bins = [[0, 100], [100, 150], [150, 500]]
        #f107_bins = [[0, 100], [105, 130], [130, 500]]
        #f107_bins = [[0, 105], [105, 125], [125, 500]]
        #f107_bins = [[0, 95], [105, 125], [130, 500]]
        f107_bins = [[0, 95], [105, 130], [140, 500]]
        #f107_bins = [[0, 110], [110, 500]]
        #f107_bins = [[0, 120], [120, 500]]

        for season in seasons:

            # create subplots
            fig, axes = plt.subplots(nrows=len(f107_bins), ncols=1, figsize=None, sharex=True)
            fig.subplots_adjust(hspace=0.3)
            if len(f107_bins) == 1:
                axes = [axes]

            #fig_dir = "./plots/velcomp/kp_l_3/data_in_geo/"
            fig_dir = "./plots/velcomp/kp_l_3/data_in_mlt/binned_by_f107/"
            if center_at_zero_mlt:
                fig_name = season + "_" + veldir+ "_vel_vs_ltm_c0" +\
                           "_lat" + str(lat_range[0]) + "_to_lat" + str(lat_range[1])
            else:
                fig_name = season + "_" + veldir+ "_vel_vs_ltm" +\
                           "_lat" + str(lat_range[0]) + "_to_lat" + str(lat_range[1])

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

                # plot the flow vector components
                if veldir == "all" :
                    title = "Velocity Magnitude, " + season[0].upper()+season[1:] + ", Kp<2"
                else:
                    title = veldir[0].upper()+veldir[1:] + " Velocities, " + season[0].upper()+season[1:] + ", Kp<3" + ", F10.7=" + str(f107_bin)
                vel_vs_lt(axes[i], data_dict, veldir=veldir, center_at_zero_mlt=center_at_zero_mlt,
                        glatc_list=glatc_list, title=title, add_err_bar=add_err_bar)

            # set axis label
            axes[-1].set_xlabel("MLT")
            axes[-1].xaxis.set_major_locator(MultipleLocator(base=3))
            if center_at_zero_mlt:
                xlabels = [item.get_text() for item in axes[-1].get_xticklabels()]
                xlabels = [str(x) for x in range(12, 24, 3) + range(0, 15, 3)]
                plt.xticks(range(-12, 15, 3), xlabels)

            # save the fig
            fig.savefig(fig_dir + fig_name + ".png", dpi=300)
            #plt.show()

    if by_imf_clock_angle:
        
        # imf clock angle bins
        #bins = [[330, 30], [150, 210]]
        #bins_txt = ["Bz+", "Bz-"]
        #fig_txt =  "_by_imf_Bz_clock_angle"
        bins = [[60, 120], [240, 300]] 
        bins_txt = ["By+", "By-"]
        fig_txt =  "_convection_by_imf_By_clock_angle_lat"


        for season in seasons:

            # create subplots
            fig, axes = plt.subplots(nrows=len(bins), ncols=1, figsize=None, sharex=True)
            fig.subplots_adjust(hspace=0.3)
            if len(bins) == 1:
                axes = [axes]

            #fig_dir = "./plots/velcomp/kp_l_3/data_in_geo/"
            fig_dir = "./plots/velcomp/kp_l_3/data_in_mlt/binned_by_imf_clock_angle/"
            if center_at_zero_mlt:
                fig_name = season + "_" + veldir+ "_vel_vs_ltm_c0" +  fig_txt +\
                           "_lat" + str(lat_range[0]) + "_to_lat" + str(lat_range[1])
            else:
                fig_name = season + "_" + veldir+ "_vel_vs_ltm" + fig_txt +\
                           "_lat" + str(lat_range[0]) + "_to_lat" + str(lat_range[1])

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

                # plot the flow vector components
                if veldir == "all" :
                    title = "Velocity Magnitude, " + season[0].upper()+season[1:] + ", Kp<2"
                else:
                    title = veldir[0].upper()+veldir[1:] + " Velocities, " + season[0].upper()+season[1:] + \
                        ", Kp<3" + ", IMF " + bins_txt[i]
                vel_vs_lt(axes[i], data_dict, veldir=veldir, center_at_zero_mlt=center_at_zero_mlt,
                        glatc_list=glatc_list, title=title, add_err_bar=add_err_bar)

            # set axis label
            axes[-1].set_xlabel("MLT")
            axes[-1].xaxis.set_major_locator(MultipleLocator(base=3))
            if center_at_zero_mlt:
                xlabels = [item.get_text() for item in axes[-1].get_xticklabels()]
                xlabels = [str(x) for x in range(12, 24, 3) + range(0, 15, 3)]
                plt.xticks(range(-12, 15, 3), xlabels)

            # save the fig
            fig.savefig(fig_dir + fig_name + ".png", dpi=300)
            #plt.show()


    return

if __name__ == "__main__":
    by_season=True
    by_f107=False
    #by_season=False
    #by_f107=True
    by_imf_clock_angle=False
    #by_imf_clock_angle=True
    main(by_season=by_season, by_f107=by_f107, by_imf_clock_angle=by_imf_clock_angle)
