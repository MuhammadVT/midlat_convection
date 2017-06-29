def cosfit_error(ax, data_dict, cmap, bounds, err_type='Magnitude', lat_min=50, title="xxx"):
    
    """ plots the flow vectors in MLT coords """

    import matplotlib.pyplot as plt
    import matplotlib as mpl

    from matplotlib.collections import PolyCollection,LineCollection
    import numpy as np

    from convection import pol2cart 

    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # plot the backgroud MLT coords
    rmax = 90 - lat_min
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

    # calculate the velocity locations in cartisian coords
    x1, y1 = pol2cart(np.deg2rad(data_dict['glonc']-90), 90-data_dict['glatc'])

    #save the param to use as a color scale
    vel_mag_err_ratio = np.abs(np.array(data_dict['vel_mag_err']) / np.array(data_dict['vel_mag']))
    vel_dir_err_ratio = np.abs(np.array(data_dict['vel_dir_err']) / np.array(data_dict['vel_dir']))
    vel_err_ratio = vel_mag_err_ratio * vel_dir_err_ratio
    # normalize
    #vel_err_ratio = vel_err_ratio / np.max(vel_err_ratio)
    if err_type=="Magnitude":
        intensities = vel_mag_err_ratio.tolist()
    if err_type=="Phase":
        intensities = vel_dir_err_ratio.tolist()
    if err_type=="both":
        intensities = vel_err_ratio.tolist()

    #do the actual overlay
    ccoll = ax.scatter(x1, y1,
                    s=7.0,zorder=10,marker='o', c=np.abs(np.array(intensities)),
                    linewidths=.5, edgecolors='face'
                    ,cmap=cmap,norm=norm)

    # add labels
    ax.set_title(title, fontsize='small')
    # add latitudinal labels
    fnts = "x-small"
    ax.annotate("80", xy=(0, -10), ha="left", va="bottom", fontsize=fnts)
    ax.annotate("60", xy=(0, -30), ha="left", va="bottom", fontsize=fnts)
    # add mlt labels
    ax.annotate("0", xy=(0, -rmax), ha="center", va="top", fontsize=fnts)
    ax.annotate("6", xy=(rmax, 0), ha="left", va="center", fontsize=fnts)
    ax.annotate("18", xy=(-rmax, 0), ha="right", va="center", fontsize=fnts)

    return  ccoll

def add_cbar(fig, coll, bounds, label="Fitting Error Ratio", cax=None):

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
        if i == len(bounds)-1:
            l.append(' ')
            continue
        l.append(str(bounds[i]))
    cbar.ax.set_yticklabels(l)
    #cbar.ax.tick_params(axis='y',direction='out')
    cbar.set_label(label)

    return

def main(by_season=True, by_f107=False, by_imf_clock_angle=False):
    import datetime as dt
    import matplotlib.pyplot as plt
    from convection import fetch_data
    import matplotlib as mpl
    import numpy as np

    # input parameters
    #nvel_min=50
    nvel_min=300
    #lat_range=[52, 59]
    lat_range=[52, 65]
    lat_min = 50
    #lat_range=[40, 60]
    #lat_min = 43
    ftype = "fitacf"
    #ftype = "fitex"


    # choose a fitting error type
    err_type = "Magnitude"
    #err_type = "Phase"
    #err_type = "both"

    # build a custom color map and bounds
    color_list = ['purple', 'b', 'dodgerblue', 'c', 'g', 'y', 'orange', 'r']
    cmap = mpl.colors.ListedColormap(color_list)
    if err_type == "both":
        bounds = np.arange(0, 0.09, 0.01).tolist()
    else:
        bounds = np.arange(0, 0.9, 0.1).tolist()

    seasons = ["winter", "summer", "equinox"]
    #seasons = ["winter", "summer"]
   
    if by_season:
        #fig_dir = "./plots/cosfit_plot/kp_l_3/data_in_geo/"
        #fig_dir = "./plots/cosfit_plot/kp_l_2/data_in_mlt/seasonal/"
        fig_dir = "./plots/cosfit_plot/kp_l_3/data_in_mlt/seasonal/"
        #fig_dir = "/home/muhammad/Dropbox/mypapers/paper_02/version_01/figures/"
        fig_name = "seasonal_cosfit_quality_" + err_type + "_" + str(lat_range[0]) +"_to_lat" + str(lat_range[1])

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
            data_dict = fetch_data(ftype=ftype, nvel_min=nvel_min, 
                            lat_range=lat_range, baseLocation=baseLocation,
                            dbName=None)

            # plot the flow vectors
            title = "Fitting Quality, " + season[0].upper()+season[1:] + ", Kp<3"
            coll = cosfit_error(axes[i], data_dict, cmap, bounds, err_type=err_type,
                                lat_min=lat_min, title=title)

        # add colorbar
        fig.subplots_adjust(right=0.80)
        cbar_ax = fig.add_axes([0.85, 0.25, 0.02, 0.5])
        if err_type == "both":
            cbar_label = "Fitted Vel Magnitude & Phase Error Ratio"
        else:
            #cbar_label = "Fitted Vel " + err_type + " Error Ratio"
            cbar_label = "Fitted Velocity Error Ratio"
        add_cbar(fig, coll, bounds, cax=cbar_ax, label=cbar_label)

        # save the fig
        fig.savefig(fig_dir + fig_name + ".png", dpi=300)
        #plt.show()

    if by_f107:

        #f107_bins = [[0, 100], [100, 175], [175, 500]]
        #f107_bins = [[0, 100], [100, 150], [150, 500]]
        f107_bins = [[0, 105], [105, 125], [125, 500]]

        for season in seasons:

            # create subplots
            fig, axes = plt.subplots(nrows=len(f107_bins), ncols=1, figsize=(6,8))
            fig.subplots_adjust(hspace=0.3)
            if len(f107_bins) == 1:
                axes = [axes]

            fig_dir = "./plots/cosfit_plot/kp_l_3/data_in_mlt/binned_by_f107/"
            fig_name = season + "_binnedby_f107_cosfit_quality_" + err_type + "_" + str(lat_range[0]) +"_to_lat" + str(lat_range[1])

            # fetches the data from db 
            baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_mlt/binned_by_f107/'
            
            # loop through the f107 bins
            for i, f107_bin in enumerate(f107_bins):
                dbName = "f107_" + str(f107_bin[0]) + "_to_" + str(f107_bin[1]) +\
                            "_" + ftype + ".sqlite"

                data_dict = fetch_data(ftype=ftype, nvel_min=nvel_min, 
                                lat_range=lat_range, baseLocation=baseLocation,
                                dbName=dbName)

                # plot the flow vectors
                title = "Fitting Quality, " + season[0].upper()+season[1:] + ", Kp<3" + ", F10.7=" + str(f107_bin)
                coll = cosfit_error(axes[i], data_dict, cmap, bounds, err_type=err_type,
                                    lat_min=lat_min, title=title)
            # add colorbar
            fig.subplots_adjust(right=0.80)
            cbar_ax = fig.add_axes([0.85, 0.25, 0.02, 0.5])
            if err_type == "both":
                cbar_label = "Fitted Vel Magnitude & Phase Error Ratio"
            else:
                cbar_label = "Fitted Vel " + err_type + " Error Ratio"
            add_cbar(fig, coll, bounds, cax=cbar_ax, label=cbar_label)

            # save the fig
            fig.savefig(fig_dir + fig_name + ".png", dpi=300)
            #plt.show()

    if by_imf_clock_angle:
        
        # imf clock angle bins
        bins = [[65, 115], [245, 295]] 

        for season in seasons:

            # create subplots
            fig, axes = plt.subplots(nrows=len(bins), ncols=1, figsize=(6,8))
            fig.subplots_adjust(hspace=0.3)
            if len(bins) == 1:
                axes = [axes]

            fig_dir = "./plots/cosfit_plot/kp_l_3/data_in_mlt/binned_by_imf_clock_angle/"
            fig_name = season + "_binned_by_imf_By_clock_angle_cosfit_quality_" + err_type + "_" + str(lat_range[0]) +"_to_lat" + str(lat_range[1])

            # fetches the data from db 
            baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_mlt/binned_by_imf_clock_angle/'
            
            # loop through the bins
            for i, bn in enumerate(bins):
                dbName = "imf_clock_angle_" + str(bn[0]) + "_to_" + str(bn[1]) +\
                            "_" + ftype + ".sqlite"

                data_dict = fetch_data(ftype=ftype, nvel_min=nvel_min, 
                                lat_range=lat_range, baseLocation=baseLocation,
                                dbName=dbName)

                # plot the flow vectors
                title = "Fitting Quality, " + season[0].upper()+season[1:] + ", Kp<3" + \
                        r", $\theta$=" + str(bn)
                coll = cosfit_error(axes[i], data_dict, cmap, bounds, err_type=err_type,
                                    lat_min=lat_min, title=title)
            # add colorbar
            fig.subplots_adjust(right=0.80)
            cbar_ax = fig.add_axes([0.85, 0.25, 0.02, 0.5])
            if err_type == "both":
                cbar_label = "Fitted Vel Magnitude & Phase Error Ratio"
            else:
                cbar_label = "Fitted Vel " + err_type + " Error Ratio"
            add_cbar(fig, coll, bounds, cax=cbar_ax, label=cbar_label)

            # save the fig
            fig.savefig(fig_dir + fig_name + ".png", dpi=300)
            #plt.show()



    return

if __name__ == "__main__":
    by_season=True
    by_f107=False
    by_imf_clock_angle=False

    main(by_season=by_season, by_f107=by_f107, by_imf_clock_angle=by_imf_clock_angle)
