def num_plot(ax, data_dict, season, cmap, bounds,
             lat_min=50, title="xxx"):
    
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

    # plot the velocity locations
    x1, y1 = pol2cart(np.deg2rad(data_dict['glonc']-90), 90-data_dict['glatc'])

    intensities = []
    vel_count = data_dict['vel_count']
    
    #save the param to use as a color scale
    intensities.extend(np.abs(vel_count))

    #do the actual overlay
    ccoll = ax.scatter(x1, y1,
                    s=7.0,zorder=10,marker='o', c=np.abs(np.array(intensities)),
                    linewidths=.5, edgecolors='face'
                    ,cmap=cmap,norm=norm)

    # add labels
    ax.set_title(title, fontsize="small")
    # add latitudinal labels
    fnts = 'x-small'
    ax.annotate("80", xy=(0, -10), ha="left", va="bottom", fontsize=fnts)
    ax.annotate("60", xy=(0, -30), ha="left", va="bottom", fontsize=fnts)
    # add mlt labels
    ax.annotate("0", xy=(0, -rmax), ha="center", va="top", fontsize=fnts)
    ax.annotate("6", xy=(rmax, 0), ha="left", va="center", fontsize=fnts)
    ax.annotate("18", xy=(-rmax, 0), ha="right", va="center", fontsize=fnts)

    return  ccoll

def add_cbar(fig, coll, bounds, label="Number of Measurements", cax=None):

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
        l.append(str(int(bounds[i])))
    cbar.ax.set_yticklabels(l)
    #cbar.ax.tick_params(axis='y',direction='out')
    cbar.set_label(label)

    return

def by_season():

    import datetime as dt
    import matplotlib.pyplot as plt
    from convection import fetch_data
    import matplotlib as mpl


    # input parameters
    #nvel_min=250
    nvel_min=300
    lat_range=[52, 59]
    #lat_range=[40, 60]
    #lat_min = 43
    lat_min = 50

    ftype = "fitacf"
    #ftype = "fitex"

    seasons = ["winter", "summer", "equinox"]
    #seasons = ["winter", "summer"]
    fig_dir = "./plots/num_measurement_points/kp_l_3/data_in_mlt/"
    fig_name = "seasonal_num_measurement_points_2"
   
    # create subplots
    fig, axes = plt.subplots(nrows=len(seasons), ncols=1, figsize=(6,8))
    fig.subplots_adjust(hspace=0.3)

    # build a custom color map and bounds
    color_list = ['purple', 'b', 'dodgerblue', 'c', 'g', 'y', 'orange', 'r']
    cmap = mpl.colors.ListedColormap(color_list)
    bounds = range(0, 4000, 500)
    bounds[0] = 300
    bounds.append(10000)


    if len(seasons) == 1:
        axes = [axes]

    for i, season in enumerate(seasons):
        # fetches the data from db 
        #baseLocation="../data/sqlite3/" + season + '/'
        #baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_geo/'
        baseLocation="../data/sqlite3/" + season + '/kp_l_3/data_in_mlt/'

        data_dict = fetch_data(ftype=ftype, nvel_min=nvel_min, 
                        lat_range=lat_range, baseLocation=baseLocation,
                        dbName=None)

        # plot the flow vectors
        coll = num_plot(axes[i], data_dict, season, cmap, bounds, lat_min=lat_min)

    # add colorbar
    fig.subplots_adjust(right=0.78)
    cbar_ax = fig.add_axes([0.83, 0.25, 0.02, 0.5])
    add_cbar(fig, coll, bounds, cax=cbar_ax, label="Number of Measurements")

    # save the fig
    fig.savefig(fig_dir + fig_name + ".png", dpi=300)

#    plt.show()

    return fig

def by_imf_clock_angle():

    import datetime as dt
    import matplotlib.pyplot as plt
    from convection import fetch_data
    import matplotlib as mpl

    # input parameters
    #nvel_min=300
    nvel_min=100
    lat_range=[52, 59]
    #lat_range=[40, 60]
    #lat_min = 43
    lat_min = 50

    ftype = "fitacf"
    #ftype = "fitex"

    #seasons = ["winter", "summer", "equinox"]
    seasons = ["winter"]
   
    # build a custom color map and bounds
    color_list = ['purple', 'b', 'dodgerblue', 'c', 'g', 'y', 'orange', 'r']
    cmap = mpl.colors.ListedColormap(color_list)
    bounds = range(100, 500, 50)
    #bounds[0] = 100
    bounds.append(10000)

        
    # imf clock angle bins
    #bins = [[65, 115], [245, 295]] 
    bins = [[300, 60], [120, 240]]
    bins_txt = ["Bz+", "Bz-"]

    for season in seasons:

        # create subplots
        fig, axes = plt.subplots(nrows=len(bins), ncols=1, figsize=(6,8))
        fig.subplots_adjust(hspace=0.3)
        if len(bins) == 1:
            axes = [axes]

        fig_dir = "./plots/cosfit_plot/kp_l_3/data_in_mlt/binned_by_imf_clock_angle/"
        fig_name = season + "_binned_by_imf_Bz_num_of_points_" + str(lat_range[0]) +"_to_lat" + str(lat_range[1])

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
            title = "Number of Measurements, " + season[0].upper()+season[1:] +\
                    r", Kp $\leq$ 2+" + ", IMF " + bins_txt[i]
            coll = num_plot(axes[i], data_dict, season, cmap, bounds,
                            lat_min=lat_min, title=title)

        # add colorbar
        fig.subplots_adjust(right=0.78)
        cbar_ax = fig.add_axes([0.83, 0.25, 0.02, 0.5])
        add_cbar(fig, coll, bounds, cax=cbar_ax, label="Number of Measurements")

        # save the fig
        #fig.savefig(fig_dir + fig_name + ".png", dpi=300)
        plt.show()

        return


if __name__ == "__main__":
    #by_season()
    by_imf_clock_angle()
