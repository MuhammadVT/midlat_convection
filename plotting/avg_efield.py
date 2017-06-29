def avg_efield_plot(efield_dir="Northward", center_at_zero_mlt=True, dmlt=0.5):
    
    """ plots the average efield vs MLT
    dmlt : float
        the mlt hour of the edge of a mlt-mlat grid to its center mlt
    """

    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection,LineCollection
    import numpy as np
    import datetime as dt
    from convection import fetch_data
    from matplotlib.ticker import MultipleLocator

    # Earth Bz component
    Bz = 40 * 1.0e-6

    # input parameters
    nvel_min=250
    lat_range=[52, 59]
    ftype = "fitacf"
    #ftype = "fitex"
    if efield_dir == "Northward":
        veldir="zonal"
    elif efield_dir == "Eastward":
        veldir="meridional"

    seasons = ["winter", "summer", "equinox"]
    #seasons = ["winter", "summer"]

    # create a subplot
    fig, ax = plt.subplots()
    # colors of different lines
    color_list = ['b', 'g', 'r']
    # MLATs
    glatc = np.arange(lat_range[0], lat_range[1], 1) + 0.5

    for i, season in enumerate(seasons):
        # fetches the data from db 
        baseLocation="../data/sqlite3/" + season + '/'
        data_dict = fetch_data(ftype=ftype, nvel_min=nvel_min, 
                        lat_range=lat_range, baseLocation=baseLocation,
                        dbName=None)
        # calculate velocity components
        vel_mag = data_dict['vel_mag']
        vel_dir = np.deg2rad(data_dict['vel_dir'])

        if veldir == "zonal":
            vel_comp = vel_mag*(-1.0)*np.sin(vel_dir)
        elif veldir == "meridional":
            vel_comp = vel_mag*(-1.0)*np.cos(vel_dir)

        vel_mlt = data_dict['glonc'] / 15.
        vel_mlat = data_dict['glatc']
        efield = -1.0 *  Bz * vel_comp * 1.0e3

        avg_efield = []
        xdata = range(24)
        for m in xdata:
            efld = [] 
            for lat in glatc:
                indx_tmp = [j for j in range(len(vel_mlt)) \
                        if (np.abs(vel_mlt[j] - m) <=dmlt) and vel_mlat[j] == lat]
                if indx_tmp:
                    efld.append(np.mean(efield[indx_tmp]))
            if efld:
                # take the latitidinal average
                avg_efield.append(np.mean(efld))
            else:
                avg_efield.append(None)

        # plot the flow vector components
        xdata.append(24)
        avg_efield.append(avg_efield[0])
        if center_at_zero_mlt:
            # center at 0 MLT
            xdata = [x if (x is None) or x <=12 else x-24 for x in xdata]
        ax.plot(xdata, avg_efield, marker='o', color=color_list[i], label=season)

    # set axis label
    ax.set_xlabel("MLT")
    ax.xaxis.set_major_locator(MultipleLocator(base=3))
    if center_at_zero_mlt:
        xlabels = [str(x) for x in range(12, 24, 3) + range(0, 15, 3)]
        plt.xticks(range(-12, 15, 3), xlabels)


    # add text
    ax.set_title("Average " + efield_dir + " Electric Field", fontsize="small")

    # add zero-line
    ax.axhline(y=0.0, color='k', linewidth=0.7)

    # set axis limits
    if center_at_zero_mlt:
        ax.set_xlim([-12, 12])
        # add legend
        ax.legend(bbox_to_anchor=(1.00, 0.96))
    else:
        ax.set_xlim([0, 24])
        # add legend
        ax.legend(bbox_to_anchor=(0.7, 0.96))
    ax.set_ylim([-0.5, 2.0])
    
    # axis labels
    ax.set_ylabel("Electric field [mV/m]")

    return fig

def main():
    
    import matplotlib.pyplot as plt
    
    efield_dir = "Northward"
    center_at_zero_mlt = True 
    #center_at_zero_mlt = False 
    dmlt=0.1

    fig = avg_efield_plot(efield_dir=efield_dir, center_at_zero_mlt=center_at_zero_mlt, dmlt=dmlt)

    # save the fig
    if center_at_zero_mlt:
        fig.savefig("./plots/efield/" + efield_dir+ "_vel_vs_mlt_dmlt" + str(dmlt) + "_c0.png", dpi=300)
    else:
        fig.savefig("./plots/efield/" + efield_dir+ "_vel_vs_mlt_dmlt" + str(dmlt) + ".png", dpi=300)
    #plt.show()

if __name__ == "__main__":
    main()
