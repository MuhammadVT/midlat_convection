def pol2cart(phi, rho):
    import numpy as np
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def plot_grid(lat_min=60, lat_max=90, dlat=1, coords='geo', half_dlat_offset=True, date_time = None):
    import numpy as np
    import matplotlib.pyplot as plt 
    import matplotlib.patches as patches
    from davitpy import utils
    import datetime as dt
    from matplotlib.collections import PolyCollection, LineCollection
    import sys
    sys.path.append("../move_to_db/")
    from bin_data import grids 

    if not date_time: 
        date_time = dt.now()

#    lats = [x + 0.5*dlat for x in range(lat_min,lat_max,dlat)] 
#    if half_dlat_offset:
#        nlons = [round(360 * np.sin(np.deg2rad(90-lat))) for lat in lats]
#    else:
#        nlons = [round(360 * np.sin(np.deg2rad(90-(lat-0.5*dlat)))) for lat in lats]
#    dlons = [360./nn for nn in nlons]

    # create grid points
    grds = grids(lat_min=lat_min, lat_max=lat_max, dlat=dlat, half_dlat_offset=half_dlat_offset)
    lats = grds.center_lats
    nlons = grds.nlons
    dlons = grds.dlons

    # flatting all lats and lons 
    lons_all = np.array([])
    lons_all_E = np.array([])
    lons_all_W = np.array([])
    lats_all = np.array([])
    lats_all_N = np.array([])
    lats_all_S = np.array([])
    for i in range(len(lats)):
        lons = [ item*dlons[i] for item in np.arange(0.5, (nlons[i]+0.5)) ]
        lons_E = [ item*dlons[i] for item in np.arange(1, (nlons[i]+1)) ]
        lons_W = [ item*dlons[i] for item in np.arange(0, (nlons[i])) ]
        lons_all = np.append(lons_all, lons)
        lons_all_E = np.append(lons_all_E, lons_E)
        lons_all_W = np.append(lons_all_W, lons_W)

        lats_all = np.append(lats_all, np.repeat(lats[i], nlons[i]))
        lats_all_N = np.append(lats_all_N, np.repeat((lats[i]+dlat*0.5), nlons[i]))
        lats_all_S = np.append(lats_all_S, np.repeat((lats[i]-dlat*0.5), nlons[i]))

    # plot grid
    fig, ax = plt.subplots(figsize=(6,6))
#    my_map = utils.mapObj(coords=coords, projection='stere', boundinglat=lat_min+10, fillContinents='w', datetime=date_time)
#    x1,y1 = my_map(lons_all_W,lats_all_S)
#    x2,y2 = my_map(lons_all_W,lats_all_N)
#    x3,y3 = my_map(lons_all_E,lats_all_N)
#    x4,y4 = my_map(lons_all_E,lats_all_S)


    x1,y1 = pol2cart(np.deg2rad(lons_all_W), 90-lats_all_S)
    x2,y2 = pol2cart(np.deg2rad(lons_all_W), 90-lats_all_N)
    x3,y3 = pol2cart(np.deg2rad(lons_all_E), 90-lats_all_N)
    x4,y4 = pol2cart(np.deg2rad(lons_all_E), 90-lats_all_S)

    verts = zip(zip(x1,y1), zip(x2,y2), zip(x3,y3), zip(x4,y4))

    # set axis limits
    rmax = 90 - lat_min - 15
    ax.set_xlim([-rmax, rmax])
    ax.set_ylim([-rmax, rmax])
    ax.set_aspect("equal")

    # remove tikcs
    ax.tick_params(axis='both', which='both', bottom='off', top='off',
                   left="off", right="off", labelbottom='off', labelleft='off')


    # plot the latitudinal circles
    for r in [10, 30, 50]:
        c = plt.Circle((0, 0), radius=r, fill=False)
        ax.add_patch(c)

    # plot the longitudinal lines
    #for l in np.deg2rad(np.array([210, 240, 270, 300, 330])):
    for l in np.deg2rad(np.array(range(0, 360, 30))):
        x1_l, y1_l = pol2cart(l, 10) 
        x2_l, y2_l = pol2cart(l, 50) 
        ax.plot([x1_l, x2_l], [y1_l, y2_l], 'k')

    pcoll = PolyCollection(np.array(verts), linewidth=0.3, facecolor='', zorder=5)
    ax.add_collection(pcoll)

    # fill in one of the bins with black
    aa = np.transpose(np.array([x1, y1]))
    indx = np.argmin(np.linalg.norm(aa - np.array([[3, -22]]), axis=1))
    vert = verts[indx:indx+1]
    pcoll = PolyCollection(np.array(vert), linewidth=0.3, facecolor='k', zorder=5)
    ax.add_collection(pcoll)

    # zoon in a mlat-mlt-az grid
    ctr = 23
    dl = 5
    xs = [ctr-dl, ctr-dl, ctr+dl, ctr+dl]
    ys = [-dl, dl, dl, -dl]
    verts2 = [tuple(zip(xs,ys))]
    pcoll = PolyCollection(np.array(verts2), linewidth=0.8, facecolor='', zorder=6)
    ax.add_collection(pcoll)

    vert_x = [vert[0][0][0], vert[0][1][0], vert[0][2][0], vert[0][3][0]]
    vert_y = [vert[0][0][1], vert[0][1][1], vert[0][2][1], vert[0][3][1]]
    for i in [0, 1, 3]:
        ax.plot([vert_x[i], xs[i]], [vert_y[i], ys[i]], 'k', linewidth=0.7)

    #axis_to_data = ax.transAxes + ax.transData.inverted()
    axis_to_data = fig.transFigure + ax.transData.inverted()
    data_to_axis = axis_to_data.inverted()
    axis_points = (data_to_axis.transform((xs[0], ys[0])),
                   data_to_axis.transform((xs[2], ys[2])))
    ax1_box = [axis_points[0][0], axis_points[0][1],
               axis_points[1][0]-axis_points[0][0],axis_points[1][1]-axis_points[0][1]]
    ax1 = fig.add_axes(ax1_box)
    ax1.set_aspect("equal")
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_xlim([-1,1])
    ax1.set_ylim([-1,1])

    # plot the azimuthal bins
    Xdat = []
    Ydat = []
    for az in range(-40, 50, 10):
        xx = 1 
        yy = xx * np.tan(np.deg2rad(az))
        Xdat.append(xx)
        Ydat.append(yy)
        ax1.plot([0, xx], [0, yy],   'k', linewidth=0.8)
        ax1.plot([0, -xx], [0, -yy], 'k', linewidth=0.8)
        ax1.plot([0, yy], [0, xx],   'k', linewidth=0.8)
        ax1.plot([0, -yy], [0, -xx], 'k', linewidth=0.8)


    # add latitudinal labels
    fnts = 'x-small'
    ax.annotate("80", xy=(0, -10), ha="left", va="bottom", fontsize=fnts)
    ax.annotate("60", xy=(0, -30), ha="left", va="bottom", fontsize=fnts)
    # add mlt labels
#    ax.annotate("0", xy=(0, -rmax), ha="center", va="top", fontsize=fnts)
#    ax.annotate("6", xy=(rmax, 0), ha="left", va="center", fontsize=fnts)
#    ax.annotate("18", xy=(-rmax, 0), ha="right", va="center", fontsize=fnts)
    ax.annotate("0", xy=(0.5, -0.005), xycoords="axes fraction", ha="center", va="top", fontsize=fnts)
    ax.annotate("6", xy=(1.005, 0.5), xycoords="axes fraction", ha="left", va="center", fontsize=fnts)
    ax.annotate("12", xy=(0.5, 1.005), xycoords="axes fraction", ha="center", va="bottom", fontsize=fnts)
    ax.annotate("18", xy=(-0.005, 0.5), xycoords="axes fraction", ha="right", va="center", fontsize=fnts)

    return fig

if __name__=="__main__":

    import matplotlib.pyplot as plt 
    import datetime as dt
    stime = dt.datetime(2014,9,12)
    lat_min=40; lat_max=90; dlat= 1
    coords = 'mlt'
    #coords = 'geo'
    fig = plot_grid(lat_min=lat_min, lat_max=lat_max, dlat=1, coords=coords,
        half_dlat_offset=False, date_time = stime)
    
    fig.savefig("./plots/grid_plot/grid_plot.png", dpi=300)
    #plt.show()

