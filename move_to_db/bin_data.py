class grids(object):
    """ This class is used to create grid points """

    def __init__(self, lat_min=50, lat_max=90, dlat=1, half_dlat_offset=False):
        '''
        "half_dlat_offset=False" implements NINT[360 sin(theta)] at theta = 89, 88, ... colatitude
        "half_dlat_offset=True" implements NINT[360 sin(theta)] at theta = 89.5, 88.5, ... colatitude
        '''

        import numpy as np

        self.lat_min = lat_min
        self.lat_max = lat_max
        self.dlat = dlat
        self.half_dlat_offset = half_dlat_offset
        self.center_lats = [x + 0.5*dlat for x in range(lat_min, lat_max, dlat)] 
        if half_dlat_offset:
            self.nlons = [round(360 * np.sin(np.deg2rad(90-lat))) for lat in self.center_lats]
        else:
            self.nlons = [round(360 * np.sin(np.deg2rad(90-(lat-0.5*dlat)))) for lat in self.center_lats]
        self.dlons = [360./nn for nn in self.nlons]

        # lat and lon bins
        self.lat_bins = [x for x in np.arange(lat_min,lat_max+dlat,dlat)] 
        self.lon_bins, self.center_lons = self._create_lonbins()
        
        # azimuthal bins and their centers
        # zero azm directs to the mag north
        self.azm_bins = [x for x in range(0, 370, 10)]
        self.center_azms = [x for x in range(5, 365, 10)]


    def _create_lonbins(self):
        import numpy as np
        lon_bins = []
        center_lons = []      # list of lists of lons
        for i in range(len(self.center_lats)):
            lon_tmp = [ round(item*self.dlons[i],2) for item in np.arange(0.5, self.nlons[i]+0.5) ]
            center_lons.append(lon_tmp)
            lon_tmp = [ item*self.dlons[i] for item in np.arange(self.nlons[i]) ]
            lon_tmp.append(360) 
            lon_bins.append(lon_tmp)

        return lon_bins, center_lons 
        

def bin_to_grid(rad, bmnum, stm=None, etm=None, ftype="fitacf",
               dbName=None, baseLocation="../data/sqlite3/"):
    """ bins the data into mlat-mlt-azm grid
    Note
    ----
        0 gazmc directs towards magnetic north. 180 gazmc directs towards south.
        gazmc spans from 5 - 355 degrees. 
    
    """

    import sys
    sys.path.append("../")
    from dbtools.db.connection.Connector import Connectors
    import numpy as np

    # create grid points
    grds = grids(lat_min=35, lat_max=90, dlat=1, half_dlat_offset=False)

    # make a db connection
    if dbName is None:
        dbName = "dopsearch_" + rad + "_" + ftype + ".sqlite"
    conn = Connectors(baseLocation = baseLocation, dbName = dbName,
                          isInMem = False, isAutoSave = False)

    table_name = rad + "_bm" + str(bmnum)

    # add new columns
    try:
        command ="ALTER TABLE {tb} ADD COLUMN glatc TEXT".format(tb=table_name) 
        conn.cursor.execute(command)
    except:
        # pass if the column latc exists
        pass
    try:
        command ="ALTER TABLE {tb} ADD COLUMN glonc TEXT".format(tb=table_name) 
        conn.cursor.execute(command)
    except:
        # pass if the column lonc exists
        pass
    try:
        command ="ALTER TABLE {tb} ADD COLUMN gazmc TEXT".format(tb=table_name) 
        conn.cursor.execute(command)
    except:
        # pass if the column latc exists
        pass

    # do the convertion to all the data in db if stm and etm are all None
    if stm is not None and etm is not None:
        command = "SELECT rowid, latc, lonc, azm, datetime FROM {tb} WHERE (DATETIME(datetime)>='{stm}' and\
                   DATETIME(datetime)<='{etm}') ORDER BY datetime".format(tb=table_name,\
                   stm=str(stm), etm=str(etm))

    # do the convertion to the data between stm and etm if they are not None
    else:
        command = "SELECT rowid, latc, lonc, azm, datetime FROM {tb} ORDER BY datetime".format(tb=table_name,\
                   sdtm=str(stm), edtm=str(etm))
    conn.cursor.execute(command)
    rows = conn.cursor.fetchall() 

    # do the conversion row by row
    if rows != []:
        for row in rows:
            rowid, latc, lonc, azm, date_time= row
            if latc:

                # convert string to a list of float
                latc = [float(x) for x in latc.split(",")]
                lonc = [float(x) for x in lonc.split(",")]
                azm = [float(x) for x in azm.split(",")]

                # grid the data
                # grid latc
                indx_latc = np.digitize(latc, grds.lat_bins)
                indx_latc = [x-1 for x in indx_latc]
                glatc = [grds.center_lats[x] for x in indx_latc]

                # grid lonc
                indx_lonc = [np.digitize(lonc[i], grds.lon_bins[indx_latc[i]]) 
                             for i in range(len(lonc))]
                indx_lonc = [x-1 for x in indx_lonc]
                glonc = [grds.center_lons[indx_latc[i]][indx_lonc[i]]\
                         for i in range(len(lonc))]

                # grid azm
                indx_azmc = np.digitize(azm, grds.azm_bins)
                indx_azmc = [x-1 for x in indx_azmc]
                gazmc = [grds.center_azms[x] for x in indx_azmc]

                # convert to comma seperated text
                glatc =",".join([str(x) for x in glatc])
                glonc =",".join([str(x) for x in glonc])
                gazmc =",".join([str(x) for x in gazmc])

                # update the table
                command = "UPDATE {tb} SET glatc='{glatc}', glonc='{glonc}', gazmc='{gazmc}'\
                          WHERE rowid=={rowid}".format(tb=table_name, glatc=glatc,\
                                                       glonc=glonc, gazmc=gazmc, rowid=rowid)
                conn.cursor.execute(command)

        # commit the data into the db
        conn._commit()
    # close the db connection
    conn.connection.close()
    return

# test code
def main():

    import datetime as dt
    # input parameters
    #rad_list = ["bks", "wal", "fhe", "fhw", "cve", "cvw"]
    rad_list = ["bks"]
    #rad_list = ["wal"]
    ftype = "fitacf"
    bmnum = 7
    stm = dt.datetime(2012,1,1)
    etm = dt.datetime(2012,1,31)
    ftype = "fitacf"

    for rad in rad_list:
        bin_to_grid(rad, bmnum, stm=stm, etm=etm, ftype=ftype, dbName=None)
    return
if __name__ == "__main__":
    main()

