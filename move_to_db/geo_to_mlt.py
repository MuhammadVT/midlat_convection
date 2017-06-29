# radar locations in geo coords
rad_loc_dict = {"bks" : (37.10, -77.95), "wal" : (37.93, -75.47), "fhe" : (38.86, -99.39),
            "fhw" : (38.86, -99.39), "cve" : (43.27, -120.36), "cvw" : (43.27, -120.36)}

def geo_to_mlt(rad, bmnum, stm=None, etm=None, ftype="fitacf",
               dbName=None, baseLocation="../data/sqlite3/",
               t_c_alt=0., stay_in_geo=False):
    """ converts latc and lonc from geo to mlt hours. Also
    calcuates the azmimuthal velocity angle relative to the magnetic pole

    stay_in_geo : bool
        if set to True no coord conversion is done. Calculation would be in geo
    """

    from davitpy.utils.coordUtils import coord_conv
    import sys
    sys.path.append("../")
    from dbtools.db.connection.Connector import Connectors
    from datetime import date
    import datetime as dt
    import pdb


    # make a db connection
    if dbName is None:
        dbName = "dopsearch_" + rad + "_" + ftype + ".sqlite"
    conn = Connectors(baseLocation = baseLocation, dbName = dbName,
                          isInMem = False, isAutoSave = False)

    table_name = rad + "_bm" + str(bmnum)

    # add new columns
    try:

        # add the azmimuthal velocity angle (azm) relative to the magnetic pole
        command ="ALTER TABLE {tb} ADD COLUMN azm TEXT".format(tb=table_name) 
        conn.cursor.execute(command)
    except:
        # pass if the column latc exists
        pass


    # do the convertion to all the data in db if stm and etm are all None
    if stm is not None and etm is not None:
        command = "SELECT rowid, latc, lonc, bmazm, datetime FROM {tb} WHERE (DATETIME(datetime)>='{sdtm}' and\
                   DATETIME(datetime)<='{edtm}') ORDER BY datetime".format(tb=table_name,\
                   sdtm=str(stm), edtm=str(etm))

    # do the convertion to the data between stm and etm if they are not None
    else:
        command = "SELECT rowid, latc, lonc, bmazm, datetime FROM {tb} ORDER BY datetime".format(tb=table_name)
    conn.cursor.execute(command)
    rows = conn.cursor.fetchall() 

    # do the conversion row by row
    if rows:
        for row in rows:
            rowid, latc, lonc, bmazm, date_time= row
            if latc:
                # convert string to a list of float
                latc = [float(x) for x in latc.split(",")]
                lonc = [float(x) for x in lonc.split(",")]

                # calculate bmazm in mag. The return value is a command seperated strings
                azm_txt = geobmazm_to_magbmazm(rad, bmazm, latc, lonc, alt=300.,
                                               time=date_time.date(), stay_in_geo=stay_in_geo)

                if stay_in_geo:
                    lonc = [x%360 for x in lonc]
                    # convert utc to local time in degrees
                    lonc_ltm = []
                    for lonc_i in lonc:
                        lonc_tmp = lonc_i if lonc_i<=180 else lonc_i-360
                        # convert utc to local time
                        local_dt = date_time + dt.timedelta(hours=lonc_tmp/15.)
                        ltm = local_dt.time()  
                        # convert local time to degrees. e.g. 0 (or 360) degree is midnight, 180 degrees is noon time. 
                        lonc_ltm.append((ltm.hour + ltm.minute/60. + ltm.second/3600.) * 15.)
                    lonc = lonc_ltm
                if not stay_in_geo:
                    # convert from geo to mlt degress
                    lonc, latc = coord_conv(lonc, latc, "geo", "mlt",
                                            altitude=t_c_alt,
                                            date_time=date_time)
                    ## convert mlt degress to mlt hours
                    #lonc = [(x%360)/15. for x in lonc]

                lonc = [(round(x,2))%360 for x in lonc]

                # convert to comma seperated text
                latc =",".join([str(x) for x in latc])
                lonc =",".join([str(round(x,2)) for x in lonc])
                
                # update into the db
                command = "UPDATE {tb} SET latc='{latc}', lonc='{lonc}', azm='{azm_txt}'\
                           WHERE rowid=={rowid}".\
                          format(tb=table_name, latc=latc, lonc=lonc, azm_txt=azm_txt, rowid=rowid)
                conn.cursor.execute(command)
            else:
                continue

    conn._commit()

    # close db connection
    conn.connection.close()
def geobmazm_to_magbmazm(rad, bmazm, latc, lonc, alt=300., time=None, stay_in_geo=False):
    """ calculates the LOS vel direction in mag coords at each range-beam cell
   
    bmazm : float
        bmazm of a certain beam in geo. 
        0 degree shows the mag north direction
        180 degree shows the mag south direction
    latc : list
    lonc : list
        center lat and on positions of range-beam cells along a beam
    alt : float
        altitude value
    stay_in_geo : bool
        if set to True no coord conversion is done. Calculation would be in geo
   
    Return
    ------
    azm_txt : string
        azm values at the positions of latc and lonc in mag coords. The values are converted
        to a commad seperated strings

    """
    
    from geomag import geomag
    from datetime import date
    import numpy as np
    import pdb
   
    rad_lat, rad_lon = rad_loc_dict[rad]
    rad_lon = rad_lon % 360
    azm_lst = []
    gm = geomag.GeoMag()
    for i in range(len(latc)):
        # calculate the los vel angle in geo using spherical trigonometry. Then angles are defined
        # in the same way as those in spherical trigonometry section in mathworld
        #B = np.deg2rad(np.abs(bmazm))
        b_prime = np.deg2rad(90. - latc[i])
        a_prime = np.deg2rad(90. - rad_lat)
        AB_dellon = np.deg2rad(np.abs(lonc[i]-rad_lon))
        c_prime = np.arccos(np.sin(np.deg2rad(rad_lat)) * np.sin(np.deg2rad(latc[i])) +\
                          np.cos(np.deg2rad(rad_lat)) * np.cos(np.deg2rad(latc[i])) * np.cos(AB_dellon))
        s_prime = 1./2 * (a_prime + b_prime + c_prime)
        if round(np.rad2deg(a_prime),5) == round(np.rad2deg(s_prime),5):
            A = np.pi
        else:
            A = 2 * np.arcsin(np.sqrt((np.sin(s_prime - b_prime) * np.sin(s_prime - c_prime)) /\
                                  (np.sin(b_prime) * np.sin(c_prime))))
        losvel_dir = np.sign(bmazm) * (180 - np.rad2deg(A))
        
        if stay_in_geo:
            azm_mag = (round(losvel_dir,2)) % 360
        else:
            # convert from geo to mag by add the magnetic declanation angle to the los vel angle in geo
            mg = gm.GeoMag(latc[i], lonc[i], h=alt, time=time)
            azm_mag = (round(losvel_dir - mg.dec,2)) % 360
#        if str(azm_mag) == "nan":
#            pdb.set_trace()
        azm_lst.append(azm_mag)

    # convert the list entries to a comma seperated strings
    azm_txt =",".join([str(x) for x in azm_lst])

    return azm_txt



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
    coords = "geo"
    ftype = "fitacf"
    stay_in_geo=False    # set this to True if you want to remain in "geo" coords

    for rad in rad_list:
        geo_to_mlt(rad, bmnum, stm=None, etm=None, ftype="fitacf",
                   dbName=None, baseLocation="../data/sqlite3/",
                   t_c_alt=0., stay_in_geo=stay_in_geo)
if __name__ == "__main__":
    main()

