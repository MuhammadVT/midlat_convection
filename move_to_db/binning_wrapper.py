def wrapper(rad, bmnum, stm, etm, ftype="fitacf", season=None, baseLocation = "../data/sqlite3/",
            stay_in_geo=False):
    """ a wrapper function that does the following to the dopsearch_xxxx.sqlite:
        1. calculate center positions (latc and lonc) in geo coords and inser them into db
        2. convert latc and lonc from geo into mlt coords, and also calculate LoS vel azm relative
           to mag north 
        3. bin the data into mlat-mlt-azm grid
    """

    from latc_lonc_to_db import latc_lonc_to_db
    from geo_to_mlt import geo_to_mlt
    from bin_data import bin_to_grid

    # calculate center positions (latc and lonc) in geo coords and inser them into db
    latclonc = latc_lonc_to_db(rad, bmnum, stm, etm, coords="geo",
                               ftype=ftype, season=season, baseLocation=baseLocation)
    latclonc.add_latclonc_to_db()
    
    # convert latc and lonc from geo into mlt coords, and also calculate LoS vel azm relative
    # to mag north 

    geo_to_mlt(rad, bmnum, stm=stm, etm=etm, ftype=ftype,
               dbName=None, baseLocation=baseLocation, t_c_alt=0., stay_in_geo=stay_in_geo)

    # bin the data into mlat-mlt-azm grid
    bin_to_grid(rad, bmnum, stm=stm, etm=etm, ftype=ftype,
                dbName=None, baseLocation=baseLocation)

    return


def worker(rads, season, baseLocation, stay_in_geo=False):

    import sys
    sys.path.append("../")
    from dbtools.db.connection.Connector import Connectors
    import datetime as dt

    # input parameters
    ftype = "fitacf"

    # set the time interval for each season
    if season == "summer":
        stms = [dt.datetime(2011,5,1), dt.datetime(2012,5,1)]
        etms = [dt.datetime(2011,9,1), dt.datetime(2012,9,1)]

    if season == "winter":
        stms = [dt.datetime(2011,1,1), dt.datetime(2011,11,1), dt.datetime(2012,11,1)]
        etms = [dt.datetime(2011,3,1), dt.datetime(2012,3,1), dt.datetime(2013,1,1)]

    if season == "equinox":
        stms = [dt.datetime(2011,3,1), dt.datetime(2011,9,1), dt.datetime(2012,3,1), dt.datetime(2012,9,1)]
        etms = [dt.datetime(2011,5,1), dt.datetime(2011,11,1), dt.datetime(2012,5,1), dt.datetime(2012,11,1)]

    # loop through time interval
    for dd in range(len(stms)):
        stm = stms[dd]
        etm = etms[dd]
                    
        # loop through radars
        for rad in rads:
            dbName = "dopsearch_" + rad + "_" + ftype + ".sqlite"

            # make db connection
            conn = Connectors(baseLocation = baseLocation, dbName = dbName,
                              isInMem = False, isAutoSave = False)

            # get all the table names
            conn.cursor.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
            tbl_names = conn.cursor.fetchall()
            tbl_names = [x[0] for x in tbl_names]

            # close db connection
            conn.connection.close()

            # get the available beam numbers 
            beam_nums = [x.split("_")[-1][2:] for x in tbl_names]
            beam_nums = sorted([int(x) for x in beam_nums])
        
            # loop through each table
            for bmnum in beam_nums:

                print "start updating " + dbName + " for beam " + str(bmnum) + " between " + str(stm) + " and " + str(etm)
                wrapper(rad, bmnum, stm, etm, ftype=ftype, season=season,
                        baseLocation=baseLocation, stay_in_geo=stay_in_geo)
                print "finish updating " + dbName + " for beam " + str(bmnum)

    return
if __name__ == "__main__":
    import multiprocessing
    
    # NOTE: run one season at a time

    #seasons = ["winter", "summer", "equinox"]
    seasons = ["winter"]
    #rads_list = [["bks", "wal", "fhe"], ["fhw", "cve", "cvw"]]
    rads_list = [["bks", "wal", "fhe", "fhw", "cve", "cvw"]]

    #stay_in_geo = False
    stay_in_geo = True
    
    jobs = []
    for season in seasons:
        if stay_in_geo:
            #baseLocation="../data/sqlite3/" + season + "/" + "data_in_geo" + "/"
            baseLocation="../data/sqlite3/" + season + "/lat_lon_az_binned_data/" + \
                         "low_vel_iscat_event_only/" + "data_in_geo" + "/"
        else:
            baseLocation="../data/sqlite3/" + season + "/lat_lon_az_binned_data/" + \
                         "low_vel_iscat_event_only/" + "data_in_mlt" + "/"
        for i in range(len(rads_list)):
#            worker(rads_list[i], season, baseLocation, stay_in_geo=stay_in_geo)
            p = multiprocessing.Process(target=worker, args=(rads_list[i], season, baseLocation, stay_in_geo))
            jobs.append(p)
            p.start()
