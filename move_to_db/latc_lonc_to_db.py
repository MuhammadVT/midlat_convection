"""
written by Muhammad Rafiq, 2016-08-13
"""

import datetime as dt
from davitpy.pydarn.radar.radFov import slantRange, calcFieldPnt, calcAzOffBore
from davitpy.pydarn.radar.radStruct import site 
from davitpy.utils.coordUtils import coord_conv

import sys
sys.path.append("../")
from dbtools.db.dao.DBAccessObject import Table as T, TableDescriptor as TD, ColumnTypes as CT
from dbtools.db.connection.Connector import Connectors, Executor
import pdb


class latc_lonc_to_db(object):

    """ calculates the center points of range-beam cells and move them into the original dopsearch db""" 

    def __init__(self, rad, bmnum, stm, etm, coords="geo",
                 ftype="fitacf", season=None, baseLocation="../data/sqlite3/"):

        rad_id_dict = {"bks":33, "wal":32, "cve":207, "cvw":206, "fhe":205, "fhw":204}
        self.rad = rad
        self.rad_id = rad_id_dict[rad]
        self.bmnum = bmnum
        self.coords = coords
        self.stm = stm
        self.etm = etm
        self.season = season
        self.baseLocation = baseLocation
        self.table_name = self.rad + "_bm" + str(self.bmnum)
        self.ftype = ftype
        self.conn = self._create_dbconn()
        self.TableDescriptor = self._create_TableDescriptor()
        self.sites = self._create_site_list()
        #self.input_dts = self._group_db_input()

    def _create_dbconn(self, dbName=None):

        # make a db connection
        if dbName is None:
            dbName = "dopsearch_" + self.rad + "_" + self.ftype + ".sqlite"
        conn = Connectors(baseLocation = self.baseLocation, dbName = dbName,
                          isInMem = False, isAutoSave = False)
        return conn


    def _create_site_list(self):

        """ creats a list of sites for a given self.rad between self.stm and self.etm """

        connector = Connectors(baseLocation = "../data/sqlite3/", dbName = "radars.sqlite",
                          isInMem = False, isAutoSave = False)
        #command = '{:s}and tval>=? and tval<=? ORDER BY tval ASC'.format(command)
        #connector.cursor.execute(command, (self.rad_id, self.stm, self.etm))
        #rows = connector.cursor.fetchall()
        command = "SELECT tval FROM hdw WHERE id=? "
        command = '{:s}and tval>=? ORDER BY tval ASC'.format(command)
        connector.cursor.execute(command, (self.rad_id, self.stm))
        tvals_stm = connector.cursor.fetchall()
        tvals_stm = [x[0] for x in tvals_stm]

        command = "SELECT tval FROM hdw WHERE id=? "
        command = '{:s}and tval>=? ORDER BY tval ASC'.format(command)
        connector.cursor.execute(command, (self.rad_id, self.etm))
        tval_etm = connector.cursor.fetchone()[0]
        indx_etm = tvals_stm.index(tval_etm)

        # select the tvals of interest
        tvals = tvals_stm[:indx_etm+1]

        site_list = []
        for tval in tvals:
            site_list.append(site(code=self.rad, dt=tval))
        return site_list

    def _create_TableDescriptor(self):
        """ creates a TableDescriptor object """

        self.conn.cursor.execute("PRAGMA table_info(" + self.table_name + ")")
        descriptions = self.conn.cursor.fetchall()
        column_map = {}
        for description in descriptions:
            column_map[description[1]] = description[2]
        td = TD(self.table_name, column_map) 

#        self.conn.cursor.execute("SELECT sql FROM sqlite_master WHERE name='{tb}'"\
#                                 .format(tb=self.table_name))
#        aa = str(self.conn.cursor.fetchone()[0])
#        sindx = aa.find("(")
#        eindx = aa.find(")")
#        aa = aa[sindx+1:eindx]
#        aa = aa.split(",")
#        column_map = {kyval.split()[0]:kyval.split()[1] for kyval in aa}
#        td = TD(self.table_name, column_map) 

        return td 

#    def _group_db_input(self):
#
#        # select variables
#        command = "SELECT datetime FROM (SELECT * FROM {tb} ORDER BY datetime ASC)\
#                   GROUP BY frang, rsep".format(tb=self.table_name)
#        #executor._execute_fetch_command()
#        self.conn.cursor.execute(command)
#        dts = self.conn.cursor.fetchall()
#        dts = [d[0] for d in dts]
#
#        return dts 
#
#
#    def add_to_db(self):
#        
#        sdtm = self.stm
#
#        for edtm in self.input_dts:
#            command = "SELECT frang, rsep FROM {tb} WHERE (DATETIME(datetime)>'{sdtm}' and\
#                       DATETIME(datetime)<='{edtm}')".format(tb=self.table_name,\
#                       sdtm=str(sdtm), edtm=str(edtm))
#            self.conn.cursor.execute(command)
#            rows = self.conn.cursor.fetchall() 
#            frang, rsep = rows[0]
#            latc, lonc = calc_latc_lonc(self.sites[0], self.bmnum, frang, rsep, altitude=300.,
#                                        elevation=None, coord_alt=0., coords="geo")
#        
#        #return rows
#        return latc, lonc

    def add_latclonc_to_db(self):
        """ calculates latc and lonc of each range-beam cell in 'geo' coordinates and update them
        into the original table """
        
        # add new columns
        try:
            command ="ALTER TABLE {tb} ADD COLUMN latc TEXT".format(tb=self.table_name) 
            self.conn.cursor.execute(command)
        except:
            # pass if the column latc exists
            pass
        try:
            command ="ALTER TABLE {tb} ADD COLUMN lonc TEXT".format(tb=self.table_name) 
            self.conn.cursor.execute(command)
        except:
            # pass if the column lonc exists
            pass

        # iterate through tvals of the self.sites
        sdtm = self.stm
        for ii, st in enumerate(self.sites):
            if ii == len(self.sites)-1:
                edtm = self.etm
            else:
                edtm = st.tval
            command = "SELECT rowid, slist, vel, frang, rsep, datetime FROM {tb} WHERE (DATETIME(datetime)>'{sdtm}' and\
                       DATETIME(datetime)<='{edtm}') ORDER BY datetime".format(tb=self.table_name,\
                       sdtm=str(sdtm), edtm=str(edtm))
            self.conn.cursor.execute(command)
            rows = self.conn.cursor.fetchall() 
            if rows != []:
                rowid, slist, vel, frang_old, rsep_old, date_time_old = rows[0]

                # calculate latc_all and lonc_all in 'geo' coords
                latc_all, lonc_all = calc_latc_lonc(self.sites[ii], self.bmnum, frang_old, rsep_old, 
                                                    altitude=300., elevation=None, coord_alt=0.,
                                                    coords="geo", date_time=None)
                for row in rows:
                    rowid, slist, vel, frang, rsep, date_time = row
                    if (frang, rsep) != (frang_old, rsep_old):
                        latc_all, lonc_all = calc_latc_lonc(self.sites[ii], self.bmnum, frang, rsep, 
                                                    altitude=300., elevation=None, coord_alt=0.,
                                                    coords="geo", date_time=None)

                        
                        frang_old, rsep_old = frang, rsep

                    # convert from string to float
                    slist = [int(float(x)) for x in slist.split(",")]
                    vel = [float(x) for x in vel.split(",")]

                    # exclude the slist values beyond maxgate and their correspinding velocities
                    vel = [vel[i] for i in range(len(vel)) if slist[i] < st.maxgate]
                    slist = [s for s in slist if s < st.maxgate]

                    # extract latc and lonc values
                    latc = [latc_all[s] for s in slist]
                    lonc = [lonc_all[s] for s in slist]

                    # convert to comma seperated text
                    slist = ",".join([str(x) for x in slist])
                    vel = ",".join([str(round(x,2)) for x in vel])
                    latc = ",".join([str(round(x,2)) for x in latc])
                    lonc = ",".join([str(round(x,2)) for x in lonc])

                    # update the table
                    command = "UPDATE {tb} SET slist='{slist}', vel='{vel}',\
                               latc='{latc}', lonc='{lonc}' WHERE rowid=={rowid}".\
                              format(tb=self.table_name, slist=slist, vel=vel,\
                              latc=latc, lonc=lonc, rowid=rowid)
                    self.conn.cursor.execute(command)

            # update sdtm
            sdtm = edtm

        # commit the data into the db
        self.conn._commit()

        # close db connection
        self.conn._close_connection()
            
        return


def calc_latc_lonc(site, bmnum, frang, rsep, altitude=300.,
                   elevation=None, coord_alt=0., coords="geo",
                   date_time=None):
    """ calculates center lat and lon of all the gates of a given bmnum """
    import numpy as np

    nbeams = site.maxbeam
    ngates = site.maxgate
    bmsep = site.bmsep
    recrise = site.recrise
    siteLat = site.geolat
    siteLon = site.geolon
    siteAlt = site.alt
    siteBore = site.boresite
    gates = np.arange(ngates)

    # Create output arrays
    lat_center = np.zeros(ngates, dtype='float')
    lon_center = np.zeros(ngates, dtype='float')

    # Calculate deviation from boresight for center of beam
    boff_center = bmsep * (bmnum - (nbeams - 1) / 2.0)

    # Calculate center slant range
    srang_center = slantRange(frang, rsep, recrise,
                              gates, center=True)

    # Calculate coordinates for Center of the current beam
    for ig in gates:
        talt = altitude
        telv = elevation
        t_c_alt = coord_alt

        # calculate projections
        latc, lonc = calcFieldPnt(siteLat, siteLon, siteAlt * 1e-3,
                                  siteBore, boff_center,
                                  srang_center[ig], elevation=telv,
                                  altitude=talt, model="IS",
                                  fov_dir="front")
        if(coords != 'geo'):
            lonc, latc = coord_conv(lonc, latc, "geo", coords,
                                    altitude=t_c_alt,
                                    date_time=date_time)

        # Save into output arrays
        lat_center[ig] = latc
        lon_center[ig] = lonc

    return lat_center, lon_center

# test code
def main():

    # input parameters
    #rad_list = ["bks", "wal", "fhe", "fhw", "cve", "cvw"]
    rad_list = ["bks"]
    #rad_list = ["wal"]
    ftype = "fitacf"
    bmnum = 7

    stm = dt.datetime(2012,1,1)
    etm = dt.datetime(2012,2,29)
    #season = "winter"
    season = None
    baseLocation = "../data/sqlite3/"

    coords = "geo"
    ftype = "fitacf"

    for rad in rad_list:
        obj = latc_lonc_to_db(rad, bmnum, stm, etm, coords="geo",
                              ftype=ftype, season=season, baseLocation=baseLocation)
    return obj
if __name__ == "__main__":
    obj = main()
    obj.add_latclonc_to_db()
