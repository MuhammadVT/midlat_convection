import datetime as dt
lonc = [-75, -105, -120, 15]
lonc = [x%360 for x in lonc]
# convert utc to local time in degrees
lonc_ltm = []
date_time = dt.datetime(2016, 9, 27, 6)
for lonc_i in lonc:
    lonc_tmp = lonc_i if lonc_i<=180 else lonc_i-360
    # convert utc to local time
    local_dt = date_time + dt.timedelta(hours= lonc_tmp/15.)
    ltm = local_dt.time()  
    # convert local time to degrees. e.g. 0 (or 360) degree is midnight, 180 degrees is noon time. 
    lonc_ltm.append((ltm.hour + ltm.minute/60. + ltm.second/3600.) * 15.)
lonc = lonc_ltm

