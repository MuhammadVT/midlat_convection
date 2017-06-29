import sqlite3
import numpy as np

rad_list = ["bks", "wal", "fhe", "fhw", "cve", "cvw"]

#season = "winter"
season = "winter"
ftype = "fitacf"
baseLocation="../data/sqlite3/" + season + "/"

dbName = "ten_min_median_" + ftype + ".sqlite"
conn = sqlite3.connect(baseLocation + "ten_min_median_" + ftype + ".sqlite")
cur = conn.cursor()

rad_dict = {}
for rad in rad_list:
    # select dates
    #cur.execute("SELECT strftime('%Y-%m-%d', datetime) FROM {tb}".format(tb=rad))
    cur.execute("SELECT COUNT(*) FROM {tb} GROUP BY glatc, glonc, gazmc, datetime ORDER BY COUNT(*)".format(tb=rad))
    rows = cur.fetchall()

#    # extract the unique dates
#    rows = np.unique(rows)

    rad_dict[rad] = rows

# close db connection
conn.close()


