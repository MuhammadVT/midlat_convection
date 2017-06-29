from davitpy.utils import plotUtils
from davitpy.pydarn.plotting.mapOverlay import overlayRadar, overlayFov
import matplotlib.pyplot as plt
import datetime as dt


dateTime = dt.datetime(2012, 1, 1)
fig = plt.figure(figsize=(8,8)) 
m = plotUtils.mapObj(lat_0=60., lon_0=-90, width=111e3*90, height=111e3*80,
                     coords='geo', dateTime=dateTime)
codes = ['bks', 'wal','fhe','fhw','cve','cvw']
# Plotting some radars
overlayRadar(m, fontSize=30, codes=codes, dateTime=dateTime)
# Plot radar fov
overlayFov(m, codes=codes, maxGate=70, dateTime=dateTime, fovColor="orange", fovAlpha=0.7)
#rcParams.update({'font.size': 12})
#plt.show()
fig.savefig("./plots/fov/fov.png", dpi=300)
