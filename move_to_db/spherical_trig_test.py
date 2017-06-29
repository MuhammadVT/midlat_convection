import numpy as np
from geomag import geomag
from datetime import date
import datetime as dt 

rad_loc_dict = {"bks" : (37.10, -77.95), "wal" : (37.93, -75.47), "fhe" : (38.86, -99.39),
            "fhw" : (38.86, -99.39), "cve" : (43.27, -120.36), "cvw" : (43.27, -120.36)}
tm = dt.datetime(2011, 1, 1).date() 
rad = "bks"
bmazm = -33.93
#bmazm = 33.93
latc = 54.73; lonc = 270.36
rad_lat, rad_lon = rad_loc_dict[rad]
rad_lon = rad_lon % 360
B = np.deg2rad(bmazm)
b_prime = np.deg2rad(90. - latc)
#sin_B = np.sin(B)
a_prime = np.deg2rad(90. - rad_lat)
#sin_A = np.sin(a_prime) * sin_B / np.sin(b_prime)
#A = np.arcsin(sin_A)
AB_dellat = np.deg2rad(np.abs(latc-rad_lat))
AB_dellon = np.deg2rad(np.abs(lonc-rad_lon))
#c_prime = 2*np.arcsin(np.sqrt((np.sin(AB_dellat/2.))**2 + np.cos(np.deg2rad(latc)) *\
#          np.cos(np.deg2rad(rad_lat)) * (np.sin(AB_dellon/2.))**2))
#c_prime = np.arccos(round(np.sin(np.deg2rad(rad_lat)) * np.sin(np.deg2rad(latc)) +\
#          np.cos(np.deg2rad(rad_lat)) * np.cos(np.deg2rad(latc)) * np.cos(AB_dellon), 8))
c_prime = np.arccos(np.sin(np.deg2rad(rad_lat)) * np.sin(np.deg2rad(latc)) +\
          np.cos(np.deg2rad(rad_lat)) * np.cos(np.deg2rad(latc)) * np.cos(AB_dellon))
s_prime = 1./2 * (a_prime + b_prime + c_prime)
#A = 2 * np.arcsin(round(np.sqrt((np.sin(s_prime - b_prime) * np.sin(s_prime - c_prime)) /\
#              (np.sin(b_prime) * np.sin(c_prime))), 6))
A = 2 * np.arcsin(np.sqrt((np.sin(s_prime - b_prime) * np.sin(s_prime - c_prime)) /\
              (np.sin(b_prime) * np.sin(c_prime))))

#sin_C = np.sin(c_prime) * sin_B / np.sin(b_prime)
#cos_B = np.cos(B)
#cos_C = 
#cos_A = -cos_B * cos_C + sin_B * sin_C * np.cos(a_prime)
#AA = np.arccos(cos_A)
#if round(AA,2) == round(A,2):
#    losvel_dir = 180 - np.rad2deg(A)
#else: 
#    losvel_dir = np.rad2deg(A)

losvel_dir = np.sign(bmazm) * (180 - np.rad2deg(A))

# convert from geo to mag by add the magnetic declanation angle to the los vel angle in geo
gm = geomag.GeoMag()
mg = gm.GeoMag(latc, lonc, h=300., time=tm)
#azm_mag = (round(bmazm - mg.dec,2)) % 360
azm_mag = (round(losvel_dir - mg.dec,2)) % 360


print "A = ", np.rad2deg(A)
print "losvel_dir = ", losvel_dir
print "azm_mag = ", azm_mag
