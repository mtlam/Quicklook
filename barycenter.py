"""
Code for calculating barycentric frequencies

Need to calculate velocity of spot on Earth!


Notes:
Consider upgrading with skyfield package?
"""
import numpy as np
from matplotlib.pyplot import *
import ephem
from astropy import units
from astropy.coordinates import EarthLocation
from astropy.time import Time
from datetime import datetime




def getLocalHour(longitude,time):


    pass




equitorial_velocity = 0.4638 #km/s
def calculateSpinVelocity(latitude,longitude,time):
    """
    Calculates the velocity vector of the observatory on the Earth due to the Earth's spin
    latitude and longitude are astropy angles, time is an astropy Time.
    Define x as the direction of the Earth's motion, y as the sunward direction, z in along the ecliptic axis
    If hour is midnight, then the observatory's velocity is in the direction of the orbital velocity.
    """
    
    pass
    
    print ( type(latitude),type(longitude),type(time) )
    print ( time ) 

    raise SystemExit

    # Convert from z=rotation axis to z=ecliptic axis
    

    
    pass







def calculateRelativeVelocity(obs_coords,psr_coords,MJD):
    """
    obs_coords in (x,y,z) from PSRFITS file
    psr_coords in (RA,DEC) from PSRFITS file (both strings)
    Returns relative velocity in km/s to first-order
    """

    t = Time(MJD, format='mjd')
    time = t.datetime.strftime("%Y/%m/%d %H:%M:%S")
    


    el = EarthLocation(x=obs_coords[0]*units.meter,y=obs_coords[1]*units.meter,z=obs_coords[2]*units.meter)#,frame='itrs')
    lon,lat,height = el.to_geodetic()

    #spin_vel = calculateSpinVelocity(lon,lat,t)

    obs = ephem.Observer()
    obs.lon = ephem.degrees(str(lon.value))
    obs.lat = ephem.degrees(str(lat.value))
    obs.date = time
    obs.pressure = 0 #no refraction

    # From pyephem docs:
    # Both hlon and hlat have a special meaning for the Sun and Moon. For a Sun body, they give the Earth's heliocentric longitude and latitude
    # Note: Sun at vernal equinox has hlong = 0, hlat = 0, Earth will return hlong = 180 degrees! This has been verified

    sun = ephem.Sun()
    sun.compute(time,epoch=2000)

    vel_ra = ephem.degrees(sun.hlong-3*np.pi/2)
    vel_dec = ephem.degrees(-1*sun.hlat)
    pulsar_ra = ephem.hours(psr_coords[0])
    pulsar_dec = ephem.degrees(psr_coords[1])

    sep = ephem.separation((vel_ra,vel_dec),(pulsar_ra,pulsar_dec))

    return 30*np.cos(sep) 



c = 2.9979e5
def calculateGamma(v): #in km/s!
    """
    Calculates Gamma, as in DM = DM_topo/Gamma (see Pennucci et al 2014)
    """
    beta = v / c
    return np.sqrt((1+beta)/(1-beta))



def convertDMtopo(DMtopo,obs_coords,psr_coords,MJD):
    v = calculateRelativeVelocity(obs_coords,psr_coords,MJD)
    Gamma = calculateGamma(v)
    return DMtopo/Gamma





if __name__=="__main__":

    #ANT_X   =            882589.65 / [m] Antenna ITRF X-coordinate (D) 
    #ANT_Y   =          -4924872.32 / [m] Antenna ITRF Y-coordinate (D)  
    #ANT_Z   =          3943729.348 / [m] Antenna ITRF Z-coordinate (D)
    #RA      = '19:09:47.448'       / Right ascension (hh:mm:ss.ssss)  
    #DEC     = '-37:44:13.920'      / Declination (-dd:mm:ss.sss) 
    #STT_IMJD=                55275 / Start MJD (UTC days) (J - long integer)    
    #STT_SMJD=                36103 / [s] Start time (sec past UTC 00h) (J) 

    ANT_X = 882589.65 #* units.meter
    ANT_Y = -4924872.32# * units.meter
    ANT_Z = 3943729.348# * units.meter
    RA = '19:09:47.448'
    DEC = '-37:44:13.920'
    STT_IMJD = 55275
    STT_SMJD = 36103

    MJD = STT_IMJD+STT_SMJD/86400.0
    MJD = 57184.5+11 #this was the test date I used approximately on the solstice

    t = Time(MJD, format='mjd')
    time = t.datetime.strftime("%Y/%m/%d %H:%M:%S")


    utc_hour = t.datetime.hour+t.datetime.minute/60.0+t.datetime.second/3600.0

    el = EarthLocation(x=ANT_X*units.meter,y=ANT_Y*units.meter,z=ANT_Z*units.meter)#,frame='itrs')
    lon,lat,height = el.to_geodetic()
    
    local_hour = utc_hour + lon.hourangle


    raise SystemExit

    #print calculateRelativeVelocity((ANT_X,ANT_Y,ANT_Z),(RA,DEC),MJD)


