import sys
import numpy as np
from scipy.ndimage import map_coordinates
from astropy.io import fits as pf
from astropy.wcs import WCS
import os
from copy import deepcopy
include_path='/Users/simon/common/python/include/'
sys.path.append(include_path)
from ImUtils.Resamp import gridding,loadfits
import ImUtils.Cube2Im  as Cube2Im
from astropy.time import Time

filename_source='/Users/simon/common/ppdisks/HD163296/DSHARP_continuum/guvmem_runs/mem_lS0.0_lL0.0_nogrid/mod_out.fits'
filename_im='mod_out_fullim.fits'
fileout='mod_out_z.fits'

side=4.0  #arcsec

######################################################################

SIMBADCenter=False # if true uses SIMBAD coordinates. 
RAICRS='15h15m48.4459023859s'
DECICRS='-37d09m16.026315179s'
pm_ra=-19.114  # mas/yr
pm_dec=-23.140 # mas/yr

######################################################################

cube0, hdr0 = loadfits(filename_source)
Cube2Im.slice0(filename_source,filename_im)
im1, hdr1 =loadfits(filename_im)
hdr2 = deepcopy(hdr1)


if SIMBADCenter:
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    #15 15 48.4459023859 -37 09 16.026315179
    c = SkyCoord(RAICRS,DECICRS, frame='icrs')
    #proper motions mas/yr : 	-19.114 -23.140 [0.110 0.081 90]

    obsdate = hdr1['DATE'] # 
    obstime = Time(obsdate, format='isot', scale='utc')
    tobs=obstime.jd
    print("tobs",tobs)
    
    obs2000='2000-1-1T00:00:00.0'
    obsJ2000 = Time(obs2000, format='isot', scale='utc')
    tobsJ2000=obsJ2000.jd
    delta_t=(tobs-tobsJ2000)/365.
    print("tobsJ2000 ",tobsJ2000,"delta_t",delta_t,"yr")
    
    RAJ2000=c.ra.degree
    DECJ2000=c.dec.degree

    #c = SkyCoord(ra=10.625*u.degree, dec=41.2*u.degree, frame='icrs')
    #c = SkyCoord(10.625, 41.2, frame='icrs', unit='deg')
    #c.ra.degree
    #c.dec.degree

    print("RAJ2000",RAJ2000)
    print("DECJ2000",DECJ2000)
    
    RAOBS=RAJ2000+(pm_ra*delta_t/(1E3*3600.))/np.cos(DECJ2000*np.pi/180.)
    DECOBS=DECJ2000+pm_dec*delta_t/(1E3*3600.)

    offset= 3600.* np.sqrt( (DECOBS - hdr2['CRVAL2'])**2 +  np.cos(np.pi * DECOBS/ 180.) * (RAOBS - hdr2['CRVAL1'])**2)
    print("offset from CRVALi (arcsec): ",offset)

    print("RAOBS",RAOBS, " vs CRVAL1 ",hdr2['CRVAL1'])
    print("DECOBS",DECOBS, "vs CRVAL2 ",hdr2['CRVAL2'] )

    hdr2['CRVAL1']=RAOBS
    hdr2['CRVAL2']=DECOBS

nx=np.rint(side/(3600.*hdr1['CDELT2']))
ny=np.rint(side/(3600.*hdr1['CDELT2']))

if ( (nx % 2) == 0):
    nx=nx+1
    ny=ny+1


print( "nx ",nx,"ny ", ny)

#hdr2['CRVAL1']=alpha_star
#hdr2['CRVAL2']=delta_star
hdr2['NAXIS1']=nx
hdr2['NAXIS2']=ny
hdr2['CRPIX1']=(nx+1)/2
hdr2['CRPIX2']=(ny+1)/2

print( "gridding  ")

import re
if re.search(r"pix",hdr1['BUNIT'],re.IGNORECASE):
    unitscale=(hdr2['CDELT2']/hdr1['CDELT2'])**2
else:
    unitscale=1.

print( "unitscale ", unitscale)

resamp=gridding(filename_im,hdr2)

resamp=resamp*unitscale


print( "fileout ", fileout)
hdrout=deepcopy(hdr2)
hdrout['CRVAL3']=hdr0['CRVAL3']

pf.writeto(fileout,resamp, hdrout, overwrite=True)


