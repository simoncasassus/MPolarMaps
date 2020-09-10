from astropy.io import fits 
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import matplotlib.colors as colors

# python map_overlay.py

# -----------------------------------------------------------
# nice fonts
# -----------------------------------------------------------
matplotlib.rc('font', family='sans-serif') 
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 12}

matplotlib.rc('font', **font)



        
plt.figure(figsize=(5, 5)) #(4,6)


#ax = plt.subplot(211)
ax = plt.subplot(111)

ax.tick_params(axis='both',length = 5, width=1., color = 'grey',direction='in',left=True, right=True,bottom=True, top=True)

ax.tick_params(axis='both',length = 5, width=0.5, color = 'black')

plt.setp(ax.get_xticklabels(),visible=True)#, fontsize=6)

######################################################################

file_b4='/Users/simon/common/ppdisks/HD34282_multifreq/polarmaps/polarmaps_b4/profile_cyl.dat'
file_b7='/Users/simon/common/ppdisks/HD34282_multifreq/polarmaps/polarmaps_b7/profile_cyl.dat'

(azim, profile_b4) = np.loadtxt(file_b4,unpack=True)
(azim, profile_b7) = np.loadtxt(file_b7,unpack=True)

Imax_b4=np.max(profile_b4)
print( "Imax_b4",Imax_b4)
Imax_b7=np.max(profile_b7)
print( "Imax_b7",Imax_b7)



i_peak_b4=np.argmax(profile_b4)
nx_profile = len(azim)                            
#delta_pix_profile = int(np.rint(i_peak_b4 - (-1+(nx_profile+1)/2)))
delta_pix_profile = int(np.rint(i_peak_b4 - ((nx_profile+1)/2)))
profile_b4=np.roll(profile_b4,-delta_pix_profile)


profile_b7=np.roll(profile_b7,-delta_pix_profile)




azim=azim-180.

plt.ylim(-0.05,1.05)
plt.xlim(-180.,180.)
plt.plot(azim,profile_b4/np.max(profile_b4),linewidth=2.5,alpha=0.6,label=r'$145\,$GHz',color='C0')
plt.plot(azim,profile_b7/np.max(profile_b7),linewidth=2.,alpha=0.8,label=r'$345\,$GHz',color='C1')

ax.set_ylabel(r'$I(\nu,\phi)$ / $I_{\rm max}(\nu)$')
ax.set_xlabel(r'$\phi$ / deg')

plt.legend(fontsize=12)


plt.subplots_adjust(hspace=0.)

fileout='fig_azimuthal_profiles.pdf'
print( fileout)
plt.savefig(fileout, bbox_inches='tight')



