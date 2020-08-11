import sys
import numpy as np

include_path='/Users/simon/common/python/include/polarmaps/'
sys.path.append(include_path)

import MPolarMaps



PA=76.
inc=32.4 # degrees



#######################################################################

a_max=1.



filename_source='V4046_ALMA_zoom.fits'
workdir='polarmaps_V4046_ALMA_default/'  # directory for products 



M=MPolarMaps.Setup(
    filename_source=filename_source,
    workdir=workdir,
    PA=PA,  # deg
    inc=inc,  # deg 
    RA=False, # set to False to use fits header info
    DEC=False, 
    dra_off=0., # offset from RA, DEC, in arcsec
    ddec_off=0.,
    XCheckInv=False, # debug tool
    DoRadialProfile=True, 
    DoAzimuthalProfile=True,
    PlotRadialProfile=True,
    a_min=0.27, #arcsec, ring radial range
    a_max=0.6, #arcsec, ring radial range
    zoomfactor=1., # use value > 1. to shrink the field of view. 
    y_label=r'$I_\mathrm{b7}$ / Jy beam$^{-1}$', # for the radial plots
    ForceCube2Im=False, # in case source is a datacube
    wBaseNoise=True, # for radial profile stats
    noise_radius=0.1, # use for base noise estimate
    Verbose=True, 
    wBaseNoiseCore=True) # use base noise at the origin of radial profiles. 



M.prep_files()
M.polar_expansions()


#M.workdir='polarmaps_V4046_ALMA_shifted/'  # directory for products 
#M.dra_off=0.1
#M.ddec_off=0.1
#M.prep_files()
#M.polar_expansions()




######################################################################

OptimM=MPolarMaps.OptimModel(M,
                             RunMCMC=False,
                             Nit=300, #MCMC iterations
                             nwalkers=20,
                             burn_in=100,
                             n_cores_MCMC=30)

rangePA=40.
rangeinc=40.
rangedra_off=0.1
rangeddec_off=0.1

M.workdir='optimpolarmaps_PA_inc_V4046_ALMA/'  # directory for products 

#OptimM.domain=( ('PA',(M.PA-rangePA/2.,M.PA+rangePA/2.)), ('inc',(M.inc-rangeinc/2.,M.inc+rangeinc/2.)), ('dra_off',(M.dra_off-rangedra_off/2.,M.dra_off+rangedra_off/2.)), ('ddec_off',(M.ddec_off-rangeddec_off/2.,M.ddec_off+rangeddec_off/2.)))

OptimM.domain=( ('PA',(M.PA-rangePA/2.,M.PA+rangePA/2.)), ('inc',(M.inc-rangeinc/2.,M.inc+rangeinc/2.)))

#from pprint import pprint as pp
#pp(OptimM)
#pp(M)

OptimM.ConjGrad(M)#


#OptimM.emcee(M)

from copy import deepcopy
Mgrid=deepcopy(M)

Mgrid.workdir='gridpolarmaps_V4046_ALMA/'  # directory for products 
OptimM.domain=( ('PA',(Mgrid.PA-rangePA/2.,Mgrid.PA+rangePA/2.)), ('inc',(Mgrid.inc-rangeinc/2.,Mgrid.inc+rangeinc/2.)))
OptimM.nmesh_grid=20
OptimM.SetOptim=False # set to True if you want to assign grid optimal values to M
OptimM.Grid(Mgrid)
print("returned from grid PA:",Mgrid.PA," inc: ",Mgrid.inc)


M.workdir='optimpolarmaps_offset_V4046_ALMA/'  # directory for products 
OptimM.domain=( ('dra_off',(M.dra_off-rangedra_off/2.,M.dra_off+rangedra_off/2.)), ('ddec_off',(M.ddec_off-rangeddec_off/2.,M.ddec_off+rangeddec_off/2.)))
OptimM.SetOptim=True # set to True if you want to assign  optimal values to M
OptimM.ConjGrad(M)#

#M.workdir='polarmaps_gridmin_V4046_ALMA/'  # directory for products 
#M.prep_files()
#M.polar_expansions()


from copy import deepcopy
Mgrid=deepcopy(M)

print("calling grid offset  PA:",Mgrid.PA," inc: ",Mgrid.inc," dra_off:",Mgrid.dra_off," ddec_off: ",Mgrid.ddec_off)
OptimM.domain=( ('dra_off',(Mgrid.dra_off-rangedra_off/2.,Mgrid.dra_off+rangedra_off/2.)), ('ddec_off',(Mgrid.ddec_off-rangeddec_off/2.,Mgrid.ddec_off+rangeddec_off/2.)))
Mgrid.workdir='gridpolarmaps_offset_V4046_ALMA/'  # directory for products 
OptimM.nmesh_grid=20
OptimM.SetOptim=False  # set to True if you want to assign grid optimal values to M 
OptimM.Grid(Mgrid)

print("returned from grid offset  PA:",Mgrid.PA," inc: ",Mgrid.inc," dra_off:",Mgrid.dra_off," ddec_off: ",Mgrid.ddec_off)

#M.workdir='polarmaps_offset_gridmin_V4046_ALMA/'  # directory for products 
#M.prep_files()
#M.polar_expansions()
