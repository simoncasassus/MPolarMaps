import sys
import numpy as np
from copy import deepcopy

include_path='/Users/simon/common/python/include/polarmaps/'
sys.path.append(include_path)

import MPolarMaps

filename_source='mod_out_z.fits'
noise_floor=1E-6 

PA=312.379492
inc=2.330638*180./np.pi  # degrees

M=MPolarMaps.Setup(
    filename_source=filename_source,
    workdir='',
    PA=PA,  # deg
    inc=inc,  # deg 
    RA=False, # set to False to use fits header info
    DEC=False, 
    dra_off=0., # offset from RA, DEC, in arcsec
    ddec_off=0.,
    DoRadialProfile=True, 
    DoAzimuthalProfile=True,
    PlotRadialProfile=True,
    a_min=0.4, #arcsec, ring radial range
    a_max=1.1, #arcsec, ring radial range
    zoomfactor=1., # use value > 1. to shrink the field of view. 
    y_label=r'$I_\mathrm{b6}$ / Jy pix$^{-1}$', # for the radial plots
    ForceCube2Im=False, # in case source is a datacube
    wBaseNoise=True, # for radial profile stats
    noise_radius=0.1, # use for base noise estimate
    Verbose=True,
    noise_floor=noise_floor,
    wBaseNoiseCore=True) # use base noise at the origin of radial profiles. 

#######################################################################
# default expansion

M.workdir='polarmaps_modout_default/'  # directory for products 
M.prep_files()
M.polar_expansions()

sys.exit()

#######################################################################
# a deliberatedly offset expansion
#M.workdir='polarmaps_shifted/'  # directory for products 
#M.dra_off=0.1
#M.ddec_off=0.1
#M.prep_files()
#M.polar_expansions()


######################################################################
# optimizations
######################################################################

OptimM=MPolarMaps.OptimModel(M,
                             RunMCMC=False,
                             Nit=500, #MCMC iterations
                             nwalkers=30,
                             burn_in=200,
                             n_cores_MCMC=30)
rangePA=40.
rangeinc=40.
rangedra_off=0.05
rangeddec_off=0.05

######################################################################
# PA, inc only

M.workdir='optimpolarmaps_PA_inc_modout/'  # directory for products 
OptimM.domain=( ('PA',(M.PA-rangePA/2.,M.PA+rangePA/2.)), ('inc',(M.inc-rangeinc/2.,M.inc+rangeinc/2.)))

#from pprint import pprint as pp
#pp(OptimM)
#pp(M)

OptimM.ConjGrad(M) # stores optimal params in M

######################################################################
# grid: PA, inc 

Mgrid=deepcopy(M)

Mgrid.workdir='gridpolarmaps_modout/'  # directory for products 
OptimM.domain=( ('PA',(Mgrid.PA-rangePA/2.,Mgrid.PA+rangePA/2.)), ('inc',(Mgrid.inc-rangeinc/2.,Mgrid.inc+rangeinc/2.)))
OptimM.nmesh_grid=20
OptimM.SetOptim=False # set to True if you want to assign grid optimal values to M
#OptimM.Grid(Mgrid)


######################################################################
# offsets only, for optimal PA, inc 

M.workdir='optimpolarmaps_offset_modout/'  # directory for products 
OptimM.domain=( ('dra_off',(M.dra_off-rangedra_off/2.,M.dra_off+rangedra_off/2.)), ('ddec_off',(M.ddec_off-rangeddec_off/2.,M.ddec_off+rangeddec_off/2.)))
OptimM.SetOptim=True # set to True if you want to assign  optimal values to M
OptimM.ConjGrad(M)#

######################################################################
# grid: offsets

Mgrid=deepcopy(M)

print("calling grid offset  PA:",Mgrid.PA," inc: ",Mgrid.inc," dra_off:",Mgrid.dra_off," ddec_off: ",Mgrid.ddec_off)
OptimM.domain=( ('dra_off',(Mgrid.dra_off-rangedra_off/2.,Mgrid.dra_off+rangedra_off/2.)), ('ddec_off',(Mgrid.ddec_off-rangeddec_off/2.,Mgrid.ddec_off+rangeddec_off/2.)))
Mgrid.workdir='gridpolarmaps_offset_modout/'  # directory for products 
OptimM.nmesh_grid=20
OptimM.SetOptim=False  # set to True if you want to assign grid optimal values to M 
#OptimM.Grid(Mgrid)


#M.workdir='polarmaps_offset_gridmin/'  # directory for products 
#M.prep_files()
#M.polar_expansions()


######################################################################
#full optimization

OptimM.domain=( ('PA',(M.PA-rangePA/2.,M.PA+rangePA/2.)), ('inc',(M.inc-rangeinc/2.,M.inc+rangeinc/2.)), ('dra_off',(M.dra_off-rangedra_off/2.,M.dra_off+rangedra_off/2.)), ('ddec_off',(M.ddec_off-rangeddec_off/2.,M.ddec_off+rangeddec_off/2.)))

M.workdir='optimpolarmaps_full_modout/'  # directory for products 
OptimM.SetOptim=True # set to True if you want to assign  optimal values to M
OptimM.ConjGrad(M)#

print("EMCEE START")
OptimM.emcee(M)

