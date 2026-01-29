import sys
import numpy as np
from copy import deepcopy

include_path='/home/simon/common/python/include/'
sys.path.append(include_path)

import MPolarMaps.Master


PA=160.2981635916787
inc=50.046591099746855


sourcedir = '/home/simon/PDS70snow/slab/proc/output_imoptim_fixq_fixamax_regultau_regulT_regulSigma/'

filename_source = 'imlogTdust.fits'

M=MPolarMaps.Master.Setup(
    filename_source=sourcedir+filename_source,
    filename_source_errormap=sourcedir+'serrimlogTdust.fits',
    workdir='',
    PA=PA,  # deg
    inc=inc,  # deg 
    dra_off=-0.002490618551486642, # offset from RA, DEC, in arcsec
    ddec_off=-0.013013964142577486,
    DoRadialProfile=True, 
    DoAzimuthalProfile=True,
    PlotRadialProfile=True,
    a_min=0.3, #arcsec, ring radial range
    a_max=1.5, #arcsec, ring radial range
    zoomfactor=1., # use value > 1. to shrink the field of view. 
    y_label=r'$I_\mathrm{b7}$ / Jy beam$^{-1}$', # for the radial plots
    ForceCube2Im=False, # in case source is a datacube
    wBaseNoise=False, # if True sets dispersion in radial profile to fixed value at noise_radius
    noise_radius=0.1, # use for base noise estimate
    Verbose=True,
    noise_floor=1E-2,  # for Chi2 stats 
    wBaseNoiseCore=True) # use base noise at the origin of radial profiles. 

#######################################################################
# default expansion

M.workdir='polarmaps_bestL_T/'  # directory for products
M.filename_source = sourcedir+'imlogTdust.fits'
M.filename_source_errormap=sourcedir+'serrimlogTdust.fits'
M.noise_floor=1E-2 # minimum noise value in error map 
M.prep_files()
M.polar_expansions()


M.workdir='polarmaps_bestL_Sigmag/'  # directory for products
M.filename_source = sourcedir+'imlogSigma_g.fits'
M.filename_source_errormap=sourcedir+'serrimlogSigma_g.fits'
M.noise_floor=1E-2 # minimum noise value in error map 
M.prep_files()
M.polar_expansions()

