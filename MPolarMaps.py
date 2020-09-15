import sys
import numpy as np
import os
import os.path
from scipy import ndimage
from astropy.io import fits as pf
import re
include_path='/Users/simon/common/python/include/'
sys.path.append(include_path)
from ImUtils.Resamp import *
from copy import deepcopy
from astropy.wcs import WCS
from pylab import *

#import pyfits
#import matplotlib.pyplot as plt
#import re
#import scipy as S
#import scipy.ndimage
#from scipy import ndimage

from funcs_MPolarMaps import * 
from funcs_OptimPolarMaps import * 


class Setup():
    def __init__(self,
                 filename_source='',
                 workdir='',
                 PA=0.,  # degrees
                 inc=0., # degrees
                 RA=False, # if False read from CRVAL1
                 DEC=False,
                 dra_off=0.,
                 ddec_off=0.,
                 fieldscale=2.,
                 XCheckInv=False,
                 DoRadialProfile=True,
                 ProfileExtractRadius=-1,
                 DoAzimuthalProfile=False,
                 PlotRadialProfile=True,
                 PlotAzimuthalProfile=True,
                 a_min=-1,
                 a_max=-1,
                 a_max_plot=-1,
                 zoomfactor=1.,
                 MinimizeRmsRing=False, # for optimizations and gridding, return chi2 of azimuthal ring radius 
                 MinimizeRmsProfile=True, # for optimizations and gridding, return chi2 of radial profile dispersion
                 Grid=False,
                 y_label='',
                 ForceCube2Im=False,
                 noise_floor=-1.,
                 wBaseNoise=False,
                 noise_radius=0.,
                 wBaseNoiseCore=False,
                 Verbose=True,
                 VerboseInit=True,
                 Hdu=False,
                 DumpAllFitsFiles=True):

        
        initlocals=locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            if VerboseInit:
                print( "MPolarMaps setting ",a_attribute," to ",initlocals[a_attribute])
            setattr(self,a_attribute,initlocals[a_attribute])
        

    def prep_files(self):
        exec_prep_files(self)

    def polar_expansions(self):
        return exec_polar_expansions(self)
            


class OptimModel():
    def __init__(self,M,
                 RunConjGrad=True,
                 RunMCMC=False,
                 Nit=100, #MCMC iterations
                 nwalkers=15,
                 burn_in=50,
                 n_cores_MCMC=2,
                 TriangleFile='cornerplot.png',
                 PrintOptimStatus=True,
                 nmesh_grid=100,
                 SetOptim=True, # set model with optimal solution in the case of Gridding - for conjgrad this is always true.  
                 domain=()):  

        initlocals=locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            print( "setting ",a_attribute," to ",initlocals[a_attribute])
            setattr(self,a_attribute,initlocals[a_attribute])





    def ConjGrad(self,M):
        from funcs_OptimPolarMaps import exec_ConjGrad
        return exec_ConjGrad(M,self)


    def emcee(self,M):
        from funcs_OptimPolarMaps import exec_emcee
        result_ml=np.load(M.workdir+'result_ml.dat.npy')
        
        retvals = exec_emcee(M,result_ml,True,self)
        return retvals

    def Grid(self,M):
        from funcs_OptimPolarMaps import exec_Grid
        exec_Grid(M,self)
        return 


    
#class OptimModel():
#    
#    #def __init__(self,M,PrintOptimStatus=True): #,DoConjGrad=False, RunMCMC=False
#    #    self.PrintOptimStatus=M.PrintOptimStatus
#    #    #self.DoConjGrad=DoConjGrad
#    #    #self.RunMCMC=RunMCMC
#
#    def __init__(self, 
#                 ConjGrad=True,
#                 RunMCMC=False,
#                 Nit=1, #MCMC iterations
#                 nwalkers=1,
#                 burn_in=50,
#                 n_cores_MCMC=1,
#                 PrintOptimStatus=True,
#                 domain=()):  
#
#        initlocals=locals()
#        initlocals.pop('self')
#        for a_attribute in initlocals.keys():
#            print( "setting ",a_attribute," to ",initlocals[a_attribute])
#            setattr(self,a_attribute,initlocals[a_attribute])
#
#
#        # print( "opening log:",M.workdir+M.filelog)
#        # fout=open(M.workdir+M.filelog,"w+")
#        # M.fout=fout
#
#
#
#
#    def ConjGrad(self,M):
#        from funcs_OptimPolarMaps import exec_ConjGrad
#        M.Verbose=False
#        M.DumpAllFitsFiles=False
#        return exec_ConjGrad(self,M)
#
#
#
#
