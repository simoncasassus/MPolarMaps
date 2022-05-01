import sys
import numpy as np
import scipy as sp
import os
import os.path
from scipy import ndimage
from astropy.io import fits as pf
import re
from copy import deepcopy
from astropy.wcs import WCS
from pylab import *

include_path = '/Users/simon/common/python/include/'
sys.path.append(include_path)
from ImUtils.Resamp import *
import ImUtils.Cube2Im as Cube2Im

import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt

#import pyfits
#import matplotlib.pyplot as plt
#import re
#import scipy as S
#import scipy.ndimage
#from scipy import ndimage

if not sys.warnoptions:
    import os, warnings
    #warnings.simplefilter("default") # Change the filter in this process
    warnings.simplefilter("ignore")  # Change the filter in this process
    #os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses
    os.environ["PYTHONWARNINGS"] = "ignore"  # Also affect subprocesses


def cartesian2polar(outcoords, inputshape, origin, fieldscale=1.):
    """Coordinate transform for converting a polar array to Cartesian coordinates. 
    inputshape is a tuple containing the shape of the polar array. origin is a
    tuple containing the x and y indices of where the origin should be in the
    output array."""

    rindex, thetaindex = outcoords
    x0, y0 = origin

    theta = thetaindex * 2 * np.pi / (inputshape[0] - 1)
    #theta = 2. * np.pi - theta

    y = rindex * np.cos(theta) / fieldscale
    x = rindex * np.sin(theta) / fieldscale

    #print( "r",rindex,"theta",theta,"x",x,"y",y)

    ix = -x + x0
    iy = y + y0

    return (iy, ix)


def polar2cartesian(outcoords, inputshape, origin, fieldscale=1.):
    yindex, xindex = outcoords
    x0, y0 = origin
    #    theta0, r0 = origin
    nx = inputshape[0] - 1
    ny = inputshape[1] - 1
    #    x0=((nx+1)/2)-1
    #    y0=((ny+1)/2)-1
    x = -float(xindex) + x0
    y = float(yindex) - y0

    #    theta = np.arctan2(-y, -x) + np.pi

    theta = np.arctan2(x, y)
    if (theta < 0):
        theta = theta + 2. * np.pi
    #    theta = np.arctan2(-y, x) + np.pi

    thetaindex = (theta * nx / (2. * np.pi))
    #theta_index = np.round((theta + np.pi) * (inputshape[1]-1) / (2 * np.pi))

    #    thetaindex = np.arctan2(yindex-y0,xindex-x0) * nx / (2. * np.pi)
    rindex = fieldscale * np.sqrt(x**2 + y**2)

    #    itheta = thetaindex + theta0
    #    ir = rindex + r0

    #print( xindex, yindex, " x ",x," y ", y," theta ",theta)

    return (rindex, thetaindex)


def exec_prep_files(M):

    filename_source = M.filename_source
    workdir = M.workdir
    PA = M.PA
    inc = M.inc
    RA = M.RA
    DEC = M.DEC
    dra_off = M.dra_off
    ddec_off = M.ddec_off
    XCheckInv = M.XCheckInv
    DoRadialProfile = M.DoRadialProfile
    ProfileExtractRadius = M.ProfileExtractRadius
    DoAzimuthalProfile = M.DoAzimuthalProfile
    PlotRadialProfile = M.PlotRadialProfile
    a_min = M.a_min
    a_max = M.a_max
    zoomfactor = M.zoomfactor
    Grid = M.Grid
    y_label = M.y_label
    ForceCube2Im = M.ForceCube2Im
    #wBaseNoise=M.wBaseNoise
    #noise_radius=M.noise_radius
    #wBaseNoiseCore=M.wBaseNoiseCore

    cosi = np.cos(inc * np.pi / 180.)

    fieldscale = M.fieldscale  # shrink radial field of view of polar maps by this factor

    os.system("rm -rf  " + workdir)

    os.system("mkdir " + workdir)
    inbasename = os.path.basename(filename_source)
    filename_fullim = re.sub('.fits', '_fullim.fits', inbasename)
    filename_fullim = workdir + filename_fullim

    hdu0 = pf.open(filename_source)
    hdr0 = hdu0[0].header
    if ((hdr0['NAXIS'] > 2) or ForceCube2Im):
        hdu = Cube2Im.slice0(filename_source, filename_fullim)
        im1 = hdu.data
        hdr1 = hdu.header
    else:
        os.system("rsync -va " + filename_source + " " + filename_fullim)
        hdu = pf.open(filename_fullim)
        im1 = hdu[0].data
        hdr1 = hdu[0].header

    hdr1.pop('CRVAL3', None)

    if (isinstance(RA, bool)):
        if (not RA):
            RA = hdr1['CRVAL1']
            DEC = hdr1['CRVAL2']
    elif (not isinstance(RA, float)):
        sys.exit("must provide a pointing center with RA, DEC in degrees")

    M.RA = RA
    M.DEC = DEC

    nx = int(hdr1['NAXIS1'] / zoomfactor)
    ny = nx
    if ((nx % 2) == 0):
        nx = nx + 1
        ny = ny + 1

    hdr2 = deepcopy(hdr1)

    hdr2['NAXIS1'] = nx
    hdr2['NAXIS2'] = ny
    hdr2['CRPIX1'] = (nx + 1) / 2  # 1 offset
    hdr2['CRPIX2'] = (ny + 1) / 2
    hdr2['CRVAL1'] = RA
    hdr2['CRVAL2'] = DEC

    if M.Verbose:
        print("zooming on  center", RA, DEC)

    resamp = gridding(filename_fullim, hdr2, fullWCS=False)

    if M.DumpAllFitsFiles:
        fileout_centered = re.sub('fullim.fits', 'centered.fits',
                                  filename_fullim)
        pf.writeto(fileout_centered, resamp, hdr2, overwrite=True)

    hdu = pf.PrimaryHDU()
    hdu.data = resamp
    hdu.header = hdr2

    M.Hdu = hdu


def exec_polar_expansions(M):
    filename_source = M.filename_source
    workdir = M.workdir
    PA = M.PA
    inc = M.inc
    RA = M.RA
    DEC = M.DEC
    dra_off = M.dra_off
    ddec_off = M.ddec_off
    XCheckInv = M.XCheckInv
    DoRadialProfile = M.DoRadialProfile
    ProfileExtractRadius = M.ProfileExtractRadius
    DoAzimuthalProfile = M.DoAzimuthalProfile
    PlotRadialProfile = M.PlotRadialProfile
    PlotAzimuthalProfile = M.PlotAzimuthalProfile
    a_min = M.a_min
    a_max = M.a_max
    a_max_plot = M.a_max_plot
    zoomfactor = M.zoomfactor
    Grid = M.Grid
    y_label = M.y_label
    ForceCube2Im = M.ForceCube2Im
    wBaseNoise = M.wBaseNoise
    noise_radius = M.noise_radius
    wBaseNoiseCore = M.wBaseNoiseCore
    fieldscale = M.fieldscale  # shrink radial field of view of polar maps by this factor

    cosi = np.cos(inc * np.pi / 180.)

    inbasename = os.path.basename(filename_source)
    filename_fullim = re.sub('.fits', '_fullim.fits', inbasename)
    filename_fullim = workdir + filename_fullim

    #hdu=deepcopy(M.Hdu)
    hdu = M.Hdu
    resamp = hdu.data
    hdr2 = hdu.header
    hdrshifted = deepcopy(hdr2)

    nx = hdr2['NAXIS1']
    ny = hdr2['NAXIS2']

    if M.Verbose:
        print("RA =", RA)
        print("DEC =", DEC)

    if ((abs(dra_off) > 0.) or (abs(ddec_off) > 0.)):
        RA = RA + ((dra_off / 3600.) / np.cos(DEC * np.pi / 180.))
        DEC = DEC + (ddec_off / 3600.)
        if M.Verbose:
            print("about to apply shit for offset")
            print("RA =", RA)
            print("DEC =", DEC)

        hdrshifted['CRVAL1'] = RA
        hdrshifted['CRVAL2'] = DEC
        resamp = gridding(hdu, hdrshifted, fullWCS=False)

        if M.DumpAllFitsFiles:
            print("PUNCHING SHIFTED CENTER")
            fileout_centered = re.sub('fullim.fits', 'centered.fits',
                                      filename_fullim)
            pf.writeto(fileout_centered, resamp, hdr2, overwrite=True)

    if M.Verbose:
        print("running polarexpansion with PA", PA, "and inc", inc)

    rotangle = PA
    im1rot = ndimage.rotate(resamp, rotangle, reshape=False)

    hdurot = pf.PrimaryHDU()
    hdurot.data = im1rot
    hdurot.header = hdr2

    if M.DumpAllFitsFiles:
        fileout_rotated = re.sub('fullim.fits', 'rotated.fits',
                                 filename_fullim)
        hdurot.writeto(fileout_rotated, overwrite=True)

    hdr3 = deepcopy(hdr2)
    if M.UVplane:
        hdr3['CDELT1'] = hdr3['CDELT1'] / np.fabs(cosi)
    else:
        hdr3['CDELT1'] = hdr3['CDELT1'] * np.fabs(cosi)

    im3 = gridding(hdurot, hdr3)
    if M.DumpAllFitsFiles:
        if M.UVplane:
            fileout_stretched = re.sub('fullim.fits', 'squeezed.fits',
                                       filename_fullim)
        else:
            fileout_stretched = re.sub('fullim.fits', 'stretched.fits',
                                       filename_fullim)

        pf.writeto(fileout_stretched, im3, hdr2, overwrite=True)

    im_polar = sp.ndimage.geometric_transform(im3,
                                              cartesian2polar,
                                              order=1,
                                              output_shape=(im3.shape[0],
                                                            im3.shape[1]),
                                              extra_keywords={
                                                  'inputshape':
                                                  im3.shape,
                                                  'fieldscale':
                                                  fieldscale,
                                                  'origin':
                                                  (((nx + 1) / 2) - 1,
                                                   ((ny + 1) / 2) - 1)
                                              })

    nphis, nrs = im_polar.shape

    hdupolar = pf.PrimaryHDU()
    hdupolar.data = im_polar
    hdrpolar = hdupolar.header
    hdrpolar['CRPIX1'] = 1
    hdrpolar['CRVAL1'] = 0.
    hdrpolar['CDELT1'] = 2. * np.pi / nphis
    hdrpolar['CRPIX2'] = 1
    hdrpolar['CRVAL2'] = 0.
    hdrpolar['CDELT2'] = (hdr3['CDELT2'] / fieldscale)
    hdupolar.header = hdrpolar

    if M.DumpAllFitsFiles:

        fileout_polar = re.sub('fullim.fits', 'polar.fits', filename_fullim)
        hdupolar.writeto(fileout_polar, overwrite=True)

    ######################################################################
    # profiles

    if (DoRadialProfile):

        Iprof = np.average(im_polar, axis=1)
        sIprof = np.std(im_polar, axis=1)
        rrs = 3600. * (np.arange(hdrpolar['NAXIS2']) - hdrpolar['CRPIX2'] +
                       1) * hdrpolar['CDELT2'] + hdrpolar['CRVAL2']
        #np.savetxt('test.out', (rrs,Iprof))

        bmaj = hdr2['CDELT2']
        if ('BMAJ' in hdr2):
            bmaj = hdr2['BMAJ']
        if M.Verbose:
            print(("bmaj = ", bmaj, "\n"))
        Nind = 2. * np.pi * rrs * np.fabs(cosi) / (bmaj * 3600.)  #cosi *
        iNind1 = np.argmin(np.fabs(Nind - 1.0))
        rNind1 = rrs[iNind1]
        if M.Verbose:
            print(("radius for Nind=1 ", rNind1))
        Nind[np.where(Nind < 1.0)] = 1.0
        dispIprof = sIprof.copy()

        if (noise_radius <= 0.):
            noise_radius = rNind1

        if (wBaseNoise):
            noise_basal = dispIprof[np.argmin(np.fabs(rrs - noise_radius))]
            if (wBaseNoiseCore):
                #dispIprof[np.where( (rrs < noise_radius) & (dispIprof<noise_basal))] = noise_basal
                dispIprof[np.where(rrs < noise_radius)] = noise_basal
            else:
                dispIprof[np.where(dispIprof < noise_basal)] = noise_basal
            if M.Verbose:
                print((">using noise floor of ", noise_basal,
                       "  from dispersion at radius ", noise_radius))
        else:
            if M.Verbose:
                print((">not applying noise floor"))

        sIprof = dispIprof / np.sqrt(Nind)

        if M.DumpAllFitsFiles:

            save_prof = np.zeros((hdrpolar['NAXIS2'], 4))
            #print( save_prof.shape)
            save_prof[:, 0] = rrs
            save_prof[:, 1] = Iprof
            save_prof[:, 2] = sIprof
            save_prof[:, 3] = dispIprof
            fileout_radialprofile = re.sub('fullim.fits', 'radial_profile.dat',
                                           filename_fullim)
            np.savetxt(fileout_radialprofile,
                       save_prof)  # x,y,z equal sized 1D arrays

        ia_min = 0
        ia_max = len(rrs) - 1

        if M.Verbose:
            print(("ia_min ", ia_min, " ia_max", ia_max))

        if (a_min > 0):
            ia_min = np.argmin(np.abs(rrs - a_min))
            if M.Verbose:
                print(("a_min ", a_min, " ia_min", ia_min))
        if (a_max > 0):
            ia_max = np.argmin(np.abs(rrs - a_max))
            if M.Verbose:
                print(("a_max ", a_max, " ia_max", ia_max))

        if (M.noise_floor < 0.):
            M.noise_floor = np.min(dispIprof[ia_min:ia_max])

        #peakval=np.max(Iprof)
        #varIprof=(dispIprof/peakval)**2

        varIprof = (dispIprof / M.noise_floor)**2

        chi2 = np.sum(varIprof[ia_min:ia_max])
        if M.Verbose:
            print("returning from funcs_PolarMaps with chi2=", chi2)
        if ((Grid) and (M.MinimizeRmsProfile)):
            return chi2

        ######################################################################
        # Polar average

        Iprof_b = Iprof.reshape(Iprof.shape[0], 1)
        #im_polar_av = im_polar.copy
        im_polar_av = np.copy(im_polar)

        im_polar_av[:, :] = Iprof_b

        if M.DumpAllFitsFiles:
            fileout_polar_av = re.sub('fullim.fits', 'polar_av.fits',
                                      filename_fullim)
            hdupolar.data = im_polar_av
            hdupolar.writeto(fileout_polar_av, overwrite=True)

        ######################################################################
        # AZIM IM

        imazim = sp.ndimage.geometric_transform(
            im_polar_av,
            polar2cartesian,
            order=0,
            output_shape=(im_polar.shape[0], im_polar.shape[1]),
            extra_keywords={
                'inputshape': im_polar.shape,
                'fieldscale': fieldscale,
                'origin': (((nx + 1) / 2) - 1, ((ny + 1) / 2) - 1)
            })

        hdu_stretch_av = pf.PrimaryHDU()
        hdu_stretch_av.data = imazim
        hdu_stretch_av.header = hdr2

        if M.DumpAllFitsFiles:
            fileout_stretched_av = re.sub('fullim.fits', 'stretched_av.fits',
                                          filename_fullim)
            hdu_stretch_av.writeto(fileout_stretched_av, overwrite=True)

        ######################################################################
        # back to sky - project

        hdr3 = deepcopy(hdr2)
        hdr3['CDELT1'] = hdr3['CDELT1'] / np.fabs(cosi)

        im4 = gridding(hdu_stretch_av, hdr3)
        if M.DumpAllFitsFiles:
            fileout_proj = re.sub('fullim.fits', 'azim_av_proj.fits',
                                  filename_fullim)
            pf.writeto(fileout_proj, im4, hdr2, overwrite=True)

        ######################################################################
        # back to sky - rotate

        im4drot = ndimage.rotate(im4, -rotangle, reshape=False)

        if M.DumpAllFitsFiles:
            fileout_drotated = re.sub('fullim.fits', 'azim_av_drot.fits',
                                      filename_fullim)
            pf.writeto(fileout_drotated, im4drot, hdr2, overwrite=True)

        ######################################################################
        # azimuthal profiles

    if ((DoAzimuthalProfile) or (ProfileExtractRadius > 0)):

        im_polar_4profiles = np.copy(im_polar)

        rs = 3600. * (
            (np.arange(hdrpolar['NAXIS2']) - hdrpolar['CRPIX2'] + 1) *
            hdrpolar['CDELT2'] + hdrpolar['CRVAL2'])
        if (a_min > 0):
            ia_min = np.argmin(np.abs(rs - a_min))
            if M.Verbose:
                print(("a_min", a_min, " ia_min", ia_min))
            im_polar_4profiles[0:ia_min, :] = 0

        if (a_max > 0):
            ia_max = np.argmin(np.abs(rs - a_max))
            if M.Verbose:
                print(("a_max", a_max, " ia_max", ia_max))
            im_polar_4profiles[ia_max:, :] = 0

        #im_polar_4profiles[0:20,:] = 0

        #import matplotlib.pyplot as plt
        #plt.imshow(im_polar_4profiles)
        #plt.show()

        Imax_profile = np.amax(im_polar_4profiles, axis=0)
        phis = (180. / np.pi) * (
            (np.arange(hdrpolar['NAXIS1']) - hdrpolar['CRPIX1'] + 1) *
            hdrpolar['CDELT1'] + hdrpolar['CRVAL1'])

        (nphis, nrs) = im_polar_4profiles.shape

        #print( "nphis ",nphis,"len phis",len(phis))
        #print( "nrs ",nrs)
        #iphis=np.arange(1,nphis+1)
        #irs=np.arange(1,nrs+1)
        #iiphis,iirs  = np.meshgrid(iphis, irs)
        #rrs =  3600.*(iirs-hdrpolar['CRPIX2']+1)*hdrpolar['CDELT2']+hdrpolar['CRVAL2']
        #pphis =  (180./np.pi)*(iiphis-hdrpolar['CRPIX1']+1)*hdrpolar['CDELT1']+hdrpolar['CRVAL1']

        ivec_rmax = np.argmax(im_polar_4profiles, axis=0)
        Rmax_profile = 3600. * (
            (ivec_rmax - hdrpolar['CRPIX2'] + 1) * hdrpolar['CDELT2'] +
            hdrpolar['CRVAL2'])
        if M.Verbose:
            print(("Rmax shape:", Rmax_profile.shape))

        WholeRing = False
        if (ProfileExtractRadius > 0):
            ringrad = ProfileExtractRadius
        elif (WholeRing):
            #Iprof = np.average(im_polar_4profiles,axis=1)
            #ringrad=np.asscalar(Iprof[np.argmax(Iprof)])
            ringrad = np.sum(
                Rmax_profile * Imax_profile) / np.sum(Imax_profile)
            if M.Verbose:
                print(("Rcav ", ringrad))
        else:
            iphimax = np.argmax(Imax_profile)
            if M.Verbose:
                print(("iphimax", iphimax))
            iringrad = ivec_rmax[iphimax]
            ringrad = Rmax_profile[iphimax]

        if M.Verbose:
            print(("iringrad ", iringrad, " ringrad ", ringrad))

        Icyl_profile = im_polar_4profiles[iringrad, :]

        #np.savetxt(workdir+'rmax_profile.dat', (rs,Rmax_profile))
        #np.savetxt(workdir+'profile_peak.dat', (rs,Imax_profile))
        #np.savetxt(workdir+'profile_cyl.dat', (rs,Icyl_profile))

        if M.DumpAllFitsFiles:
            save_prof = np.zeros((hdrpolar['NAXIS2'], 2))
            save_prof[:, 0] = phis
            save_prof[:, 1] = Rmax_profile
            np.savetxt(workdir + 'rmax_profile.dat', save_prof)

            save_prof = np.zeros((hdrpolar['NAXIS2'], 2))
            save_prof[:, 0] = phis
            save_prof[:, 1] = Imax_profile
            np.savetxt(workdir + 'profile_peak.dat', save_prof)

            save_prof = np.zeros((hdrpolar['NAXIS2'], 2))
            save_prof[:, 0] = phis
            save_prof[:, 1] = Icyl_profile
            np.savetxt(workdir + 'profile_cyl.dat', save_prof)

        avring = np.sum(Rmax_profile * Imax_profile) / np.sum(Imax_profile)
        #rmsring=np.std(Rmax_profile*Imax_profile/np.sum(Imax_profile))
        rmsring2 = np.sqrt(
            np.sum((Rmax_profile - avring)**2 * Imax_profile /
                   np.sum(Imax_profile)))
        #print( "rmsring ",rmsring,"rmsring2 ", rmsring2 )

        if ((Grid) and (M.MinimizeRmsRing)):
            return rmsring2

    ######################################################################
    # CROSS CHECK INVERSE TRANSFORM

    if (XCheckInv):
        im_x = sp.ndimage.geometric_transform(im_polar,
                                              polar2cartesian,
                                              order=0,
                                              output_shape=(im_polar.shape[0],
                                                            im_polar.shape[1]),
                                              extra_keywords={
                                                  'inputshape':
                                                  im_polar.shape,
                                                  'fieldscale':
                                                  fieldscale,
                                                  'origin':
                                                  (((nx + 1) / 2) - 1,
                                                   ((ny + 1) / 2) - 1)
                                              })

        fileout_stretched_x = re.sub('fullim.fits', 'stretched_x.fits',
                                     filename_fullim)
        pf.writeto(fileout_stretched_x, im_x, hdr2, overwrite=True)

        hdr3 = deepcopy(hdr2)
        hdr3['CDELT1'] = hdr3['CDELT1'] / np.fabs(cosi)
        fileout_proj_x = re.sub('fullim.fits', 'x_proj.fits', filename_fullim)
        im4_x = gridding(fileout_stretched_x, hdr3)
        pf.writeto(fileout_proj_x, im4_x, hdr2, overwrite=True)

        im4_x_drot = ndimage.rotate(im4_x, -rotangle, reshape=False)
        fileout_drotated_x = re.sub('fullim.fits', 'x_drot.fits',
                                    filename_fullim)
        pf.writeto(fileout_drotated_x, im4_x_drot, hdr2, overwrite=True)

        fileout_diff_x = re.sub('fullim.fits', 'x_diff.fits', filename_fullim)
        pf.writeto(fileout_diff_x, resamp - im4_x_drot, hdr2, overwrite=True)

    ######################################################################
    # profile plotting

    if (PlotAzimuthalProfile):

        # -----------------------------------------------------------
        # nice fonts
        # -----------------------------------------------------------
        matplotlib.rc('font', family='sans-serif')
        matplotlib.rcParams.update({'font.size': 12})
        #matplotlib.rcParams.update({'figsize': (7,5)})

        plt.figure(figsize=(7, 5))
        axprofile = plt.subplot(111)

        print("loading ", workdir + 'rmax_profile.dat')
        (phis, Rmax_profile) = np.loadtxt(workdir + 'rmax_profile.dat',
                                          unpack=True)

        #plt.setp(axprofile.get_xticklabels(),visible=True) #, fontsize=6)
        #plt.setp(axprofile.get_yticklabels(),visible=True) #, fontsize=6)

        plt.xlim(0., 360.)
        plt.ylim(0.9 * np.min(Rmax_profile), 1.1 * np.max(Rmax_profile))

        plt.plot(phis, Rmax_profile, linewidth=0.1, linestyle='solid')

        #plt.fill_between(rrs, Iprof+sIprof, Iprof-sIprof, lw=0.1,color='r', alpha=0.3, interpolate=True, step='mid')
        #plt.fill_between(rrs, Iprof+dispIprof, Iprof-dispIprof, lw=0.1,color='b', alpha=0.3, interpolate=True)

        if (y_label == ''):
            plt.ylabel(r'$r_\mathrm{max}$ / arcsec')
        else:
            plt.ylabel(y_label)

        plt.xlabel(r'$\phi$ / rads')

        fileout_fig = re.sub('fullim.fits', 'fig_azprofile.pdf',
                             filename_fullim)
        plt.savefig(fileout_fig, bbox_inches='tight')

    if (PlotRadialProfile):

        # -----------------------------------------------------------
        # nice fonts
        # -----------------------------------------------------------
        matplotlib.rc('font', family='sans-serif')
        matplotlib.rcParams.update({'font.size': 12})
        #matplotlib.rcParams.update({'figsize': (7,5)})

        plt.figure(figsize=(7, 5))
        axprofile = plt.subplot(111)

        rmax = np.max(rrs)

        plt.setp(axprofile.get_xticklabels(), visible=True)  #, fontsize=6)
        plt.setp(axprofile.get_yticklabels(), visible=True)  #, fontsize=6)
        if (a_max_plot > 0.):
            plt.xlim(0., a_max_plot)

        #plt.ylim(np.min(Iprof),1.1*np.max(Iprof))
        plt.ylim(-0.1 * np.max(Iprof), 1.1 * np.max(Iprof))

        plt.plot(rrs,
                 rrs * 0.,
                 color='black',
                 linewidth=0.1,
                 linestyle='solid')

        #plt.plot(rrs,Icutav,color='blue',linewidth=1,linestyle='solid')
        plt.plot(rrs, Iprof, color='grey', linewidth=0.1, linestyle='solid')

        plt.fill_between(rrs,
                         Iprof + sIprof,
                         Iprof - sIprof,
                         lw=0.1,
                         color='r',
                         alpha=0.3,
                         interpolate=True,
                         step='mid')

        plt.fill_between(rrs,
                         Iprof + dispIprof,
                         Iprof - dispIprof,
                         lw=0.1,
                         color='b',
                         alpha=0.3,
                         interpolate=True)

        if (y_label == ''):
            plt.ylabel(r'$\langle I(r) \rangle$ / mJy pix$^{-1}$')
        else:
            plt.ylabel(y_label)

        plt.xlabel(r'$r$ / arcsec')

        fileout_fig = re.sub('fullim.fits', 'fig_profile.pdf', filename_fullim)
        plt.savefig(fileout_fig, bbox_inches='tight')
