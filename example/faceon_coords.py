import numpy as np


def foncoords(dalpha_star, ddelta_star, PA, inc, dx_0, dy_0):
    PArad = np.pi * PA / 180.
    incrad = np.pi * inc / 180.

    #orientation tests:
    #dalpha_star =  0.1*np.sin(PArad-np.pi/2)-dx_0
    #ddelta_star =  0.1*np.cos(PArad-np.pi/2)-dy_0
    #print("dalpha_star",dalpha_star,"ddelta_star",ddelta_star)

    da_star_ring0 = dalpha_star + dx_0
    dd_star_ring0 = ddelta_star + dy_0

    da_star_rot = da_star_ring0 * np.cos(PArad) - dd_star_ring0 * np.sin(PArad)
    dd_star_rot = da_star_ring0 * np.sin(PArad) + dd_star_ring0 * np.cos(PArad)
    da_star_fon = da_star_rot / np.fabs(np.cos(incrad))
    dd_star_fon = dd_star_rot
    print("da_star_fon", da_star_fon)
    print("dd_star_fon", dd_star_fon)
    return da_star_fon, dd_star_fon


#pm_astropy.py
dalpha_star = 0.00818441734188446
ddelta_star = -0.00423169342127494

#LB19SB16

#PA (160.4093138812665, 0.09629171227726374, 0.09640281949765495) 160.41974841260688
#inc (129.98890808056711, 0.07450643168454008, 0.0702485409296969) 129.99158057657752
#dra_off (-0.0022451509046847243, 0.00030731763895720386, 0.00034338828832711655) -0.002255584250509137
#ddec_off (-0.0031291613691573775, 0.0003880973750100287, 0.00041186849457229774) -0.0030864188866638826

PA = 160.4093138812665
inc = 129.98890808056711
# offset of ring centroid:
dx_0 = -0.0022451509046847243
dy_0 = -0.0031291613691573775

da_star_fon, dd_star_fon = foncoords(dalpha_star, ddelta_star, PA, inc, dx_0,
                                     dy_0)

