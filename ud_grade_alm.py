import healpy as hp
import numpy as np

def ud_grade_alm(map_in, nside_out):
    alm = hp.map2alm(map_in, lmax=3*(int(nside_out))-1)
    map_out = hp.alm2map(alm, int(nside_out))
    return map_out


