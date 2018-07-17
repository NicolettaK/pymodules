import healpy as hp
import numpy as np

def high_pass_filter(l1, l2, lmax):
    ell = np.arange(lmax)+1
    wl = np.zeros(len(ell))
    wl[np.where(ell>=l2)] = 1.
    for l in range(l1, l2):
        wl[l] = 0.5*(1-np.cos(np.pi*(l-l1)/(l2-l1)))
    return wl
    
def low_pass_filter(l1, l2, lmax):
    ell = np.arange(lmax+1)
    wl = np.zeros(len(ell))
    wl[np.where(ell<=l1)] = 1.
    for l in range(l1, l2):
        wl[l] = 0.5*(1-np.cos(np.pi*(l2-l)/(l2-l1)))
    return wl

def apply_filter(hmap, filter, nside=256):
    print 'apply_filter', nside
    lmax = len(filter)
    alm = hp.map2alm(hmap, lmax=lmax)
    almf = hp.almxfl(alm, filter)
    hmap_out = hp.alm2map(almf, nside)
    return hmap_out
