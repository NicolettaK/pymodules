import healpy as hp
import numpy as np
import camb
from camb import model, initialpower

def make_cl(r=0, return_cl=True, unlensed=False, lmax=2500):
    pars = camb.CAMBparams()
    pars.WantTensors = True
    pars.set_cosmology(H0=67.86682, ombh2=0.02227716, omch2=0.1184293, mnu=0.06, omk=0, tau=0.06)
    pars.InitPower.set_params(ns=0.9682903, r=r)
    pars.set_for_lmax(lmax, lens_potential_accuracy=0)
    results = camb.get_results(pars)
    powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
    if unlensed:
        totCL = powers['unlensed_total']
    else:
        totCL=powers['total']
    totCL = totCL.T
    if return_cl:
        l = np.arange(len(totCL[0]))
        totCL = 2*np.pi*totCL/l/(l+1)
        totCL[:,0:2] = np.NaN
    return totCL

def make_map(r=0, nside=256):
    cl = make_cl(r=r)
    cmb = hp.synfast(cl, nside, new=True)
    return cmb
     
