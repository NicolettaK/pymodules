import healpy as hp
import numpy as np
import camb
from camb import model, initialpower

def make_cl(r=0, return_cl=True, Alens=1., lmax=2500, tau=0.06):
    lmax = float(lmax)
    pars = camb.CAMBparams()
    pars.WantTensors = True
    pars.set_cosmology(H0=67.86682, ombh2=0.02227716, omch2=0.1184293, mnu=0.06, omk=0, tau=tau)
    pars.InitPower.set_params(ns=0.9682903, r=r)
    pars.set_for_lmax(lmax-50, lens_potential_accuracy=0)
    results = camb.get_results(pars)
    powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
    if Alens=0:
        totCL = powers['unlensed_total'].T
    else:
        tot_tesor = powers['tensor'].T
        tot_lens = powers['lensed_scalar'].T*Alen
        totCL = tot_tensor+tot_lens
    if return_cl:
        l = np.arange(len(totCL[0]))
        totCL = 2*np.pi*totCL/l/(l+1)
        totCL[:,0] = 0 
    return totCL

def make_map(r=0, nside=256, Alens=1., lmax=2500, seed=None):
    cl = make_cl(r=r, Alens=Alens, lmax=lmax)
    if seed:
        np.random.seed(seed)
    else:
        np.random.seed()
    cmb = hp.synfast(cl, nside, new=True)
    return cmb
     
