import healpy as hp
import numpy as np 

def Kcmb_Krj(nu):
    nu = np.array(nu, dtype=np.float64)
    x1 = nu/56.78
    S1 = x1**2.*np.exp(x1)/(np.exp(x1)-1)**2.
    return S1

def scaling_dust(freq1, freq2, sp_index=1.59, d_sp_index=0, T=19.6, spectrum=False, rj=False):
    freq1 = float(freq1)
    freq2 = np.array(freq2, dtype=np.float64)
    x1 = freq1/56.78
    x2 = freq2/56.78
    S1 = x1**2.*np.exp(x1)/(np.exp(x1)-1)**2.
    S2 = x2**2.*np.exp(x2)/(np.exp(x2)-1)**2.
    vd = 375.06/18.*T
    A = (np.exp(freq1/vd)-1)/(np.exp(freq2/vd)-1)
    scaling_factor = A*(freq2/freq1)**(sp_index+1)
    err_scaling_factor = np.sqrt((A*(freq2/freq1)**(sp_index+1)*np.log(freq2/freq1)*d_sp_index)**2.)
    if rj is False:
        scaling_factor = S1/S2*scaling_factor
        err_scaling_factor = S1/S2*err_scaling_factor
    if d_sp_index==0 and spectrum is False:
        return scaling_factor
    elif d_sp_index==0 and spectrum is True:
        return scaling_factor**2.
    elif d_sp_index!=0 and spectrum is False:
        return scaling_factor, err_scaling_factor
    elif d_sp_index!=0 and spectrum is True:
        return scaling_factor**2., 2*scaling_factor*err_scaling_factor

def scaling_synch(freq1, freq2, sp_index=-3, d_sp_index=0, spectrum=False, rj=False):
    freq1 = float(freq1)
    freq2 = np.array(freq2, dtype=np.float64)
    x1=freq1/56.78
    x2=freq2/56.78
    S1=x1**2.*np.exp(x1)/(np.exp(x1)-1)**2.
    S2=x2**2.*np.exp(x2)/(np.exp(x2)-1)**2.
    scaling_factor = (freq2/freq1)**(sp_index)
    err_scaling_factor = np.sqrt(((freq2/freq1)**(sp_index)*np.log(freq2/freq1)*d_sp_index)**2.)
    if rj is False:
        scaling_factor = S1/S2*scaling_factor
        err_scaling_factor = S1/S2*err_scaling_factor
    if d_sp_index==0 and spectrum is False:
        return scaling_factor
    elif d_sp_index==0 and spectrum is True:
        return scaling_factor**2.
    elif d_sp_index!=0 and spectrum is False:
        return scaling_factor, err_scaling_factor
    elif d_sp_index!=0 and spectrum is True:
        return scaling_factor**2., 2*scaling_factor*err_scaling_factor

def ud_grade_alm(map_in, nside_out):
    alm = hp.map2alm(map_in, lmax=3*(int(nside_out))-1)
    map_out = hp.alm2map(alm, int(nside_out))
    return map_out

def remove_mono(hmap):
    map_nomono = hmap-np.mean(hmap)
    return map_nomono

def gauss_map(cl, nside=256):
    print 'gauss_map', nside
    beta_sim_gauss = hp.synfast(cl, nside)
    return beta_sim_gauss

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

def apply_filter(hmap, filt, nside=256):
    import healpy as hp

    print 'apply_filter', nside
    lmax = 3*nside
    filt = filt[0:lmax]
    alm = hp.map2alm(hmap, lmax=lmax)
    almf = np.array([hp.almxfl(alm[i], filt) for i in xrange(len(alm))])
    hmap_out = hp.alm2map(almf, nside)
    return hmap_out
