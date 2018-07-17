import healpy as hp
import numpy as np

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

