import numpy as np 

def Kcmb_Krj(nu):
    nu = np.array(nu, dtype=np.float64)
    x1 = nu/56.78
    S1 = x1**2.*np.exp(x1)/(np.exp(x1)-1)**2.
    return S1
