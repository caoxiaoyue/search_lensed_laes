import numpy as np
def filter_spec(kernel,flux,invvar):
    con1 = np.correlate(flux*invvar,kernel,mode='same')
    con1[0] = 0
    con1[-1] = 0
    
    con2 = np.correlate(invvar,kernel**2,mode='same')
    con2[0] = 0
    con2[-1] = 0

    lflux = con1 / np.clip(con2,1e-10,None)
    sn = con1 / np.clip( np.sqrt(np.clip(con2,1e-10,None)),1e-10,None ) 
    noise = 1.0 / np.clip( np.sqrt(np.clip(con2,1e-10,None)),1e-10,None ) 
    return lflux,sn,noise
