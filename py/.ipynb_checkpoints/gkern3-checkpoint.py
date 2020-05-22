import numpy as np
from scipy.stats import norm
def gkern3(sigma=1,length=3):
    hw=int((length-1)/2)
    x=np.arange(-hw,hw+1,1)
    y=np.zeros(x.shape)
    for i,item in enumerate(x):
        y[i] = norm.cdf(item+0.5,loc=0, scale=sigma) \
                - norm.cdf(item-0.5,loc=0, scale=sigma)
    
    return y/y.sum()
