import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
import h5py
import copy
from scipy.signal import find_peaks
from scipy.special import erf

def gkern3(sigma=1,length=3):
    hw=int((length-1)/2)
    x=np.arange(-hw,hw+1,1)
    y=np.zeros(x.shape)
    for i,item in enumerate(x):
        y[i] = norm.cdf(item+0.5,loc=0, scale=sigma) \
                - norm.cdf(item-0.5,loc=0, scale=sigma)
    return y/y.sum()

def find_blocks(value, a):
    """
    a: the searching array
    value: the value
    searching value a in array a;
    """
    # Create an array that is 1 where a is `value`, and pad each end with an extra 0.
    isvalue = np.concatenate(([0], np.equal(a, value).view(np.int8), [0]))
    absdiff = np.abs(np.diff(isvalue))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

def filter_spec(kernel,flux,invvar):
    con1 = np.correlate(flux*invvar,kernel,mode='same')
    con1[0] = 0.0
    con1[-1] = 0.0
    #print('kern-----',kernel)
    #print('flux-----',flux)
    #print('invvar---',invvar)
    #print('con1-------',con1)
    con2 = np.correlate(invvar,kernel**2,mode='same')
    #print('con2-------',con1)
    con2[0] = 0.0
    con2[-1] = 0.0

    lflux = con1 / np.clip(con2,1e-10,None)
    sn = con1 / np.clip( np.sqrt(np.clip(con2,1e-10,None)),1e-10,None ) 
    noise = 1.0 / np.clip( np.sqrt(np.clip(con2,1e-10,None)),1e-10,None ) 
    return lflux,sn,noise


def findhits_lya(resflux=None, resivar=None, zans=None, wave=None,lsigma=2.2,lksize=15):
    sn_min = 6.0 #minimum SNR for candidate Lyman-alpha
    maxhits = 5  #maximum number of hits in a spectrum
    
    #Fake identification
    lya_lam = 1215.67
    o2_lam = np.sqrt(3727.092 * 3729.875)
    hb_lam = 4862.683
    o3_lam = 5008.239
    ha_lam = 6564.614
    fake_list = np.array([hb_lam, o3_lam, ha_lam])
    
    fake_lam = lya_lam * fake_list / o2_lam
    
    if resflux.ndim ==1:
        npix = resflux.shape[0]
        nspec = 1
    else:
        nspec,npix = resflux.shape
        
    wmask = np.logical_and(wave > 3700.0, 
                           wave < 4800.0)
    
    mask5577 = np.logical_or(wave < 5575.0, 
                              wave > 5584.0)
    
    #lsigma = 2.2  #correspond to 150km/s
    #lksize = 15  #kernel size
    lkern = gkern3(lsigma, lksize)  #no offset term compare to idl version
    #print(lkern)
    
    #To convert to per-pixel baseline
    dlam = wave * 1.0e-4 * np.log(10.0)
    #To store the output line-flux noise levels
    out_noise = np.zeros(resflux.shape)
    
    fakesn = np.copy(0.0 * fake_lam).reshape(1,fake_lam.size) #fake s/n
    peak_id = np.array([-1],dtype='int')
    spec_id = np.array([-1],dtype='int')
    peak_wave = np.array([-1.0],dtype=np.float)
    peak_sn = np.array([0.0],dtype=np.float)
    hw = 15
    peak_sn_array = np.zeros(2*hw+1,dtype='float').reshape(1,2*hw+1)
    
    #print('-------------',nspec)
    
    for i in range(nspec):
        thisflux = np.copy(resflux[i,:] * dlam)
        thisinvvar = np.copy(mask5577 * resivar[i,:] / dlam**2)
            
        fil_flux, sn,thisnoise = filter_spec(kernel=lkern, flux=thisflux, invvar=thisinvvar)
        
        out_noise[i,:] = thisnoise[:]
        #Find things that qualify as "hits":
        htest = ((sn * wmask) > sn_min).astype('int')
        bi = find_blocks(1,htest)
        
        #fig,ax = plt.subplots(2,2)
        #print(thisflux.shape)
        #ax[0,0].plot(thisflux)
        #ax[0,1].plot(thisinvvar)
        #ax[1,0].plot(sn)
        #ax[1,1].plot(htest)
        #fig.savefig('test.png')
        #print('fil_flux max',fil_flux.max())
        #print('S/N max',sn.max())
        
        
        #save hits info
        nhits = bi.shape[0]
        #print('number of hits',nhits)
        if ((nhits > 0) and (nhits < maxhits)):
            #print('------have hits------')
            sn_this = np.zeros(nhits)  #s/n of hits
            sn_this_id = np.zeros(nhits,dtype='int')  #id of s/n hits
            sn_this_array = np.zeros((nhits,2*hw+1),dtype='float')
            for j,ids in enumerate(bi):
                snsub = sn[ids[0]:ids[1]+1]
                maxid = np.argmax(snsub)
                sn_this_id[j] = ids[0]+maxid
                sn_this[j] = snsub[maxid]
                sn_this_array[j,:] = sn[sn_this_id[j]-hw:sn_this_id[j]+hw+1]
                
            # Compute fake-wavelength SNRs:
            fakethis = np.zeros((nhits,fake_list.size))
            #print('faketthis shape',fakethis.shape)
            for jj in range(nhits):
                #print(type(sn_this_id[jj]))
                obsfakelam = wave[sn_this_id[jj]] / lya_lam * fake_lam
                #print('interpo----',np.interp(obsfakelam, wave, sn))
                fakethis[jj,:] = np.interp(obsfakelam, wave, sn) * (obsfakelam > wave.min()) * (obsfakelam < wave.max())
            
            fakesn = np.concatenate((fakesn,fakethis),axis=0)
            peak_sn_array = np.concatenate((peak_sn_array, sn_this_array),axis=0)
            
            peak_id = np.concatenate((peak_id,sn_this_id))
            #print('peak id',peak_id,'sn_this_id',sn_this_id)
            spec_id = np.concatenate((spec_id,i+np.zeros(nhits,dtype='int')))
            peak_wave = np.concatenate((peak_wave,wave[sn_this_id]))
            peak_sn = np.concatenate((peak_sn,sn_this))

    peak_sn_array = peak_sn_array[1:,:]
    fakesn = fakesn[1:,:]
    peak_id = peak_id[1:]
    spec_id = spec_id[1:]
    peak_wave = peak_wave[1:]
    peak_sn = peak_sn[1:]
        
    plate_hits = np.array([zans[i][0] for i in spec_id],dtype='int')
    mjd_hits =  np.array([zans[i][2] for i in spec_id],dtype=np.int)
    fiber_hits = np.array([zans[i][3] for i in spec_id],dtype='int')
    #redshift of lens, no_qso z
    z_hits = np.array([zans[i][-9] for i in spec_id],dtype='float') #or 12
        
    ret_dict = {'fakesn':fakesn,'peak_id':peak_id,'spec_id':spec_id,
               'peak_wave':peak_wave,'peak_sn':peak_sn, #'peak_sn_array':peak_sn_array,
               'plate_hits':plate_hits,'mjd_hits':mjd_hits,
               'fiber_hits':fiber_hits,'z_hits':z_hits}
        
    return ret_dict

def search_hits(h5_file=None):
    with h5py.File(h5_file, 'r') as f:
        arr0 = f['flux']
        flux = arr0[:]
        
        arr1 = f['loglam']
        loglam = arr1[:]
        
        arr2 = f['zans']
        zans = arr2[:]
        
        arr3 = f['newsynflux']
        synflux = arr3[:]
        
        arr4 = f['res_ivar']
        res_ivar = arr4[:]
    res = findhits_lya(resflux=flux-synflux,resivar=res_ivar,zans=zans,wave=10**loglam)
    return res
        
def concatenate_dict(d,d_sub):
    for key in d:
        if isinstance(d[key],np.ndarray):
            d[key] = np.concatenate((d[key],d_sub[key]))
        if isinstance(d[key],list):
            d[key] = d[key] + d_sub[key]
    return d
            
def cut_dict(dicts,ids):
    d = copy.deepcopy(dicts)
    for key in d:
        d[key] = d[key][ids]
    return d

def smooth(array,kernal_size):
    """
    simulate the idl smooth function
    """
    kernal = np.ones(kernal_size)/kernal_size
    res = np.correlate(array,kernal,mode='same')
    res[0] = array[0]
    res[-1] = array[-1]
    return res

def convert_vdisp_to_dpix(v):
    '''
    units of v is in km/s
    '''
    c=3*1e5
    d_lnlam = v/c
    d_log10lam = d_lnlam/np.log(10)
    dpix = d_log10lam/1e-4
    return dpix

def convert_dpix_to_vdisp(dpix):
    '''
    units of v is in km/s
    '''
    d_log10lam = dpix*1e-4
    d_lnlam = d_log10lam*np.log(10)
    v = d_lnlam*3*1e5
    return v

#def return_mask_region_id(res_ivar):
#    id_mask = []
#    mask = np.zeros(res_ivar.shape)
#    for i in range(mask.shape[0]):
#        mask[i,:] = np.equal(res_ivar[i,:],0)
#        id_mask.append(find_blocks(1,mask[i,:] )) #the begin and end id of mask 
#    return id_mask

def return_mask_region_id(res_ivar):
    mask = np.equal(res_ivar,0)
    id_mask = find_blocks(1,mask) #the begin and end id of mask 
    return id_mask


def filter_sn(resflux, resivar, wave, vel_sigma=None, pix_sigma=2.2, kern_size=15):
    if vel_sigma is not None:
        sigma = convert_vdisp_to_dpix(vel_sigma)
    elif pix_sigma is not None:
        sigma = pix_sigma
    else:
        print('you must set the sigma of kernel')
    #print('sigma value=',sigma)
    
    kernel = gkern3(sigma=sigma,length=kern_size)
    #print(kernel)
    #To convert to per-pixel baseline
    dlam = wave * 1.0e-4 * np.log(10.0)
 
    thisflux = np.copy(resflux * dlam)
    thisinvvar = np.copy(resivar / dlam**2)
    fil_flux, sn, thisnoise = filter_spec(kernel=kernel, flux=thisflux, invvar=thisinvvar)
    return fil_flux, sn, thisnoise


def return_peak(wave, sn, height=2.5, distance=2.0):
    indexes, _ = find_peaks(sn, height=height, distance=distance)
    return wave[indexes],sn[indexes]


def skew_gauss(npix , A, x, omega, s, hw_sub=None):
    '''
    lam is wavelength variable, which is dependant variable in fitting
    '''
    if hw_sub is None:
        lam = np.arange(0,npix,1)
    else:
        lam = np.arange(0,npix)[:,None]+np.arange(-hw_sub,hw_sub+1,1)*1.0/hw_sub/2
        
    factor = np.sqrt(np.pi)/2.0
    gauss_comp = np.exp(-0.5*np.power((lam-x)/omega,2)) 
    erf_comp = factor * ( erf(s*(lam-x)/omega)+1.0 )
    
    if hw_sub is None:
        return  A*gauss_comp*erf_comp
    else:
        return np.mean(A*gauss_comp *erf_comp,axis=1)
    
    
def double_peak(npix , f_ratio, offset, omega1, s1, A2, x2, omega2, s2, hw_sub=None):
    if hw_sub is None:
        lam = np.arange(0,npix,1)
    else:
        lam = np.arange(0,npix)[:,None]+np.arange(-hw_sub,hw_sub+1,1)*1.0/hw_sub/2
    
    factor = np.sqrt(np.pi)/2.0
    
    x1 = x2 + offset
    A1 = f_ratio*A2
    gauss_comp1 = np.exp(-0.5*np.power((lam-x1)/omega1,2)) 
    erf_comp1 = factor * ( erf(s1*(lam-x1)/omega1)+1.0 )
    subfunc1 = A1*gauss_comp1*erf_comp1

    gauss_comp2 = np.exp(-0.5*np.power((lam-x2)/omega2,2)) 
    erf_comp2 = factor * ( erf(s2*(lam-x2)/omega2)+1.0 )
    subfunc2 = A2*gauss_comp2*erf_comp2
    
    subfunc = subfunc1 + subfunc2
    if hw_sub is None:
        return  subfunc
    else:
        return np.mean(subfunc,axis=1)
    
def loglike_single(cube,ndim, nparams):
    mock = skew_gauss(npix , cube[0], cube[1], cube[2], cube[3],2)
    return -0.5*np.sum( (obs-mock)**2*weight)

def prior_single(cube,ndim, nparams):
    cube[0] = 100*cube[0] #norm.ppf(cube[0],loc=10,scale=1) #norm.ppf(cube[0],loc=0.3,scale=0.1) #normalization
    cube[1] = norm.ppf(cube[1],loc=wave_fit[np.argmax(obs)],scale=10) #norm.ppf(cube[1],loc=3050,scale=5)  #center x
    cube[2] = 20*cube[2]#norm.ppf(cube[2],loc=10,scale=3) #omega,scatter
    cube[3] = 20.0*cube[3] - 10.0 #norm.ppf(cube[3],loc=1.8,scale=0.3) #skewness [-20,20]
    
def loglike_double(cube,ndim, nparams):
    mock = double_peak(npix , cube[0], cube[1], cube[2], cube[3], cube[4], cube[5], cube[6], cube[7],hw_sub=2)
    return -0.5*np.sum( (obs-mock)**2*weight)

def prior_double(cube,ndim, nparams):
    cube[0] = 1.0*cube[0]  #flux ratio
    cube[1] = 30.0*cube[1] - 30.0 #norm.ppf(cube[1],loc=wave_fit[np.argmax(obs)],scale=15) #offset
    cube[2] = 20.0*cube[2] #norm.ppf(cube[2],loc=10,scale=3) #omega,scatter
    cube[3] = 20.0*cube[3] - 10.0 #norm.ppf(cube[3],loc=1.8,scale=0.3) #skewness [-20,20]
    cube[4] = 100.0*cube[4] #flux 2
    cube[5] = norm.ppf(cube[5],loc=wave_fit[np.argmax(obs)],scale=15) #center of red peak
    cube[6] = 20.0*cube[6] #omega2
    cube[7] = 20.0*cube[7] - 10.0 #skewness2
    
    
from astropy import constants as const
from astropy.cosmology import Planck15 as cosmo
def cal_einstein_radius_in_arcsec(z_l,z_s,sigma):
    #b_sie = 4*pi * sigma^2/c^2 * D_ls/D_s
    c = const.c.value / 1000 #to km/s 
    v_fac = (sigma/c)**2
    d_fac = cosmo.angular_diameter_distance_z1z2(z_l, z_s).value/cosmo.angular_diameter_distance(z_s).value
    radian_to_arcsec = 180.0/np.pi*60*60
    return 4.0*np.pi*v_fac*d_fac*radian_to_arcsec