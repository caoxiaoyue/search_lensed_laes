import numpy as np
from block_id import find_blocks
from filter_spec import filter_spec
from gkern3 import gkern3

def findhits_lya(resflux=None, resivar=None, zans=None, wave=None):
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
    
    mask5577 = np.logical_and(wave < 5575.0, 
                              wave > 5584.0)
    
    lsigma = 2.2  #correspond to 150km/s
    lksize = 15  #kernel size
    lkern = gkern3(lsigma, lksize)  #no offset term compare to idl version
    
    #To convert to per-pixel baseline
    dlam = wave * 1.0e-4 * np.log(10.0)
    #To store the output line-flux noise levels
    out_noise = np.zeros(resflux.shape)
    
    fakesn = np.copy(0.0 * fake_lam) #fake s/n
    peak_id = -1
    spec_id = -1
    peak_wave = -1.0
    peak_sn = 0.0
    
    for i in range(nspec):
        thisflux = np.copy(resflux[i,:] * dlam)
        thisinvvar = np.copy(mask5577 * resivar[i,:] / dlam**2)
        fil_flux, sn,thisnoise = filter_spec(kernel=lkern, flux=thisflux, invvar=thisinvvar)
        out_noise[i,:] = thisnoise[:]
        #Find things that qualify as "hits":
        htest = ((sn * wmask) > snmin).astype('int')
        bi = find_blocks(1,htest)
        
        #save hits info
        nhits = bi.shape[0]
        if ((nhits > 0) and (nhits < maxhits)):
            sn_this = np.zeros(nhits)  #s/n of hits
            sn_this_id = np.zeros(nhits)  #id of s/n hits
            for j,ids in enumerate(bi):
                snsub = sn[ids]
                maxid = np.argmax(snsub)
                sn_this_id[j] = ids[0]+maxid
                sn_this[j] = snsub[maxid]
                
            # Compute fake-wavelength SNRs:
            fakethis = np.zeros((nhits,fake_list.size))
            for jj in range(nhits):
                obsfakelam = wave[sn_this_id[jj]] / lya_lam * fake_lam
                fakethis[jj,:] = np.interp(obsfakelam, wave, sn) * \
                                (obsfakelam > wave.min()) * (obsfakelam < wave.max())
            
            fakesn = np.concatenate((fakesn,fakethis),axis=0)
            peak_id = np.concatenate((peak_id,sn_this_id))
            spec_id = np.concatenate((spec_id,i+np.zeros(nhits,dtype='int')))
            peak_wave = np.concatenate((peak_wave,wave[sn_this_id]))
            peak_sn = np.concatenate((peak_sn,sn_this))
        
    fakesn = fakesn[1:]
    peak_id = peak_id[1:]
    spec_id = spec_id[1:]
    peak_wave = peak_wave[1:]
    peak_sn = peak_sn[1:]
        
    plate_hits = np.array([zans[i][0] for i in spec_id],dtype='np.int8')
    mjd_hits =  np.array([zans[i][2] for i in spec_id],dtype='np.int8')
    fiber_hits = np.array([zans[i][3] for i in spec_id],dtype='np.int8')
    #redshift of lens, no_qso z
    z_hits = np.array([zans[i][-9] for i in spec_id],dtype='float') 
    print('-----------------------------')
    print('ssssssnnnn',sn)
    ret_dict = {'fakesn':fakesn,'peak_id':peak_id,'spec_id':spec_id,
               'peak_wave':peak_wave,'peak_sn':peak_sn,
               'plate_hits':plate_hits,'mjd_hits':mjd_hits,
               'fiber_hits':fiber_hits,'z_hits':z_hits}
        
    return ret_dict
        
        
        
        
        
        
            
            
            
                
            
            
            
                
        
        
        
    
    
    


