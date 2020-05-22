import subprocess
import os
import utils as ul
import copy
import pickle
import gzip
from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import medfilt
import h5py
import matplotlib.gridspec as gridspec

in_dir = '/data/inspur_disk03/userdir/caoxy/eboss_lya/data/'
lya_hits = pickle.load(gzip.open(in_dir+'lya_hits_cut.pzip'))

dir1 = os.getcwd()
os.chdir('/data/inspur_disk03/userdir/caoxy/eboss_lya/data/2')
dir2 = os.getcwd()

plate_list = lya_hits['plate_hits']
fiber_list = lya_hits['fiber_hits']
mjd_list = lya_hits['mjd_hits']
spec_id_list = lya_hits['spec_id']
peak_id_list = lya_hits['peak_id']
peak_wave_list = lya_hits['peak_wave']
peak_sn_list = lya_hits['peak_sn']


for i in range(plate_list.size):
    this_plate = plate_list[i]
    this_mjd = mjd_list[i]
    this_fiber = fiber_list[i]
    this_spec_id = spec_id_list[i]
    this_peak_id = peak_id_list[i]
    this_peak_wave = peak_wave_list[i]
    this_peak_sn = peak_sn_list[i]
    h5_file = '{}_{}.h5'.format(this_plate,this_mjd)
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
    
    
    this_wave = 10**loglam
    this_flux = flux[this_spec_id,:]
    this_synflux = synflux[this_spec_id,:]
    this_resflux = this_flux - this_synflux
    this_res_ivar = res_ivar[this_spec_id,:]
    this_sn = np.sqrt(this_resflux**2 * this_res_ivar)
    
    
    hw_sub = 25 #empirically choice
    wave_sub = this_wave[this_peak_id-hw_sub:this_peak_id+hw_sub+1]
    resflux_sub =  this_resflux[this_peak_id-hw_sub:this_peak_id+hw_sub+1]
    res_ivar_sub = this_res_ivar[this_peak_id-hw_sub:this_peak_id+hw_sub+1]
    sn_sub = this_sn[this_peak_id-hw_sub:this_peak_id+hw_sub+1]
    
    #make sure reading process is correct
    #print(zans[this_spec_id])
    #print(this_plate,this_mjd,this_fiber)
    
    #visulization
    name = r'ID: {}-{}-{}'.format(this_plate,this_mjd,this_fiber)
    hit_wave = r'Hit $\lambda$: {:.2f}'.format(this_peak_wave)
    hit_sn = r'Hit SN: {:.2f}'.format(this_peak_sn)
    textstr = '\n'.join((name,hit_wave,hit_sn))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    
    #fig, ax = plt.subplots(2,2,figsize=(20,15)) #,sharey=True
    fig = plt.figure(figsize=(30,20),tight_layout=True)
    gs = gridspec.GridSpec(3, 2)
    
    #--------------------------subplot 0,:
    ax = fig.add_subplot(gs[0, :])
    ax.plot(this_wave,this_flux,color='deepskyblue')
    ax.plot(this_wave, this_synflux,color='red')
    ax.axvline(x=this_peak_wave,color='magenta')
    id_mask = ul.return_mask_region_id(this_res_ivar)
    for this_mask_id in id_mask:
        low = this_mask_id[0]
        high = this_mask_id[1]
        #print(this_wave[low],this_wave[high-1])
        ax.axvspan(this_wave[low], this_wave[high-1], facecolor='g', alpha=0.5)
    #ax[0,0].set_xlim(3700,4800)
    up_lim = this_flux.max()*1.2
    if this_flux.min() < 0:
        low_lim = this_flux.min()*1.2
    else:
        low_lim = -0.1
    ax.set_ylim(low_lim,up_lim)
    
    #--------------------------subplot 1,0
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(this_wave,this_flux,color='deepskyblue')
    ax.plot(this_wave, this_synflux,color='red')
    ax.axvline(x=this_peak_wave,color='magenta')
    id_mask = ul.return_mask_region_id(this_res_ivar)
    for this_mask_id in id_mask:
        low = this_mask_id[0]
        high = this_mask_id[1]
        #print(this_wave[low],this_wave[high-1])
        ax.axvspan(this_wave[low], this_wave[high-1], facecolor='g', alpha=0.5)
    ax.set_xlim(3700,4800)
    id_tmp = np.where(np.logical_and(this_wave<4800,this_wave>3700))
    up_lim = this_flux[id_tmp].max()*1.2
    if this_flux[id_tmp].min() < 0:
        low_lim = this_flux[id_tmp].min()*1.2
    else:
        low_lim = -0.1
    ax.set_ylim(low_lim,up_lim)
    
    
    #--------------------------subplot 1,1
    ax = fig.add_subplot(gs[1, 1])
    ax.bar(wave_sub,sn_sub,alpha=0.5)
    ax.plot(wave_sub,sn_sub,color='red')
    ax.axvline(x=this_peak_wave,color='magenta')
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    wave_, sn_ = ul.return_peak(wave_sub,sn_sub)
    ax.plot(wave_,sn_,'*',color='black')
    
    
    #--------------------------subplot 2,0
    ax = fig.add_subplot(gs[2, 0])
    _,sn_filter_0,_ = ul.filter_sn(resflux=resflux_sub,resivar=res_ivar_sub,wave=wave_sub,pix_sigma=2.2,kern_size=15)
    ax.bar(wave_sub,sn_filter_0,alpha=0.5)
    ax.plot(wave_sub,sn_filter_0,color='red')
    ax.axvline(x=this_peak_wave,color='magenta')
    wave_, sn_ = ul.return_peak(wave_sub,sn_filter_0)
    ax.plot(wave_, sn_,'*',color='black')
    
    
    #--------------------------subplot 2,1
    ax = fig.add_subplot(gs[2, 1])
    _,sn_filter_1,_ = ul.filter_sn(resflux=resflux_sub,resivar=res_ivar_sub,wave=wave_sub,pix_sigma=1.1,kern_size=15)
    ax.bar(wave_sub,sn_filter_1,alpha=0.5)
    ax.plot(wave_sub,sn_filter_1,color='red')
    ax.axvline(x=this_peak_wave,color='magenta')
    wave_, sn_ = ul.return_peak(wave_sub,sn_filter_1)
    ax.plot(wave_, sn_,'*',color='black')
    
    
    #plt.show()
    png_name = '{}-{}-{}'.format(this_plate,this_mjd,this_fiber)
    out_png = '/data/inspur_disk03/userdir/caoxy/eboss_lya/data/4/{}.png'.format(png_name)
    fig.savefig(out_png,dpi=150,bbox_inches='tight')