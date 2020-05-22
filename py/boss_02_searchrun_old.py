import subprocess
import os
import utils as ul
import copy

dir1 = os.getcwd()
os.chdir('/data/inspur_disk03/userdir/caoxy/eboss_lya/data/2')
dir2 = os.getcwd()

output = subprocess.Popen(['ls','-l'],stdout=subprocess.PIPE,shell=True)
lines = output.stdout.read().splitlines()


res_tot = {}
tot_file_num =len(lines)
for i,line in enumerate(lines):
    res_sub = ul.search_hits(line)
    if not res_tot:
        res_tot = copy.deepcopy(res_sub)
    else:
        res_tot = ul.concatenate_dict(res_tot,res_sub)
    if i%30==0:
        print(i,'of',tot_file_num)
        print('number of hits: ',res_tot['plate_hits'].size)
        print('-------------------------')

        
os.chdir('/data/inspur_disk03/userdir/caoxy/eboss_lya/data/')
with h5py.File('hits_02.h5', 'w') as f:
    dset = f.create_dataset("fakesn", data=res_tot['fakesn'])
    dset = f.create_dataset("fiber_hits", data=res_tot['fiber_hits'])
    dset = f.create_dataset("mjd_hits", data=res_tot['mjd_hits'])
    dset = f.create_dataset("peak_id", data=res_tot['peak_id'])
    
    dset = f.create_dataset("peak_sn", data=res_tot['peak_sn'])
    dset = f.create_dataset("peak_wave", data=res_tot['peak_wave'])
    dset = f.create_dataset("plate_hits", data=res_tot['plate_hits'])
    dset = f.create_dataset("spec_id", data=res_tot['spec_id'])
    
    dset = f.create_dataset("z_hits", data=res_tot['z_hits'])

    