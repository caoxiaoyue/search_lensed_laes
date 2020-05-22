import subprocess
import os
import utils as ul
import copy
import pickle
import gzip

dir1 = os.getcwd()
os.chdir('/data/inspur_disk03/userdir/caoxy/eboss_lya/data/2')
dir2 = os.getcwd()

output = subprocess.Popen(['ls','-l'],stdout=subprocess.PIPE,shell=True)
lines = output.stdout.read().splitlines()
#######line.decode("utf-8")

dir3 = '/data/inspur_disk03/userdir/caoxy/eboss_lya/data/3/'
nhit = 0
res_tot = {}
tot_file_num =len(lines)
for i,line in enumerate(lines):
    #if (i < 1500) or (i>1520): continue ####fallback
    print(i,'of',tot_file_num)
    print(line)
    try:
        res_sub = ul.search_hits(line)
        nhit += res_sub['fiber_hits'].size
    except Exception as e:
        print(e)
        print('--------bug happen-----------------')
    else:
        print('number of hit',nhit)
        print('--------sucess-----------------')
        out_fname = dir3+line[0:-3]+'_hits'+'.pzip'
        if not os.path.exists(out_fname):
            pickle.dump(res_sub,gzip.open(out_fname,'wb'))
            
        if not res_tot:
            res_tot = copy.deepcopy(res_sub)
        else:
            res_tot = ul.concatenate_dict(res_tot,res_sub)
            

pickle.dump(res_tot,gzip.open(dir3+'../'+'lya_hits.pzip','wb'))
    
    