out_dir = '/data/inspur_disk03/userdir/caoxy/eboss_lya/data/cross_talk/'
hits_list_file = './data/lya_hits_cut.csv'

readcol, hits_list_file, plate,mjd,fiber,format='(I4, L5, I4)', stringskip='#', DELIMITER=','

for i=0L, n_elements(plate)-1 do begin
    readspec, plate[i], fiber[i], mjd=mjd[i], flux=flux_0, wave=wave_0,zans=zans_0
    readspec, plate[i], fiber[i]-1, mjd=mjd[i], flux=flux_1, wave=wave_1,zans=zans_1
    readspec, plate[i], fiber[i]+1, mjd=mjd[i], flux=flux_2, wave=wave_2,zans=zans_2

    ;write hdf5
    out_file = STRCOMPRESS(out_dir + string(plate[i])+'_'+string(mjd[i])+'_'+string(fiber[i])+'_0.h5',/remove_all)
    mg_h5_putdata, out_file, 'flux', flux_0
    mg_h5_putdata, out_file, 'wave', wave_0
    mg_h5_putdata, out_file, 'zans', zans_0
    print,out_file+' hdf5 file written finished!'

    out_file = STRCOMPRESS(out_dir + string(plate[i])+'_'+string(mjd[i])+'_'+string(fiber[i])+'_1.h5',/remove_all)
    mg_h5_putdata, out_file, 'flux', flux_1
    mg_h5_putdata, out_file, 'wave', wave_1
    mg_h5_putdata, out_file, 'zans', zans_1
    print,out_file+' hdf5 file written finished!'
    
    out_file = STRCOMPRESS(out_dir + string(plate[i])+'_'+string(mjd[i])+'_'+string(fiber[i])+'_2.h5',/remove_all)
    mg_h5_putdata, out_file, 'flux', flux_2
    mg_h5_putdata, out_file, 'wave', wave_2
    mg_h5_putdata, out_file, 'zans', zans_2
    print,out_file+' hdf5 file written finished!'
endfor

end