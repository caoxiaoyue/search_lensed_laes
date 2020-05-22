;pro boss_02_searchrun, dosky=dosky

; Incorporates calls to both the single-line and the many-line
; line finding routines, to save time on the use of the preparatory
; calculations.
;
; Updated for DR9 = v5_4_45 in Feb 2012 for HST Cycle 20,
; based upon earlier version for HST Cycle 18.
;
; Updated for DR10 = v5_5+12 in Jan 2013 for HST Cycle 21,
; based upon earlier version for HST Cycle 20, added call
; to new Lyman-alpha-specific code.

; Updated for DR12, Nov 2014

; export BOSS_SPECTRO_REDUX=/Volumes/DR10
; export RUN2D=v5_5_12
; export RUN1D=v5_5_12

; For up-to-the-minute:
; export BOSS_SPECTRO_REDUX=/Volumes/PHY_Bolton_Data3/bolton_data3/SDSS3/BOSS/bossredux_lbl
; export RUN2D=v5_6_0
; export RUN1D=v5_6_0


; Get sky-residual rescaling info:
in_file = FILEPATH('boss_skyrenorm_' + getenv('RUN2D') + '.fits', ROOT_DIR='./data/')
skystruc = mrdfits(in_file,1)

; Get the plates from the skystruc:
plate = skystruc.platelist
mjd = skystruc.mjdlist

nmain = n_elements(plate)

splog, 'Found ', nmain, ' good plates'

; Get the emission lines:
readcol, 'emlines.txt', emwave, emname, format='F,A'

; Get the empirical absorption-line info:
readcol, 'abslines.txt', abswave, absname, format='F,A', comment='#'

; Get empirical skyline info:
readcol, 'skylines.txt', skywave, skyname, format='F,A', comment='#'

; Use Yanmei's templates:
gtfile = getenv('IDLSPEC2D_DIR') + '/templates/spEigenGal-55740.fits'
;gtfile = 'spMLpcaGal-55331-7temp.fits'
galtemp = mrdfits(gtfile,0,galhdr)
ntpix = (size(galtemp))[1]
ntemp = (size(galtemp))[2]
gtloglam = sxpar(galhdr, 'COEFF0') + sxpar(galhdr, 'COEFF1') * findgen(ntpix)


; Loop over plates:
for i = 0L, nmain-1 do begin
;for i = 1001L, 1016 do begin
    filename_h5 = STRCOMPRESS(string(plate[i])+'_'+string(mjd[i])+'.h5',/remove_all)
    ;if filename_h5 ne '5017_55715.h5' then continue
    out_file = FILEPATH(filename_h5, ROOT_DIR='./data/2/')
    ;if file_test(out_file) then begin
    ;    print,'plate',i+1,out_file + 'had been written, Omitted!'
    ;    ;if (plate[i] ne 6759) and (mjd[i] ne 56416) then begin
    ;    continue
    ;    ;endif
    ;endif
    
    ;if (plate[i] ne 6759) or (mjd[i] ne 56416) then begin
    ;    continue
    ;endif
    splog, 'plate number ===', plate[i], 'mjd number =',mjd[i]
    splog, 'Scanning plate ', i+1, ' of ', nmain
; get data:

    readspec, plate[i], mjd=mjd[i], flux=flux, $
      synflux=synflux, $
      invvar=invvar, plug=plug, wave=wave, loglam=loglam, $
      zans=zans, zhdr=zhdr, $
      andmask=andmask, ormask=ormask, /align

; Apply skymask::
;    newivar = skymask(invvar, andmask, ormask)
    newivar = skymask(invvar, andmask)

; Cut down to spectro targeted AND classed galaxies
; with ZWARNING of zero:
    whgal = where( $
     (strtrim(zans.objtype,2) eq 'GALAXY') and $
     (strtrim(zans.class_noqso,2) eq 'GALAXY') and $
     (zans.zwarning_noqso eq 0), ngal)

; Use sky fibers instead, if requested:
    if keyword_set(dosky) then whgal = where(strtrim(zans.objtype,2) eq 'SKY', ngal)

    ;print,'----ngal-----',ngal
    if (ngal gt 0) then begin

       ; Structure to hold info about the plates under consideration:
        fibermask = keyword_set(sdss1) ? replicate(0B, 640) : replicate(0B, 1000)
        fibermask[whgal] = 1B
        plate_struc = {plate: plate[i], mjd: mjd[i], fibermask: fibermask, $
                       wave: 10.^skystruc.skyloglam, med_lnoise: 0. * skystruc.skyloglam}

        flux = flux[*,whgal]
        synflux = synflux[*,whgal]
        invvar = invvar[*,whgal]
        newivar = newivar[*,whgal]
        andmask = andmask[*,whgal]
        ormask = ormask[*,whgal]
        zans = zans[whgal]
        plug = plug[whgal]

        if keyword_set(dosky) then zans.z_noqso *= 0.

; Generate emission/absorption-line mask:
        emask = make_emline_mask(loglam=loglam, z=zans.z_noqso, linelist=alog10([emwave, abswave]), hw=10)

; Kluge observed-frame sky feature masks onto this (zeroing out redshifts):
        emask = emask * make_emline_mask(loglam=loglam, z=0.*zans.z_noqso, linelist=alog10(skywave), hw=10)

; Fit new synflux models and zero out inverse variance where models
; are undefined:
        newsynflux = 0. * synflux
        if (not keyword_set(dosky)) then begin
           for j = 0L, ngal-1 do begin
              newsynflux[*,j] = $
                 galtempfit(objflux=flux[*,j], objivar=(newivar*emask)[*,j], objloglam=loglam, $
                            z=zans[j].z_noqso, tempset=galtemp, temploglam=gtloglam, npoly=3L, tempmask=tempmask)
              invvar[*,j] = invvar[*,j] * tempmask
              newivar[*,j] = newivar[*,j] * tempmask
           endfor
        endif


; Rescale errors:
        res_ivar = rescale_ivar(resflux=(flux-newsynflux), invvar=(invvar*emask), $
              loglam=loglam, skyloglam=skystruc.skyloglam, skyrenorm=skystruc.skyrenorm)

; Now apply the skymasks to this:
        res_ivar = skymask(res_ivar, andmask)
        
        ;write hdf5
        mg_h5_putdata, out_file, 'flux', flux
        mg_h5_putdata, out_file, 'newsynflux', newsynflux
        mg_h5_putdata, out_file, 'res_ivar', res_ivar
        mg_h5_putdata, out_file, 'loglam', loglam
        mg_h5_putdata, out_file, 'zans', zans
        mg_h5_putdata, out_file, 'plug', plug
        
        ;mg_h5_putdata, out_file, 'info', "some other info, see atrribute"
        ;mg_h5_putdata, out_file, 'info.plate', plate[i]
        ;mg_h5_putdata, out_file, 'info.mjd', mjd[i]
        print,out_file+' hdf5 file written finished!'
        

; Scan for hits: multi-line, single-line [OII], Lyman-alpha:
        ;hitsub_lya = findhits_lya(resflux=(flux-newsynflux), resivar=res_ivar, zans=zans, wave=wave, out_noise=out_noise)
        ;if (ngal gt 1) then out_noise = median(out_noise, dimen=2, /even)
        ;plate_struc.med_lnoise = interpol(out_noise, wave, plate_struc.wave)
;
        ;skystr = keyword_set(dosky) ? 'SKY' : ''
;
        ;if keyword_set(hitsub_lya) then begin
        ;   hitmaster_lya = keyword_set(hitmaster_lya) ? [hitmaster_lya, hitsub_lya] : hitsub_lya
        ;endif
        ;
        ;platemaster_lya = keyword_set(platemaster_lya) ? [platemaster_lya, plate_struc] : plate_struc

    endif

endfor

; Write out the master files:

;ofile_lya = 'hitmaster_boss_lya_' + getenv('RUN2D') + skystr + '.fits'
;out_file1 = FILEPATH(ofile_lya, ROOT_DIR='./data/')
;mwrfits, hitmaster_lya, out_file1, /create
;
;ofile_plate = 'platemaster_boss_lya_' + getenv('RUN2D') + skystr + '.fits'
;out_file2 = FILEPATH(ofile_plate, ROOT_DIR='./data/')
;mwrfits, platemaster_lya, out_file2, /create

return
end

