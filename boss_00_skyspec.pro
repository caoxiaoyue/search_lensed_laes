pro boss_00_skyspec

; reconfigured for BOSS from SDSS1-DR7, program by program.
; re-reconfigured for BOSS DR9 & HST Cycle 20
; re-re-reconfigured for BOSS DR10++ & HST Cycle 21
; And yet again for DR12, 2014 November...

; For DR10:
; export BOSS_SPECTRO_REDUX=/Volumes/DR10
; export RUN2D=v5_5_12
; export RUN1D=v5_5_12

; For up-to-the-minute:
; export BOSS_SPECTRO_REDUX=/Volumes/PHY_Bolton_Data3/bolton_data3/SDSS3/BOSS/bossredux_lbl
; export RUN2D=v5_6_0
; export RUN1D=v5_6_0

; For DR12, use v5_7_0
;This work use v5_13_0, eboss data

platelist, plist=plist
wh = where(strtrim(plist.platequality, 2) eq 'good')
plist = plist[wh]

nmain = n_elements(plist)
splog, nmain, ' good plates'

plate = plist.plate
mjd = plist.mjd

readspec, plate[0], mjd=mjd[0], objhdr=objhdr
dloglam = sxpar(objhdr, 'COEFF1')

npwave = 6000L ; a conservative number of pixels
skytot2 = replicate(0., npwave, nmain) ; total squared scaled sky residual vector for plate
skynpix = replicate(0L, npwave, nmain) ; total number of sky pixels contributing
skyzloglam = replicate(0., nmain) ; zero-pixel loglam of this plate's sky vectors in the arrays

; Loop over plates and total up chi^2 as well as
; number of contributiong pixes as a f'n of wavelength
; for each plate.
splog, 'Looping over plates to read sky spectra...'
for i = 0L, nmain-1 do begin
    print, i
    readspec, plate[i], mjd=mjd[i], plug=plug, /silent
    whsky = where(strtrim(plug.objtype, 2) eq 'SKY', nsky)
    if (nsky gt 0) then begin
        readspec, plate[i], plug[whsky].fiberid, mjd=mjd[i], $
          flux=flux, invvar=invvar, loglam=loglam, /align, /silent
        npthis = n_elements(loglam)
        skytot2[0:npthis-1,i] = total(flux^2 * (invvar > 0.), 2)
        skynpix[0:npthis-1,i] = total((invvar gt 0), 2)
        skyzloglam[i] = loglam[0]
    endif
endfor

; Sort out the anchor pixels and shift the individuals within the arrays:

minloglam = min(skyzloglam)
pshift = round((skyzloglam - minloglam) / dloglam)
for i = 0L, nmain-1 do skytot2[*,i] = shift(skytot2[*,i], pshift[i])
for i = 0L, nmain-1 do skynpix[*,i] = shift(skynpix[*,i], pshift[i])

; Identify truly horrendous plates and exclude them from this
; analysis.  The rchi2 cut value is ad-hoc for DR9:
rchi2_plate = total(skytot2, 1) / total(skynpix, 1)
rchi2max = 1.4
useplate = where(rchi2_plate lt rchi2max, n_use)
splog, 'For cut value chi^2_r < ', rchi2max, ':'
splog, 'Using ', n_use, ' of ', nmain, ' plates'
; 

; Compute the overall RMS statistics:
skytot2 = total(skytot2[*,useplate], 2, /double)
skynpix = total(skynpix[*,useplate], 2, /double)
skyloglam = minloglam + dloglam * dindgen(npwave)

skynmax = max(skynpix)
skynmin = 0.05 * skynmax ; threshold at 5% of the maximum total

; (Restrict to the range of pixels with decent coverage):
wh_enough = where(skynpix ge skynmin)
pixlo = min(wh_enough)
pixhi = max(wh_enough)
skytot2 = skytot2[pixlo:pixhi]
skynpix = skynpix[pixlo:pixhi]
skyloglam = skyloglam[pixlo:pixhi]

skyrms = sqrt(skytot2 / skynpix)

ostruc = {skyloglam: skyloglam, skyrms: skyrms, skynpix: long(skynpix), $
 platelist: plate[useplate], mjdlist: mjd[useplate]}
 
ofile = 'boss_skyrms_'+getenv('RUN2D')+'.fits'
out_file = FILEPATH(ofile, ROOT_DIR='./data/')

mwrfits, ostruc, out_file, /create

return
end

