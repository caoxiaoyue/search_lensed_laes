function fit_lya_hitprofile, hm_lya, dosky=dosky

nhits_3 = n_elements(hm_lya)

readspec, hm_lya.plate, hm_lya.fiberid, mjd=hm_lya.mjd, $
 ;tsobj=tsobj_lya, plug=plug_lya, $
 zans=zans_lya, /silent

; Now let's do some fitting of the profiles:
fstruc = replicate({lflux: 0., lflux_err: 0., $
                    lwave: 0., lwave_err: 0., $
                    lsigma: 0., lsigma_err: 0., $
                    lalpha: 0., lalpha_err: 0., $
                    lfit_chi2: 0., lfit_ndof: 0L, $
                    dchi2_skew: 0.}, nhits_3)
hw = 20 ; half-width for line fitting
; Parinfo structures to control skew parameter:
parinfo = replicate({fixed: 0B}, 4)
parinfo[3].fixed = 1B
parinfo2 = replicate({limited: [0B, 0B], limits: [0., 0.]}, 4)
parinfo2[3].limited[*] = 1B
parinfo2[3].limits[0] = -25.
parinfo2[3].limits[1] = 25.
; pixel baseline:
pbase = findgen(2*hw+1)

for i = 0L, nhits_3-1 do begin ;& $
;for i = 0L, 10 do begin ;& $
   print, 'Fitting ', i+1, ' of ', nhits_3 ;& $
   readspec, hm_lya[i].plate, hm_lya[i].fiberid, mjd=hm_lya[i].mjd, $
             flux=flux, synflux=synflux, invvar=invvar, $
             wave=wave, loglam=loglam, znum=zans_lya[i].znum_noqso, /silent ;& $
;   if (n_elements(flux) gt 1) then begin ;& $
   junk = min(abs(wave - hm_lya[i].wave), cpix) ;& $
   wave = wave[cpix-hw:cpix+hw] ;& $
   synflux = synflux[cpix-hw:cpix+hw] ;& $
   if keyword_set(dosky) then synflux = 0. * synflux
   flux = flux[cpix-hw:cpix+hw] - synflux ;& $
   invvar = invvar[cpix-hw:cpix+hw] ;& $
   ; Convert to per-pixel units:
   dlam = wave * 1.0e-4 * alog(10.0) ;& $
   flux = flux * dlam ;& $
   invvar = invvar / dlam^2 ;& $
   spar = [flux[hw] > 1., float(hw), 2., 0.] ;& $
   ft = {flux: flux, invvar: invvar, deviates: 1B} ;& $
   ; First fit with no skewness:
   fpar = mpfit('mpfit_skewnormpix', spar, functargs=ft, parinfo=parinfo, bestnorm=chi2_gauss, /quiet) ;& $
   ; Now free up skewness:
   fpar[3] = 0.1 ;& $
   fpar2 = mpfit('mpfit_skewnormpix', fpar, functargs=ft, parinfo=parinfo2, perror=perror, bestnorm=chi2_skew, /quiet) ;& $
   fstruc[i].lfit_chi2 = chi2_skew ;& $
   fstruc[i].lfit_ndof = long(total(invvar gt 0)) - n_elements(fpar2) ;& $
   fstruc[i].dchi2_skew = chi2_gauss - chi2_skew ;& $
   fstruc[i].lflux = fpar2[0] ;& $
   fstruc[i].lflux_err = perror[0] ;& $
   fstruc[i].lwave = interpol(wave, pbase, fpar[1]) ;& $
   fstruc[i].lwave_err = interpol(dlam, pbase, fpar[1]) * perror[1] ;& $
   fstruc[i].lsigma = fpar2[2] ;& $
   fstruc[i].lsigma_err = perror[2] ;& $
   fstruc[i].lalpha = fpar2[3] ;& $
   fstruc[i].lalpha_err = perror[3] ;& $
;   endif ;& $
endfor

return, fstruc
end
