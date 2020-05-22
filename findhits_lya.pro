function findhits_lya, resflux=resflux, resivar=resivar, zans=zans, wave=wave, sdss1=sdss1, $
                       out_noise=out_noise

; Find candidate Lyman-alpha emission lines


; Wrapper routine to targ_linescan_2, implementing parameter
; choices, loops, and write-out packaging.

; This version does a single-line search kluge to the multi-line
; routine, and appends information on veto SNR values for weeding
; out misidentifications of [OII].

;;;;;;;;;
; Up-front choices:
snmin = 6. ; minimum SNR for candidate Lyman-alpha
maxhits = 5 ; maximum number of hits in a spectrum
;;;;;;;;;

; Fake-wavelength stuff for mis-ID with [OII]:
lyalam = 1215.67
o2lam = sqrt(3727.092 * 3729.875)
hblam = 4862.683
o3lam = 5008.239
halam = 6564.614

fakelam = lyalam * [hblam, o3lam, halam] / o2lam

npix = (size(resflux))[1]
nspec = (size(resflux, /n_dimen) eq 2) ? (size(resflux))[2] : 1L

hit_element = {plate: 0L, mjd: 0L, fiberid: 0, wave: 0., $
 ztarg: 0., sn: 0., snfake: fltarr(n_elements(fakelam))}

; Mask to confine attention to the spectral region of interest:
wmask = (wave gt 3600.) * (wave lt 4800.)

; Mask for 5577:
mask5577 = (wave lt 5575.) or (wave gt 5584.)

; Build convolution kernel:
lsigma = 2.2
lksize = 15
lkern = gkern3(lsigma, lksize, 0.)

; Initialize variables to accumulate the output:
blockid = -1L 
specid = -1L
waveid = 0.
snhit = 0.
fakesn = 0. * fakelam

; To convert to per-pixel baseline:
dlam = wave * 1.0e-4 * alog(10.0)
; To store the output line-flux noise levels:
out_noise = 0. * resflux

for i = 0L, nspec-1 do begin
;   print, i
   ; Do the convolution filtering:
   thisflux = resflux[*,i] * dlam
   thisinvvar = mask5577 * resivar[*,i] / dlam^2
   filter_spec, kernel=lkern, flux=thisflux, invvar=thisinvvar, sn=sn, noise=thisnoise
   out_noise[*,i] = thisnoise
   ; Find things that qualify as "hits":
   htest = (sn * wmask) ge snmin
   ; Identify contiguous hit blocks:
   bi = block_id(htest, lengths=lengths, nblocks=nhits)
   if ((nhits gt 0) and (nhits le maxhits)) then begin
      sn_this = fltarr(nhits)
      ; Reduce hit-blocks to most significant pixel:
      for j = 0, nhits - 1 do begin
         snsub = sn[bi[j]:bi[j]+lengths[j]-1]
         maxval = max(snsub, maxid)
         bi[j] = bi[j] + maxid
         sn_this[j] = maxval
      endfor
; Compute fake-wavelength SNRs:
      fakethis = fltarr(n_elements(fakelam), n_elements(bi))
      for jj = 0, n_elements(bi)-1 do begin
          obsfakelam = (wave[bi[jj]] / lyalam) * fakelam
          fakethis[*,jj] = interpol(sn, wave, obsfakelam) $
            * (obsfakelam gt min(wave)) * (obsfakelam lt max(wave))
      endfor
      fakesn = [[fakesn], [fakethis]]
      blockid = [blockid, bi]
      specid = [specid, i+0*bi]
      waveid = [waveid, wave[bi]]
      snhit = [snhit, sn_this]
   endif
endfor

;i=0
;i++
;splot, wave, flux[*,specid[i]]
;soplot, wave, newsynflux[*,specid[i]], color=3
;soplot, waveid[i] * [1.,1.], [0.,5.], color=6

nhits = n_elements(blockid) - 1
if (nhits gt 0) then begin
  hitsub = replicate(hit_element, nhits)
  hitsub.fiberid = zans[specid[1:*]].fiberid
  hitsub.plate = zans[specid[1:*]].plate
  hitsub.mjd = zans[specid[1:*]].mjd
;  hitsub.z = zzs[1:*]
  hitsub.wave = waveid[1:*]
  hitsub.ztarg = keyword_set(sdss1) ? zans[specid[1:*]].z : zans[specid[1:*]].z_noqso
  hitsub.sn = snhit[1:*]
  hitsub.snfake = fakesn[*,1:*]
endif else hitsub = 0

return, hitsub

end

