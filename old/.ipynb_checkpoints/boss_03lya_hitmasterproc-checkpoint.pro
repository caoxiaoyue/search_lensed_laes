pro boss_03lya_hitmasterproc, plotwait=plotwait, noplot=noplot

; Program to make cuts on the "hitmaster" files and supplement them
; with other photometric and spectroscopic info.
;
; Hack version to deal with Lyman-alpha candidates specifically.
; bolton@utah 2013jan

; Read in the two files.
; Eventually some version-stamp-chceking should be implemented.
;
; Updated for DR12 and general version handling, Jan 2015
; Note: x-forwarding from chpc computers now seems to need idl 8.1...
;
splog, 'Trimming emission-line hit lists...'
hitfile_lya = 'hitmaster_boss_lya_' + getenv('RUN2D') + '.fits'
hitfile_lya_sky = 'hitmaster_boss_lya_' + getenv('RUN2D') + 'SKY.fits'
hm_lya = mrdfits(hitfile_lya,1)
hm_lya_sky = mrdfits(hitfile_lya_sky,1)

; Wavelength-histogramming cuts on the one-line list:
obswave = hm_lya.wave
obswave_sky = hm_lya_sky.wave
restwave = obswave / (1. + hm_lya.ztarg)
restwave_sky = obswave_sky / (1. + hm_lya_sky.ztarg)
if (not keyword_set(noplot)) then splot, alog10(obswave), alog10(restwave), ps=3, color=2, $
 title='One-line hits', xtitle='log10(obswave/Ang)', $
 ytitle='log10(restwave/Ang)'
if keyword_set(plotwait) then wait, plotwait

; Let's do the restwave histo first:
dloglam = 3. * 0.0001 ; three times the SDSS redux pixel unit
minloglam = min(alog10(restwave))
maxloglam = max(alog10(restwave))
nbins = ceil((maxloglam - minloglam) / dloglam)
restbase = dloglam * (findgen(nbins) + 0.5) + minloglam
resthist = histogram(alog10(restwave), min=minloglam, $
 nbins=nbins, binsiz=dloglam)

; Wide median filter:
histwid = 50
histfilt = median(resthist, 2*histwid+1, /even) > 1.
histfilt[0:histwid] = histfilt[histwid+1]
histfilt[nbins-histwid-1:nbins-1] = histfilt[nbins-histwid-2]
if (not keyword_set(noplot)) then splot, restbase, resthist, ps=10, color=4, $
 title='Target galaxy rest-frame hitwave histogram', $
 xtitle='bin number', ytitle='number of hits'
if (not keyword_set(noplot)) then soplot, restbase, histfilt, ps=10

; To identify wavelengths of bad locations:
;splot, 10.^restbase, resthist, ps=10, color=4, $
; title='Target galaxy rest-frame hitwave histogram', $
; xtitle='bin number', ytitle='number of hits'

; I'll use root-n even though this is Poisson not Gauss.
sigthresh = 4.
histmask = resthist gt (histfilt + sigthresh*sqrt(histfilt))
if (not keyword_set(noplot)) then soplot, restbase, resthist * (1B - histmask), ps=10, color=1
if keyword_set(plotwait) then wait, plotwait

; Grow the mask a bit:
ngrow = 1
histmask = smooth(float(histmask), 2*ngrow+1) gt 0.

; Figure out who gets cut based on this mask:
cutit = interpol(float(histmask), restbase, alog10(restwave)) gt 0.5
cutit_sky = interpol(float(histmask), restbase, alog10(restwave_sky)) gt 0.5

;hm_one = hm_one[where(cutit eq 0)]
hm_lya = hm_lya[where(cutit eq 0)]
nhits_1 = n_elements(hm_lya)
hm_lya_sky = hm_lya_sky[where(cutit_sky eq 0)]
nhits_1_sky = n_elements(hm_lya_sky)

; Make sure this worked:
obswave = hm_lya.wave
obswave_sky = hm_lya_sky.wave
restwave = obswave / (1. + hm_lya.ztarg)
restwave_sky = obswave_sky / (1. + hm_lya_sky.ztarg)
if (not keyword_set(noplot)) then splot, alog10(obswave), alog10(restwave), ps=3, color=2, $
 title='One-line hits', xtitle='log10(obswave/Ang)', $
 ytitle='log10(restwave/Ang)'
if keyword_set(plotwait) then wait, plotwait

; OK, now cut based on observed-frame wavelength:
dloglam = 3. * 0.0001 ; three times the SDSS redux pixel unit
minloglam = min(alog10(obswave))
maxloglam = max(alog10(obswave))
nbins = ceil((maxloglam - minloglam) / dloglam)
obsbase = dloglam * (findgen(nbins) + 0.5) + minloglam
obshist = histogram(alog10(obswave), min=minloglam, $
 nbins=nbins, binsiz=dloglam)

; Wide median filter:
histwid = 50
histfilt = median(obshist, 2*histwid+1, /even) > 1.
histfilt[0:histwid] = histfilt[histwid+1]
histfilt[nbins-histwid-1:nbins-1] = histfilt[nbins-histwid-2]

if (not keyword_set(noplot)) then splot, obsbase, obshist, ps=10, color=4, $
 title='Observer-frame hitwave histogram', $
 xtitle='bin number', ytitle='number of hits'
if (not keyword_set(noplot)) then soplot, obsbase, histfilt, ps=10

; To identify offensive wavelengths::
;splot, 10.^obsbase, obshist, ps=10, color=4, $
; title='Observer-frame hitwave histogram', $
; xtitle='bin number', ytitle='number of hits'

; I'll use root-n even though this is Poisson not Gauss.
sigthresh = 4.
histmask = obshist gt (histfilt + sigthresh*sqrt(histfilt))
if (not keyword_set(noplot)) then soplot, obsbase, obshist * (1B - histmask), ps=10, color=1
if keyword_set(plotwait) then wait, plotwait

; Grow the mask a bit:
ngrow = 1
histmask = smooth(float(histmask), 2*ngrow+1) gt 0.

; Figure out who gets cut based on this mask:
cutit = interpol(float(histmask), obsbase, alog10(obswave)) gt 0.5
cutit_sky = interpol(float(histmask), obsbase, alog10(obswave_sky)) gt 0.5

;hm_one = hm_one[where(cutit eq 0)]
hm_lya = hm_lya[where(cutit eq 0)]
hm_lya_sky = hm_lya_sky[where(cutit_sky eq 0)]

; Make sure this worked:
obswave = hm_lya.wave
obswave_sky = hm_lya_sky.wave
restwave = obswave / (1. + hm_lya.ztarg)
restwave_sky = obswave_sky / (1. + hm_lya_sky.ztarg)
if (not keyword_set(noplot)) then splot, alog10(obswave), alog10(restwave), ps=3, color=2, $
 title='One-line hits', xtitle='log10(obswave/Ang)', $
 ytitle='log10(restwave/Ang)'
if keyword_set(plotwait) then wait, plotwait

; Now cut (aggressively!) based upon various mis-identifications:
nhits_2 = n_elements(hm_lya)
nhits_2_sky = n_elements(hm_lya_sky)
sntest = 2.5
is_o2 = total(hm_lya.snfake ge sntest, 1) gt 0
is_o2_sky = total(hm_lya_sky.snfake ge sntest, 1) gt 0

; Cut it down:
hm_lya = hm_lya[where(is_o2 eq 0)]
hm_lya_sky = hm_lya_sky[where(is_o2_sky eq 0)]

nhits_3 = n_elements(hm_lya)
nhits_3_sky = n_elements(hm_lya_sky)

; More conservative SNR cut:
;sn_newmin = -1.e6 ; (currently not enabled)
;hm_one = hm_one[where(hm_one.sno2 gt sn_newmin)]
;hm_lya = hm_lya[where(hm_lya.sn gt sn_newmin)]

fstruc = fit_lya_hitprofile(hm_lya)
fstruc_sky = fit_lya_hitprofile(hm_lya_sky, /dosky)

; Append the fit info:
hm_lya = struct_addtags(hm_lya, fstruc)
hm_lya_sky = struct_addtags(hm_lya_sky, fstruc_sky)

; Output files:
; (probably want to automate output-file generation based upon
; hitmaster filenames eventually.)
splog, 'Writing output files...'

ofile_lya = 'boss_lyaline_' + getenv('RUN2D') + '.fits'
mwrfits, hm_lya, ofile_lya, /create
ofile_lya_sky = 'boss_lyaline_' + getenv('RUN2D') + 'SKY.fits'
mwrfits, hm_lya_sky, ofile_lya_sky, /create

return
end
