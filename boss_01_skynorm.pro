pro boss_01_skynorm, noplot=noplot

; This routine "continuum"-normalizes the sky residual RMS spectrum,
; to derive the spectrum by which error estimates should be rescaled.
; This procedure allows us to capture the full significance difference
; between on-skyline and off-skyline without any undesired reduction
; of the error estimates in the "continuum" region.
;
; Copyright 2008 Adam S. Bolton
; Trivial mods for BOSS DR9: 2011 Adam S. Bolton
; Trivial mods for BOSS DR10: 2013 Adam S. Bolton
; Trivial mods for BOSS DR12: Nov. 2014 Adam S. Bolton
in_file = FILEPATH('boss_skyrms_' + getenv('RUN2D') + '.fits', ROOT_DIR='./data/')
st = mrdfits(in_file,1)

dloglam = (st.skyloglam[1] - st.skyloglam[0])
logzreff = 3.5791 ; reference zero pixel for when we actually did this
offset = round((st.skyloglam[0] - logzreff) / dloglam)
; positive "offset" means the current zero pixel is redder than
; the zero pixel of the original case.  That means we should
; subtract "offset" from the dialed-in indices, because as they stand
; they will point to wavelengths that are too red.


; Choices have been slightly re-kludged for BOSS.
; Do
;   splot, findgen(npix) + offset, st.skyrms
; to see the sky noise scale spectrum against the index baseline
; referenced in the mask definition given here.
npix = n_elements(st.skyloglam)
masklo = ([-290L, -183, 593, 1311, 1575, 1663, 1903, 2162, 2326, 2543, 2796, 3070, 3377, 3627, 3880, 4067, 4328] - offset) > 0
maskhi = ([-240L, -168, 611, 1360, 1592, 1688, 1923, 2259, 2414, 2701, 3025, 3305, 3536, 3820, 4046, 4313, npix-1] - offset) > 0

nmask = n_elements(masklo)
masklo[0] = 0L
maskhi[nmask-1] = npix-1

smask = replicate(1B, npix)
for i = 0L, nmask-1 do smask[masklo[i]:maskhi[i]] = 0B
;try to do not fit the region with wavelength < 3700.
wh_3700 = where(st.skyloglam lt alog10(3700))
smask[wh_3700] = 0B

bkpt = [0L, 1000, 1833, 1875, 1925, 1981, 2040, 2100, 3000, 3500, npix-1] - offset
nbkpt = n_elements(bkpt)
bkpt[0] = 0L
bkpt[nbkpt-1] = npix-1

bkpt = st.skyloglam[bkpt]

sset = bspline_iterfit(st.skyloglam, st.skyrms, invvar=float(smask), nord=3, $
 bkpt=bkpt, yfit=rmsfit, maxiter=0)

;if (not keyword_set(noplot)) then splot, findgen(n_elements(st.skyrms))+offset, st.skyrms, yrange=[0,3]
;if (not keyword_set(noplot)) then soplot, findgen(n_elements(st.skyrms))+offset, st.skyrms*smask, color=1
;if (not keyword_set(noplot)) then soplot, findgen(n_elements(st.skyrms))+offset, rmsfit, color=2

renorm = (st.skyrms / rmsfit) > 1

if (not keyword_set(noplot)) then splot, 10.^st.skyloglam, st.skyrms, ps=10, yrange=[0,2]
if (not keyword_set(noplot)) then soplot, 10.^st.skyloglam, st.skyrms*smask, ps=10, color=2
if (not keyword_set(noplot)) then for i = 0L, nbkpt-1 do soplot, (10.^bkpt[i])*[1.,1.], [.8,1.2], color=6
if (not keyword_set(noplot)) then soplot, 10.^st.skyloglam, rmsfit, color=1
if (not keyword_set(noplot)) then soplot, 10.^st.skyloglam, renorm, ps=10, color=4

st = struct_addtags(st, {rmsfit: rmsfit, skyrenorm: renorm, smask: smask})


ofile = 'boss_skyrenorm_' + getenv('RUN2D') + '.fits'
out_file = FILEPATH(ofile, ROOT_DIR='./data/')
mwrfits, st, out_file, /create

return
end
