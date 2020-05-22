;+
;
; NAME: galtempfit
;
; PURPOSE: perform fit of rest-frame template set to observed galaxy
;   spectrum
;
; USAGE:
;   yfit = galtempfit(objflux=objflux, objivar=objivar, objloglam=objloglam, z=z, $
;                     tempset=tempset, temploglam=temploglam, npoly=npoly)
;
; ARGUMENTS (all keywords, all mandatory):
;   objflux: object flux vector
;   objivar: object inverse-variance vector
;   objloglam: object log10-wavelength vector
;   z: galaxy redshift
;   tempset: ntpix X ntemp array of eigen-galaxy spectra
;   temploglam: template log10-wavelength (rest-frame) baseline vector
;   npoly: number of polynomial terms to incorporate in addition
;
; RETURNS:
;   Fitted model to the spectrum (that's all for now).
;
; Copyright 2008 Adam S. Bolton
;
; Added zeroing-out of non-covered red-end rest wavelengths
; via tempmask (also output), ASBfeb2010
;
;-

function galtempfit, objflux=objflux, objivar=objivar, objloglam=objloglam, z=z, $
 tempset=tempset, temploglam=temploglam, npoly=npoly, tempmask=tempmask

ntemp = (size(tempset))[2]
npix = n_elements(objloglam)

; Interpolate the templates to the frame of the galaxy:

thistemp = fltarr(npix, ntemp)
restloglam = objloglam - alog10(1. + z)
maxrest = max(temploglam)
for i = 0L, ntemp-1 do thistemp[*,i] = $
 interpol(tempset[*,i], temploglam, restloglam)

; Normalize the templates:
for i = 0L, ntemp-1 do thistemp[*,i] = thistemp[*,i] / sqrt(mean(thistemp[*,i]^2))

; Add on some polynomial terms:
pbase = findgen(npix) / float(npix-1)
polyset = fpoly(pbase, npoly)
thistemp = [[thistemp], [polyset]]

; Do the fitting:
tempmask = (restloglam le maxrest)
ithistemp = thistemp * ((objivar * tempmask) # replicate(1., ntemp+npoly))
alpha = transpose(ithistemp) # thistemp
beta = transpose(ithistemp) # objflux
svdc, alpha, w, u, v
coeff = svsol(u, w, v, beta)
;ialpha = invert(alpha)
;coeff = ialpha # beta
yfit = thistemp # coeff

return, yfit
end
