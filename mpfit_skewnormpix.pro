function mpfit_skewnormpix, par, npix=npix, flux=flux, invvar=invvar, $
   subsamp=subsamp, deviates=deviates

; par[0]: normalization
; par[1]: centroid
; par[2]: sigma
; par[3]: alpha

if (not keyword_set(npix)) then npix = n_elements(flux)
if (not keyword_set(subsamp)) then subsamp = 10L
subbase = (dindgen(subsamp * npix) + 0.5d0) / double(subsamp) - 0.5d0
subbase = reform(subbase, subsamp, npix)
subbase = (subbase - par[1]) / par[2]
subfunc = par[0] * skew_normal(subbase, alpha=par[3]) / par[2]
y_out = total(subfunc, 1) / double(subsamp)

if keyword_set(deviates) then y_out = (flux - y_out) * sqrt(invvar)

return, y_out
end
