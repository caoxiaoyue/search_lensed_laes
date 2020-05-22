pro filter_spec, kernel=kernel, flux=flux, invvar=invvar, $
 lflux=lflux, sn=sn, noise=noise

; IDL program to fit a normalized kernel-filter to 
; a spectrum (the residual spectrum, for the "rogue
; line" project) with flux array "flux".
; flux is assumed to have units of energy/cm^2/s/PIXEL,
; not /Angstrom, so make sure you convert the standard
; SDSS flux and invvar before passing them, perhaps as
; follows:
;   dlam = wave * 1.0e-4 * alog(10.0)
;   flux = flux * dlam
;   invvar = invvar / dlam^2
; kernel size should be odd, and kernel should be normalized!
; lflux, sn, and noise are the best-fit line fluxes,
; s/n on these line fluxes, and noises
; of the fit at each position.


; Perform the necessary convolutions:
con1 = convol((flux * invvar), kernel)
con2 = convol(invvar, (kernel^2))

; The line flux and line-flux s/n arrays:
lflux = con1 / (con2 > 1.e-10)
sn = con1 / (sqrt(con2 > 1.e-10) > 1.e-10)
noise = 1.0 / (sqrt(con2 > 1.e-10) > 1.e-10)

end
