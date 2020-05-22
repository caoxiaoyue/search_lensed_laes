pro lya_hitplot, hm

readspec, hm.plate, hm.fiberid, mjd=hm.mjd, zans=zans
plotspec, hm.plate, hm.fiberid, mjd=hm.mjd, znum=zans.znum_noqso

print, ' ' + strtrim(string(hm.plate, format='(i4)'),2) + '-' + $
       strtrim(string(hm.mjd, format='(i5)'),2) + '-' + $
       strtrim(string(hm.fiberid, format='(i4)'),2)
print, '  Detection S/N: ', hm.sn
print, '  Model Fit S/N: ', hm.lflux / hm.lflux_err
print, '  Line Flux:     ', hm.lflux, '   +/-', hm.lflux_err
print, '  Wavelength:    ' , hm.lwave, '   +/-', hm.lwave_err
print, '  Sigma (pix):   ', hm.lsigma, '   +/-', hm.lsigma_err
print, '  Alpha_skew:    ', hm.lalpha, '   +/-', hm.lalpha_err
print, '  Fit rchi^2:    ', hm.lfit_chi2 / float(hm.lfit_ndof)
print, '  Skew dchi^2:   ', hm.dchi2_skew

if (hm.fiberid gt 1L) then begin
   readspec, hm.plate, hm.fiberid-1, mjd=hm.mjd, flux=flux, wave=wave
   soplot, wave, flux, color=4
endif
if (hm.fiberid lt 1000L) then begin
   readspec, hm.plate, hm.fiberid+1, mjd=hm.mjd, flux=flux, wave=wave
   soplot, wave, flux, color=4
endif

readspec, hm.plate, hm.fiberid, mjd=hm.mjd, flux=flux, wave=wave, synflux=synflux, znum=zans.znum_noqso
soplot, wave, flux, thick=2
soplot, wave, synflux, color=3

rlam = [1216., 1549., 1909., 2799.]
for i = 1L, 3 do soplot, rlam[i] * hm.wave / rlam[0], 0., ps=2, color=5, syms=2

soplot, hm.wave, 0., ps=6, color=6, thick=2, syms=2
flam = [3728., 4863., 5008., 6565.]
for i = 1L, 3 do soplot, flam[i] * hm.wave / flam[0], 0., ps=2, color=6, syms=2

return
end
