function rdname, ra_deg, dec_deg, raw=raw, texpm=texpm
; rdname: Gives two ra decimal places and one dec.
; Same as rdconv2, except it truncates instead of rounding.
; Set the "raw" keyword if you want to get back one string
; without colons.  Set texpm keyword to get + and -
; between $$ for tex-ing style.

centiseconds = long64(24000.d0 * ra_deg)
hours = centiseconds / 360000
centiseconds = centiseconds - 360000 * hours
minutes = centiseconds / 6000
centiseconds = centiseconds - 6000 * minutes
seconds = centiseconds / 100
centiseconds = centiseconds - 100 * seconds

if (keyword_set(raw)) then begin
  out_string = string(hours, format='(i2.2)') $
   + string(minutes, format='(i2.2)') $
   + string(seconds, format='(i2.2)') + '.' $
   + string(centiseconds, format='(i2.2)')
endif else begin
  ra_string = string(hours, format='(i2.2)') + ':' $
   + string(minutes, format='(i2.2)') + ':' $
   + string(seconds, format='(i2.2)') + '.' $
   + string(centiseconds, format='(i2.2)')
endelse

dec_sign = (dec_deg ge 0) ? '+' : '-'
dec_sign = keyword_set(texpm) ? '$' + dec_sign + '$' : dec_sign
deciarcsec = long64(36000.d0 * abs(dec_deg))
degrees = deciarcsec / 36000
deciarcsec = deciarcsec - 36000 * degrees
arcmin = deciarcsec / 600
deciarcsec = deciarcsec - 600 * arcmin
arcsec = deciarcsec / 10
deciarcsec = deciarcsec - 10 * arcsec

if (keyword_set(raw)) then begin
  out_string = out_string + dec_sign $
   + string(degrees, format='(i2.2)') $
   + string(arcmin, format='(i2.2)') $
   + string(arcsec, format='(i2.2)') + '.' $
   + string(deciarcsec, format='(i1.1)')
endif else begin
  dec_string = dec_sign + string(degrees, format='(i2.2)') + ':' $
   + string(arcmin, format='(i2.2)') + ':' $
   + string(arcsec, format='(i2.2)') + '.' $
   + string(deciarcsec, format='(i1.1)')
endelse

if (keyword_set(raw)) then return, out_string
return, [ra_string, dec_string]
end
