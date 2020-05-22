;
; Script for inspecting Ly-a candidates...
;

hm = mrdfits('boss_lyaline_' + getenv('RUN2D') + '.fits',1)
wh = where(hm.sn ge 8.0)
hm = hm[wh]
n_new = n_elements(wh)


; Generate a new file for recording inspections:
grade = replicate('_', n_new)
comm = replicate('_', n_new)
nsys = n_new
ofile = 'lya_inspect_notes.txt'
spawn, 'rm -f ' + ofile
openw, 11, ofile
for i = 0L, nsys-1 do printf, 11, hm[i].plate, hm[i].mjd, hm[i].fiberid, hm[i].wave, $
  grade[i], comm[i], format = '(1x,i4,2x,i5,2x,i4,2x,f6.1,2x,a,2x,a)'
close, 11

; Make sure you rename this file before recording inspections in it!!
; (To avoid inadventent clobbering.)

; Manual-exection loop for inspecting:
i = -1L
i ++ & lya_hitplot, hm[i]

