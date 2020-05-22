; The inspection note file
; Change it to your own file name
note_file='lya_inspect_notes_Shu_v1.txt'

; Restrict to detection significance greater than 8
hm = mrdfits('boss_lyaline_' + getenv('RUN2D') + '.fits',1)
wh = where(hm.sn ge 8.0)
hm = hm[wh]

readcol, note_file, plate, mjd, fiber, wave, grade, comment, format='(I4, L5, I4, F6.1, A, A)', stringskip='#'

; Select all the yes'es from visual inspection
wh_yes=where(strmatch(grade, '*Yes*')) ; 189 yes'es

plate=plate[wh_yes]
mjd=mjd[wh_yes]
fiber=fiber[wh_yes]
wave=wave[wh_yes]
grade=grade[wh_yes]
comment=comment[wh_yes]

hm=hm[wh_yes]

; Get ra and dec information
readspec, plate, fiber, mjd=mjd, zans=zans
ra=zans.plug_ra
dec=zans.plug_dec

; Remove duplicates with identical ra, dec, and wave
ind=[]

for i=0L, n_elements(plate)-1 do begin &$
	print, i &$
	j=i+1 &$
	flag=0 &$
	while (flag eq 0 and j lt n_elements(plate)-1) do begin &$
		if (abs(ra[i]-ra[j]) lt 10.^(-6) and abs(dec[i]-dec[j]) lt 10.^(-6) and abs(wave[i]-wave[j]) lt 10.^(-6)) then flag=1 &$
		j=j+1 &$
	endwhile &$
	if (flag eq 1) then continue else ind=[ind, i] &$
endfor

help, ind ; 186 unique emission lines

plate_unique=plate[ind]
fiber_unique=fiber[ind]
mjd_unique=mjd[ind]
hm_unique=hm[ind]

i=0
; Check 1-D spectra from all subexposures
; real signal should appear in all of them
splog, 'The detected emission line is at ', hm[i].wave
plotspec, plate_unique[i], fiber_unique[i], mjd=mjd_unique[i], xrange=[hm[i].wave-10, hm[i].wave+10], /allexp
; Check 2-D spectra from all subexposures
; to see if there is any CCD defects
atvspec, plate_unique[i], fiber_unique[i], wave=hm[i].wave, mjd=mjd_unique[i]

i=i+1

; All 186 candidates pass the above two tests.

