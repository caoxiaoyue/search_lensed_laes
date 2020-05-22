pro aptingest

infile = 'boss_lyaline_v5_7_0_hst_candidates.fits'
st = mrdfits(infile,1)
;p_st = mrdfits(infile, 2)
ncan = n_elements(st)

srt = sort(st.ra)
st = st[srt]
;p_st = p_st[srt]

i_magni = st.cmodelmag[3]

; Dump these to an ingest table (which we still need to hand-edit):
oname = strarr(ncan)
for i = 0L, ncan-1 do oname[i] = 'SDSSJ'+rdname(st[i].ra, st[i].dec, /raw)
onum = string(lindgen(ncan)+1, format='(i2)')
o_ra = string(st.ra, format='(f10.6)')
o_dec = string(st.dec, format='(f10.6)')
o_iband = '"i = ' + string(i_magni, format='(f4.1)') + '"'
o_string = onum + ', ' + oname + ', ' + o_ra + ', ' + o_dec + ', , , ' + o_iband
;o_string = onum + ', ' + oname + ', ' + oname + ', , , ' + o_iband
;ingest_file = 'apt_ingest.txt'
ingest_file = 'apt_ingest_c23.targets'
spawn, 'rm -f ' + ingest_file
openw, 11, ingest_file
for i = 0L, ncan-1 do printf, 11, o_string[i]
close, 11

return
end
