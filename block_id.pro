function block_id, vec, lengths=lengths, nblocks=nblocks

; Function to identify contiguous blocks of non-zero entries within a
; 1D array "vec" and return the indices of their beginnings.  Lengths
; of the blocks are returned in "lengths" and the number of blocks in
; "nblocks".

veclen = n_elements(vec)

blockid = -1L
lengths=-1L

inblock = 0
nblocks=0L
for i = 0L, veclen-1 do begin
  if (vec[i] eq 0.) then begin
    inblock = 0
  endif else begin
    if (inblock eq 0) then begin
      inblock = 1
      nblocks = nblocks + 1
      blockid = [blockid, i]
      lengths = [lengths, 1]
    endif else begin
      lengths[nblocks] = lengths[nblocks] + 1
    endelse
  endelse
endfor

if (nblocks gt 0) then begin
  blockid = blockid[1:nblocks]
  lengths = lengths[1:nblocks]
endif

return, blockid

end
