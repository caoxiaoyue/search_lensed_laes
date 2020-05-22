; gkern3(sigma, length, offset, nonorm=nonorm)
;
; A simple, finite-range, normalized Gaussian kernel generator.
; This one integrates over the bin, unlike g_kern.
; Also allows for offset of profile center from kernel center.
;
; Arguments:
;    sigma:   standard deviation of the kernel, in units of the
;             array index
;    length:  total length of the kernel array, which should be
;             odd (or it will be made odd)
;
; Output is an array with an odd number of elements, symmetric
; about the middle element, and normalized.
; Set keyword "nonorm" to prevent the kernel from being
; renormalized over its extent.

function gkern3, sigma, length, offset, nonorm=nonorm


   use_length = abs(fix(length))
   if ((use_length MOD 2) EQ 0) then use_length = use_length + 1
   kern_array = findgen(use_length) - offset
   kern_array = kern_array - float((use_length - 1) / 2)
   print,kern_array
   kern_array = gaussint((kern_array + 0.5) / sigma) - $
                 gaussint((kern_array - 0.5) / sigma)
   if (not keyword_set(nonorm)) then $
    kern_array = kern_array / total(kern_array)
   return, kern_array

end
