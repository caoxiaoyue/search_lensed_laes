;+
;
; NAME:
;  skew_normal
;
; PURPOSE:
;  Evaluate skew normal distribution
;
; USAGE:
;  y = skew_normal(x, alpha=alpha)
;
;-

function skew_normal, x, alpha=alpha
  alpha_loc = keyword_set(alpha) ? double(alpha) : 0.d0
  gaussfact = exp(-0.5d0 * x^2) / sqrt(2.d0 * !pi)
  erffact = 0.5d0 * (1.d0 + erf(alpha_loc * x / sqrt(2.d0)))
  y = 2.d0 * gaussfact * erffact
  dtype = size(x, /type)
  if (dtype eq 4) then y = float(y)
  return, y
end

