;+
;PURPOSE
;	make an array
;SYNTAX
;	res=mkarr(low, hi, step)
;INPUTS
;	low: lowest value of array
;	hi: highest value of array
;	step: step of array
;KEYWORDS
;	/nbins: if set then step corresponds
;		to the number of bins
;-

function mkarr, low1, hi, step, nbins=nbins
if N_PARAMS() EQ 1 then begin
low=low1[0]
hi=low1[1]
step=low1[2]
endif else low=low1
if keyword_set(nbins) EQ 0 then begin
nsteps=(hi-low)/step+1
arr=dindgen(round(nsteps))*step+low
endif else begin

nsteps=step
step0=double(hi-low)/(nsteps-1)
arr=dindgen(round(nsteps))*step0+low

endelse

return, arr
end
