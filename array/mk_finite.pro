;+
;PURPOSE
;	to set any non finite values to 0
;SYNTAX
;	res=mk_finite(arr)
;INPUTS
;	arr: an array
;OUPUT
;	res: the same array but non finite values
;		(found via finite function) are set to 0
;Written by R. da Silva, UCSC, 4-16-10
;-
FUNCTION mk_finite, arr
whbad=where(finite(arr) EQ 0, ct) ;find non finite values
if ct NE 0 then arr[whbad]=0 ;replace with zero
return, arr
end
