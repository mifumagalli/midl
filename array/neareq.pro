;+
;PURPOSE
;	sometimes when dealing with roundoff errors
;	numbers arent equal when they really are meant to be
;	This applies a tolerance cut
;SYNTAX
;	bool=neareq(var1, var2[ tol=tol])
;INPUTS
;	var1, var2: the value of a number to be compared
;	tol: the tolerance [default to 1d-5
;
;Written by R. da Silva, UCSC, 11-30-09
;-

FUNCTION neareq, var1, var2, tol=tol
if not keyword_set(tol) then tol=1d-5

return, abs(var1-var2) LT tol
end
