;+
;PURPOSE
;	to produce a compressed string
;	returns  strcompress(string(var, _extra=_extra), /remove_all)
;SYNTAX
;	str=rstring(arr, _extra=_extra)
;INPUTS
;	arr: array you wish to convert
;	_extra: extra keywords for string()
;
;Written by R. da Silva, UCSC, Fall 2009
;-
function rstring, var, _extra=_extra

return, strcompress(string(var, _extra=_extra), /remove_all)


end
