;+
;PURPOSE
;	to do what uniq should do... return a list of sorted uniq values
;	NOT the indices
;SYNTAX
;	u=unique(a)
;KEYWORDS:
;	/fin: set to only return finite values
;
;Written by R.da Silva, 2-12-09, UCSC
;-


FUNCTION unique, a1, fin=fin, index=index
a=a1

sa=sort(a)
a=a[sa]
una=uniq(a)
a=a[una]
index=sa[una]

if keyword_set(fin) then a=a[where(finite(A))]

return, a
end
