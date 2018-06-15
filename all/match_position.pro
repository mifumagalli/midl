;per un set di coordinate (X,Y) and (Xbis,Ybis), 
;ind, indbis restituisce gli indici che hanno match
;tol set the tolerance in the comparison


PRO match_position, X , Y, Xbis, Ybis, ind, indbis,  TOL=tol 


Num=N_ELEMENTS(X)
Numbis=N_ELEMENTS(Xbis)

ind=MAKE_ARRAY(Numbis,/INTEGER,VALUE=-99)
indbis=MAKE_ARRAY(Numbis,/INTEGER,VALUE=-99)

ibis=0
match=0

WHILE(ibis LT Numbis) DO BEGIN

i=0
WHILE(i LT Num) DO BEGIN
IF(ABS(Xbis[ibis]-X[i]) LT tol AND ABS(Ybis[ibis]-Y[i]) LT tol) THEN BEGIN
ind[match]=i 
indbis[match]=ibis 
match=match+1
ENDIF

i=i+1
ENDWHILE

ibis=ibis+1
ENDWHILE

ind=ind[0:match-1]
indbis=indbis[0:match-1]

;print, "Number match: ", match

END
