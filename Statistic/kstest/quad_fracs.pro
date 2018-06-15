;+
; Compute the fraction of points in each quadrant given the center
; point xc,yc.  There's an ambiguity which quadrant to stick
; things if xc,yc lands on a data point.  If that happens, quad_fracs
; returns an array where it computes the fractions putting the xc,yc point in each quadrant. 
; the equations are not accurate if N < 20 or the returned probability
; is > 0.2.  if the returned is greater than 0.2, then it is a
; sign that you can treat them as drawn from the same distribution.
;
; Stolen from http://www.astro.washington.edu/users/yoachim/code.php
;
;
;-


function quad_fracs, xc,yc,xx,yy

;;find in which quadrant you are
r1=where(xx gt xc and yy gt yc)
r2=where(xx lt xc and yy gt yc)
r3=where(xx lt xc and yy lt yc)
r4=where(xx gt xc and yy lt yc)

nt=double(n_elements(xx))

;;assign the number
if max(r1) eq -1 then n1=0 else n1=n_elements(r1)
if max(r2) eq -1 then n2=0 else n2=n_elements(r2)
if max(r3) eq -1 then n3=0 else n3=n_elements(r3)
if max(r4) eq -1 then n4=0 else n4=n_elements(r4)


;;make the fractions
fracs=[n1,n2,n3,n4]/nt

;;check if xc,yc hit a point
if total(fracs) ne 1.D then begin 
    fracs=rebin([n1,n2,n3,n4],4,4)
    diag=[[1.,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
    fracs=(fracs+diag)/nt
endif

return, fracs


end

