;compute the offset in arcsec between two objects  1-->2
;(positive axis N,E)

;pos1 --> position in hh:mm:ss and dd:mm:ss or x,y  
;pos2 --> position in hh:mm:ss and dd:mm:ss or x,y 


function arcsec_offset, pos1, pos2

;generate deg values

x_radec, pos1[0], pos1[1], rag1, deg1
x_radec, pos2[0], pos2[1], rag2, deg2


;bring at the same delta of pos2

deg_off=3600.*(deg2-deg1)

;get the delta correction
delta=cos(deg2*!PI/180.)
rag_off=3600.*(rag2-rag1)*delta          

return, [rag_off,deg_off]
end
