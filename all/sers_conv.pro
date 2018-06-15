;Sersic function convolved with a moffat function. 
;No analitic derivative
;4 parametri per sersic, tre per moff e uno per normalizzazione convoluzione

PRO sers_conv, X, P, Funct, Deriv

COMMON Moffat_paramteri

;calcola le funzioni
Sers=P[0]*EXP(-P[1]*((X/P[2])^(1./P[3])-1))
Moff=PAR_psf[0]/(1+(X/PAR_psf[1])^2)^PAR_psf[2]

;fai la convoluzione
;Funct=CONVOL(Sers,Moff,/NORMALIZE,CENTER=0,/EDGE_TRUNCATE )
Funct=CONVOL(Sers,Moff,CENTER=0,/EDGE_TRUNCATE )
Deriv=[0]
END
 

 
