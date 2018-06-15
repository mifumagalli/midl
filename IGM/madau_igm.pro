;correct the magnitude in one filter according to IGM blanket (Madau 1995)
;This is the approximated treatement following the integral in eq (17 e 18)
;lambda---> wavelenght array of a filter trasmission curve
;transm---> filter transmission curve 
;zemi---> redhsift at which you want the IGM correction
;cIGM --> on output, Flux=Observedflux*cIGM


PRO madau_igm, lambda, transm, zemi, cIGM


;define frequency
ll_obs=912*(1+zemi)
lalpha_obs=1216*(1+zemi)
lbeta_obs=1026*(1+zemi)


;calcolo coefficienti DA e DB

A2=3.6D-3
A3=1.7D-3
A4=1.2D-3
A5=9.3D-4

lalpha=1216
lbeta=1026
lgamma=973
ldelta=950


;DA   
lobs=mkarr(1050.*(1+zemi),1170.*(1+zemi),0.01)
deltalobs=lobs-SHIFT(lobs,1)
deltalobs[0]=lobs[1]-lobs[0]
Argum=EXP(-A2*(lobs/lalpha)^(3.46))
;calcolo 1-DA (eq 17, Madau 1995)
CoeffDA=(1./(120.*(1+zemi)))*TOTAL(Argum*deltalobs)




;DB  
lobs=mkarr(920.*(1+zemi),1015.*(1+zemi),0.01)
deltalobs=lobs-SHIFT(lobs,1)
deltalobs[0]=lobs[1]-lobs[0]
Sum1=A3*(lobs/lbeta)^(3.46)
Sum2=A4*(lobs/lgamma)^(3.46)
Sum3=A5*(lobs/ldelta)^(3.46)
Argum=EXP(-Sum1-Sum2-Sum3)
;calcolo 1-DB (eq 18, Madau 1995)
CoeffDB=(1./(95.*(1+zemi)))*TOTAL(Argum*deltalobs)


;make mask with IGM coefficient IGM
Mask=FLTARR(N_ELEMENTS(lambda))

;prima di 912
LLzero=where(lambda LE ll_obs, nll)
IF(nll GT 0) THEN Mask[LLzero]=0.
;da 912 a Lbeta
LLbeta=where(lambda GT ll_obs AND lambda LE lbeta_obs, nlb)
IF(nlb GT 0) THEN Mask[LLbeta]=CoeffDB
;da Lbeta a Lalpha
LLalpha=where(lambda GT lbeta_obs AND lambda LE lalpha_obs, nla)
IF(nla GT 0) THEN Mask[LLalpha]=CoeffDA
;dopo Lalpha
Lnoigm=where(lambda GT lalpha_obs, nigm)
IF(nigm GT 0) THEN Mask[Lnoigm]=1

;compute area
Delta=lambda-SHIFT(lambda,1)
Delta[0]=lambda[1]-lambda[0]

Integrando=Delta*transm*Mask
NormTrans=TOTAL(Delta*transm)

;compute final coefficient
Integra=TOTAL(Integrando)/NormTrans
IF (Integra NE 0) THEN cIGM=1./Integra ELSE cIGM=1D30

END
