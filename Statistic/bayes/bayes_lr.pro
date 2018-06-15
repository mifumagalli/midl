;procedure that compute the likelihood ratio

; Several PRIOR available (now only b)
; SB3    bimpact only from simulation z=3
; OB3    bimpact from observation, scaled to z=3
; OB25   bimpact from observation, scaled to z=2.5
; OB2    bimpact from observation, scaled to z=2
; OB0_3    bimpact from observation, NON evolution, applied to z=3
; OB0_25   bimpact from observation, NON evolution, applied to z=2.5
; OB0_2    bimpact from observation, NON evolution, applied to z=2


;LR-->out value of likelihood ratio
;arrayin---> B impact in kpc for SB3,OB3, OB25,OB2 and OB0_


PRO  bayes_lr, arrayin, LR, PRIOR=prior, SILENT=SILENT
;checked and working!!

IF (PRIOR EQ 'SB3') THEN BEGIN
IF ~keyword_set(SILENT) THEN print, "Using prior SB3"

;problem with extrapolation... switch to analytical....
;;read f(b)
;readcol, "sb3_table.dat", x, y, skipline=2, /silent
;;suppress extrapolation
;inside=where(arrayin LT MAX(x) and arrayin GT MIN(x),quan,complement=out,ncomplement=nout)
;IF(quan GT 0) THEN extract_interpol, x, y, arrayin[inside], fb, PREC=0.01
;;set LR
;LR=FLTARR(N_ELEMENTS(arrayin))
;IF(quan GT 0) THEN LR[inside]=Fb
;IF(nout GT 0) THEN LR[out]=0.

IF ~keyword_set(SILENT) THEN print, "Using prior SB3"
;set parameters [A,alp,B,beta] and compute f(b)
P=[0.234,0.68,0.37,1.05]
Fb=(P[0]*arrayin^P[1])*EXP(-P[2]*arrayin^P[3])
;set LR
LR=Fb



ENDIF
  
 
IF (PRIOR EQ 'OB3') THEN BEGIN
IF ~keyword_set(SILENT) THEN print, "Using prior OB3"
;set parameters [A,alp,B,beta] and compute f(b)
P=[0.166,0.37,0.14,1.29]
Fb=(P[0]*arrayin^P[1])*EXP(-P[2]*arrayin^P[3])
;set LR
LR=Fb
ENDIF

IF (PRIOR EQ 'OB25') THEN BEGIN
IF ~keyword_set(SILENT) THEN print, "Using prior OB25"
;set parameters and compute f(b)
P=[0.138,0.37,0.12,1.29]
Fb=(P[0]*arrayin^P[1])*EXP(-P[2]*arrayin^P[3])
;set LR
LR=Fb
ENDIF

IF (PRIOR EQ 'OB2') THEN BEGIN
IF ~keyword_set(SILENT) THEN print, "Using prior OB2"
;set parameters and compute f(b)
P=[0.112,0.37,0.10,1.29]
Fb=(P[0]*arrayin^P[1])*EXP(-P[2]*arrayin^P[3])
;set LR
LR=Fb
ENDIF


IF (PRIOR EQ 'OB0_3' OR PRIOR EQ 'OB0_25' OR PRIOR EQ 'OB0_2') THEN BEGIN
IF ~keyword_set(SILENT) THEN print, "Using prior", prior
;set parameters and compute f(b)
P=[0.064,0.37,0.057,1.29]
Fb=(P[0]*arrayin^P[1])*EXP(-P[2]*arrayin^P[3])
;set LR
LR=Fb
ENDIF


END



