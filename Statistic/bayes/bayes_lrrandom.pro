;procedure that compute the likelihood ratio from random fields
;calls bayes_lr


;PRIOR    --> those supported by bayes_lr
;catname  --> name of a input catalogue 
;numexp   --> number of an experiment  
;LRout    --> output LR distribution 
;qsopos   --> approx x,y position of the quasar

PRO bayes_lrrandom,  catname, numexp, LRout, qsopos, PRIOR=prior
;check and running!!!!

common colori

;set Nexp per each catalogue
;note that in each frame there are ~150 indpendent realisation 
;(for search radius of 10")
Nexp=numexp


readcol, catname, x, y, mag, errmag, star 

;kill stars
nostar=where(star LT 0.85, nums) 
IF(nums GT 0) THEN BEGIN
x=x[nostar]
y=y[nostar]
ENDIF


;find x and y edges
xmin=MIN(x)
xmax=MAX(x)
ymin=MIN(y)
ymax=MAX(y)

;do not considere the gap now.. maybe later.
PS=0.135
;gap=10./0.135


;generate a grid of random positions
Xrand=make_random(Nexp,xmin,xmax)
Yrand=make_random(Nexp,ymin,ymax)


;ESCLUDERE QUASAR POSITION
near=where(abs(qsopos[0]-Xrand) LT 40. and abs(qsopos[1]-Yrand) LT 40.,qq)
IF(qq GT 0.) THEN BEGIN
Xrand[near]=Xrand[near]+40.
Yrand[near]=Yrand[near]-40.
print, qq, " near position updated!"
ENDIF 


;start loop for each experiment
Search=10./PS
all_LR=0

FOR i=0, Nexp-1 DO BEGIN
;find elements within 10"
bimp=SQRT((Xrand[i]-x)^2+(Yrand[i]-y)^2)
ind=where(bimp LT search,quant)
IF (quant GT 0) THEN BEGIN
;make in arcsec
bimp=bimp[ind]*PS
;convert bimp in kpc using right  distance 
IF (PRIOR EQ 'SB3')  THEN impact_calculator, 3.0, bimp, bproper, /impact, /proper, /DEF
IF (PRIOR EQ 'OB3')  THEN impact_calculator, 3.0, bimp, bproper, /impact, /proper, /DEF
IF (PRIOR EQ 'OB25') THEN impact_calculator, 2.5, bimp, bproper, /impact, /proper, /DEF
IF (PRIOR EQ 'OB2')  THEN impact_calculator, 2.0, bimp, bproper, /impact, /proper, /DEF
IF (PRIOR EQ 'OB0_3')  THEN impact_calculator, 3.,  bimp, bproper, /impact, /proper, /DEF
IF (PRIOR EQ 'OB0_25') THEN impact_calculator, 2.5, bimp, bproper, /impact, /proper, /DEF
IF (PRIOR EQ 'OB0_2')  THEN impact_calculator, 2.,  bimp, bproper, /impact, /proper, /DEF
;get LR (NOTE!!!!!!!!!!!!!!! Only for prior SB3,OB3,OB25, OB2 and OB0_ you can pass array!!!)
bayes_lr, bproper, LR, PRIOR=prior, /silent
;append lr
all_LR=[all_LR,LR]
ENDIF
ENDFOR

;return
mm=N_ELEMENTS(all_LR)
LRout=all_LR[1:mm-1]

END






