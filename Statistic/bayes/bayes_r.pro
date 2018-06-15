;compute the raliability using montecarlo approach

;calls bayes_lr
;array --> input data for bayes_lr
;PRIOR --> those supported by bayes_lr

PRO bayes_r, array, Rout, PRIOR=prior

;compute LR for all candidates
bayes_lr, array, LR, PRIOR=prior, /SILENT

;compute Rout with empirical approach
Rout=FLTARR(N_ELEMENTS(array))

;load LR random
IF (PRIOR EQ 'SB3')  THEN str=MRDFITS("/home/mikifuma/idl/midl/bayes/LR_rand_SB3.fits",1)
IF (PRIOR EQ 'OB3')  THEN str=MRDFITS("/home/mikifuma/idl/midl/bayes/LR_rand_OB3.fits",1)
IF (PRIOR EQ 'OB25') THEN str=MRDFITS("/home/mikifuma/idl/midl/bayes/LR_rand_OB25.fits",1)
IF (PRIOR EQ 'OB2')  THEN str=MRDFITS("/home/mikifuma/idl/midl/bayes/LR_rand_OB2.fits",1)
IF (PRIOR EQ 'OB0_3')  THEN str=MRDFITS("/home/mikifuma/idl/midl/bayes/LR_rand_OB0_3.fits",1)
IF (PRIOR EQ 'OB0_25') THEN str=MRDFITS("/home/mikifuma/idl/midl/bayes/LR_rand_OB0_25.fits",1)
IF (PRIOR EQ 'OB0_2')  THEN str=MRDFITS("/home/mikifuma/idl/midl/bayes/LR_rand_OB0_2.fits",1)


Ntotal=N_ELEMENTS(str.LRrand)


FOR i=0, N_ELEMENTS(array)-1 DO BEGIN
;find number LR>LR_i
ii=where(str.LRrand GE LR[i],ngreat)
Rout[i]=1-1D0*ngreat/Ntotal
ENDFOR


END

