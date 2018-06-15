;procedure that make a full likelihood ratio from random fields
;calls bayes_lrrandom

;numexp --> total number of experiment
;PRIOR  --> those supported by bayes_lr

PRO bayes_makerandom, numexp, PRIOR=prior


Nexp=fix(numexp/3.)+1

;open catalogue dla1
catname='/home/mikifuma/PROGETTI/DLAIMAGING/FINALE/lumfunct/dla1_clean.dat'
;call bayes_lrrandom
qsopos=[785.,152.]
bayes_lrrandom,  catname, Nexp, LR1,qsopos, PRIOR=prior



;open catalogue dla2
catname='/home/mikifuma/PROGETTI/DLAIMAGING/FINALE/lumfunct/dla2_clean.dat'
;call bayes_lrrandom
qsopos=[1128.,712.]
bayes_lrrandom,  catname, Nexp, LR2, qsopos, PRIOR=prior


;open catalogue dla3
catname='/home/mikifuma/PROGETTI/DLAIMAGING/FINALE/lumfunct/dla3_clean.dat'
;call bayes_lrrandom
qsopos=[125.,674.]
bayes_lrrandom,  catname, Nexp, LR3, qsopos, PRIOR=prior
totLR=[LR1,LR2,LR3]


;replace 0.
zero=where(totLR LE 0.,num)
lgLR=totLR
if(num GT 0) then lgLR[zero]=1D-40 

;makehisto
histLR=HISTOGRAM(ALOG10(lgLR),BINSIZE=1.,locations=x_totlr)
histLR=histLR/TOTAL(histLR*1)
;make cumulative
cumulative, histLR, cumLR, BIN=1.


!P.MULTI=[0,2,1]
plot,  x_totlr, histLR/3., psym=10, xtitle="LR rand", ytitle="Freq." 
plot,  x_totlr, cumLR, psym=10, xtitle="LR rand", ytitle="Freq. Cumul."  
!P.MULTI=0

;dump into structure
str={LRrand:totLR}
outff=string("LR_rand_",prior,".fits")
mwrfits, str, outff, /create





END





