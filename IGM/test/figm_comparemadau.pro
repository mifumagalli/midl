;a comparison with madau calculation

PRO figm_comparemadau


common colori

path="/a/miki/PROGETTI/IGM/result/full_evol/"
tail="_c91.dat"

;read my simulation 
figm_readigmcalc, path+"IGM_transmit_z4"+tail, lambda4, trans4, /noplot
figm_readigmcalc, path+"IGM_transmit_z3"+tail, lambda3, trans3, /noplot
figm_readigmcalc, path+"IGM_transmit_z2"+tail, lambda2, trans2, /noplot
figm_readigmcalc, path+"IGM_transmit_z1"+tail, lambda1, trans1, /noplot

figm_readigmcalc, path+"IGM_transmit_z5"+tail, lambda5, trans5, /noplot
figm_readigmcalc, path+"IGM_transmit_z6"+tail, lambda6, trans6, /noplot
figm_readigmcalc, path+"IGM_transmit_z7"+tail, lambda7, trans7, /noplot
figm_readigmcalc, path+"IGM_transmit_z8"+tail, lambda8, trans8, /noplot



;now read  madau
read_madau95, 4., lam_mad4, tau_mad4, igm_mad4
read_madau95, 3., lam_mad3, tau_mad3, igm_mad3
read_madau95, 2., lam_mad2, tau_mad2, igm_mad2
read_madau95, 1., lam_mad1, tau_mad1, igm_mad1


read_madau95, 5., lam_mad5, tau_mad5, igm_mad5
read_madau95, 6., lam_mad6, tau_mad6, igm_mad6
read_madau95, 7., lam_mad7, tau_mad7, igm_mad7
read_madau95, 8., lam_mad8, tau_mad8, igm_mad8


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;from 0-4
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;make a fine grid for integration
prec=0.001
fine_lambda=mkarr(1500,6500,prec)

;interplolate and integrate my function
fine_trans=interpol(trans4,lambda4,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+4) AND fine_lambda LT 1215.67*(1+4),num)
Minte4=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "My integral 4.: ", Minte4

fine_trans=interpol(trans3,lambda3,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+3) AND fine_lambda LT 1215.67*(1+3),num)
Minte3=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "My integral 3.: ", Minte3

fine_trans=interpol(trans2,lambda2,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+2) AND fine_lambda LT 1215.67*(1+2),num)
Minte2=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "My integral 2.: ", Minte2

fine_trans=interpol(trans1,lambda1,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+1) AND fine_lambda LT 1215.67*(1+1),num)
Minte1=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "My integral 1.: ", Minte1


;interplolate and integrate madau
fine_trans=interpol(igm_mad4,lam_mad4,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+4) AND fine_lambda LT 1215.67*(1+4),num)
Pinte4=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "Madau integral 4.: ", Pinte4

fine_trans=interpol(igm_mad3,lam_mad3,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+3) AND fine_lambda LT 1215.67*(1+3),num)
Pinte3=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "Madau integral 3.: ", Pinte3

fine_trans=interpol(igm_mad2,lam_mad2,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+2) AND fine_lambda LT 1215.67*(1+2),num)
Pinte2=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "Madau integral 2.: ", Pinte2

fine_trans=interpol(igm_mad1,lam_mad1,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+1) AND fine_lambda LT 1215.67*(1+1),num)
Pinte1=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "Madau integral 1.: ", Pinte1



;compute ratios
print, "M/P 4: ", Minte4/Pinte4
print, "M/P 3: ", Minte3/Pinte3
print, "M/P 2: ", Minte2/Pinte2
print, "M/P 1: ", Minte1/Pinte1


;some plot here

;m_psopen, "compare_madau_full_evol.eps", /encapsulate, /maxs
window, 0

!y.style=1
!x.style=1

;mia (1500-7000)
plot, lambda4, trans4, xrange=[1500,7000], xtitle=Textoidl("\lambda (A)"),$
         ytitle=Textoidl("exp(-\tau_{eff})")
oplot, lambda3, trans3
oplot, lambda2, trans2
oplot, lambda1, trans1


;madau
oplot, lam_mad4, igm_mad4, color=rosso
oplot, lam_mad3, igm_mad3, color=rosso
oplot, lam_mad2, igm_mad2, color=rosso
oplot, lam_mad1, igm_mad1, color=rosso
;m_psclose


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;from 5-8
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;make a fine grid for integration
prec=0.001
fine_lambda=mkarr(5000,10000,prec)

;interplolate and integrate my function
fine_trans=interpol(trans5,lambda5,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+5) AND fine_lambda LT 1215.67*(1+5),num)
Minte5=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "My integral 5.: ", Minte5

fine_trans=interpol(trans6,lambda6,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+6) AND fine_lambda LT 1215.67*(1+6),num)
Minte6=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "My integral 6.: ", Minte6

fine_trans=interpol(trans7,lambda7,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+7) AND fine_lambda LT 1215.67*(1+7),num)
Minte7=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "My integral 7.: ", Minte7

fine_trans=interpol(trans8,lambda8,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+8) AND fine_lambda LT 1215.67*(1+8),num)
Minte8=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "My integral 8.: ", Minte8


;interplolate and integrate madau
fine_trans=interpol(igm_mad5,lam_mad5,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+5) AND fine_lambda LT 1215.67*(1+5),num)
Pinte5=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "Madau integral 5.: ", Pinte5

fine_trans=interpol(igm_mad6,lam_mad6,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+6) AND fine_lambda LT 1215.67*(1+6),num)
Pinte6=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "Madau integral 6.: ", Pinte6

fine_trans=interpol(igm_mad7,lam_mad7,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+7) AND fine_lambda LT 1215.67*(1+7),num)
Pinte7=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "Madau integral 7.: ", Pinte7

fine_trans=interpol(igm_mad8,lam_mad8,fine_lambda)
interval=where(fine_lambda GT 912.6*(1+8) AND fine_lambda LT 1215.67*(1+8),num)
Pinte8=TOTAL(fine_trans[interval]*prec)/(num*prec)
print, "Madau integral 8.: ", Pinte8


;compute ratios
print, "M/P 5: ", Minte5/Pinte5
print, "M/P 6: ", Minte6/Pinte6
print, "M/P 7: ", Minte7/Pinte7
print, "M/P 8: ", Minte8/Pinte8


;some plot here

;m_psopen, "compare_madau_high_full_evol.eps", /encapsulate, /maxs
window, 1
!y.style=1
!x.style=1

;mia (4500-10000)
plot, lambda5, trans5, xrange=[5000,10000], xtitle=Textoidl("\lambda (A)"),$
         ytitle=Textoidl("exp(-\tau_{eff})")
oplot, lambda6, trans6
oplot, lambda7, trans7
oplot, lambda8, trans8


;madau
oplot, lam_mad5, igm_mad5, color=rosso
oplot, lam_mad6, igm_mad6, color=rosso
oplot, lam_mad7, igm_mad7, color=rosso
oplot, lam_mad8, igm_mad8, color=rosso
;m_psclose


END
