

readcol, "U_LRIS.res", lamKeck, TraKeck
readcol, "U_SDSS.res", lamSDSS, TraSDSS

psopen, "Keck_SDSS.ps"

!p.multi=[0,1,2]

plot, lamKeck, TraKeck, xrange=[2800,4200], ytitle="Trans", xtitle="Lambda (AA)"
oplot, lamSDSS, TraSDSS, line=1

;difference 
TraInterSDSS=interpol(TraSDSS,lamSDSS,lamKeck) 
plot, lamKeck, TraInterSDSS/TraKeck, xrange=[2800,4200], ytitle="SDSS/Keck", xtitle="Lambda (AA)"
 
!p.multi=0

psclose


;normalise transmissions

DLkeck=lamKeck-shift(lamKeck,1)
DLkeck[0]=lamkeck[1]-lamkeck[0]
TraKeckNorm=TraKeck/TOTAL(TraKeck*DLkeck)

DLSDSS=lamSDSS-shift(lamSDSS,1)
DLSDSS[0]=DLSDSS[1]-DLSDSS[0]
TraSDSSNorm=TraSDSS/TOTAL(TraSDSS*DLSDSS)

;compute shift in magnitude
keckSDSS=-2.5*ALOG10(TOTAL(TraKeckNorm*DLKeck/lamKeck^2)/TOTAL(TraSDSSNorm*DLSDSS/lamSDSS^2))
print, "Mag Keck-SDSS ", keckSDSS


end
