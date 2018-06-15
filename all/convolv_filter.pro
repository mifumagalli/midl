;for a given SED and filter transmission in AA, compute the observed luminosity

PRO convolv_filter, Lambda, Flux, Lfilter, Transm, Luminosity 

;compute bin size
Delta=Lfilter-SHIFT(Lfilter,1)
Delta[0]=Lfilter[1]-Lfilter[0]
Nmax=N_ELEMENTS(Lfilter)-1
Delta[0]=Lfilter[1]-Lfilter[0]

;normalizza filtro
TransNorm=TOTAL(Delta*Transm)

;interpola SED su filtro
FSED=INTERPOL(Flux,Lambda,Lfilter)

;convolvi
Luminosity=TOTAL(FSED*Transm*Delta)/TransNorm
END




