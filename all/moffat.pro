;procedure that returns a moffat function

PRO moffat, X, P, Funct, Deriv

Funct=P[0]/(1+(X/P[1])^2)^P[2]

DFDa=1/(1+(X/P[1])^2)^P[2] 
DFDb=(2*P[0]*P[2]*X^2.)/(P[1]^3*((1+(X/P[1])^2))^(1+P[2]))
DFDc=-P[0]*(1+(X/P[1])^2)^(-P[2])*ALOG(1+(X/P[1])^2)

Deriv=[[DFDa],[DFDb],[DFDc]]

END
