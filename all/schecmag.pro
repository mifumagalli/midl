;procedure that returns a schecter function for appearent magnitude
;for absolute magnite, just invert the sign of m and P[1]

;P[0] phi*
;P[1] m*
;P[2] alpha


PRO schecmag, X, P, Funct, Deriv

Funct=0.921034*P[0]*10^(0.4*(P[2]+1)*(X-P[1]))*EXP(-10^(0.4*(X-P[1])))

DFDa=0.921034*10^(0.4*(P[2]+1)*(X-P[1]))*EXP(-10^(0.4*(X-P[1]))) 
DFDb=0.8483*10^(0.4*(-P[1]+X)+0.4*(1+P[2])*(-P[1]+X))*EXP(-10^(0.4*(-P[1]+X)))*P[0]-0.84830*10^(0.4*(1+P[2])*(-P[1]+X))*EXP(-10^(0.4*(-P[1]+X)))*P[0]*(1+P[2]) 
DFDc=0.848304*10^(0.4*(1+P[2])*(-P[1]+X))*EXP(-10^(0.4*(-P[1]+X)))*P[0]*(-P[1]+X)

Deriv=[[DFDa],[DFDb],[DFDc]]


END
