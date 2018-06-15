;procedure that returns a sersic function

PRO sersic, X, P, Funct, Deriv

Funct=P[0]*EXP(-P[1]*((X/P[2])^(1./P[3])-1))

DFDa=EXP(-P[1]*((X/P[2])^(1./P[3])-1)) 
DFDb=EXP(-P[1]*(-1+(X/P[2])^(1./P[3])))*P[0]*(1-(X/P[2])^(1./P[3]))
DFDc=(EXP(-P[1]*(-1+(X/P[2])^(1./P[3])))*P[0]*P[1]*X*(X/P[2])^(-1+1./P[3]))/(P[2]^2*P[3])
DFDd=(EXP(-P[1]*(-1+(X/P[2])^(1./P[3])))*P[0]*P[1]*(X/P[2])^(1./P[3])*ALOG(X/P[2]))/P[3]^2

Deriv=[[DFDa],[DFDb],[DFDc],[DFDd]]

END
