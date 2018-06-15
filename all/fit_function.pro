;procedura che fitta una coppia di vettori X and Y
;restitisce parametri fit
;In input param e' un guess dei valori, in output e' il fit

;sersic-> sersic radial profile
;schecmag-> schecter function for appearant magnitude  
;moffat-> moffat profile
;gauss -> gaussian
;sers_conv -> sersic convolved with a moffat (fix parameter for moffat)

PRO fit_function, X, Y, PARAM, ERR, FIT, CHI, WEIGHTS=Weights,SERSIC=sersic, SCHEMAG=schemag,$
    MOFFAT=moffat, GAUSS=gauss, SERS_CONV=sers_conv

IF KEYWORD_SET(sersic) THEN Fit=CURVEFIT(X,Y,Weights,PARAM,CHISQ=chi,FUNCTION_NAME="sersic",YERROR=err,ITMAX=100) 
IF KEYWORD_SET(schecmag) THEN Fit=CURVEFIT(X,Y,Weights,PARAM,CHISQ=chi,FUNCTION_NAME="schecmag",YERROR=err,ITMAX=100)
IF KEYWORD_SET(moffat) THEN Fit=CURVEFIT(X,Y,Weights,PARAM,CHISQ=chi,FUNCTION_NAME="moffat",YERROR=err,ITMAX=100)
IF KEYWORD_SET(gauss) THEN Fit=GAUSSFIT(X,Y,PARAM,CHISQ=chi,NTERMS=3,SIGMA=ERR) 
IF KEYWORD_SET(sers_conv) THEN Fit=CURVEFIT(X,Y,Weights,PARAM,CHISQ=chi,FUNCTION_NAME="sers_conv",YERROR=err,ITMAX=500,/NODERIVATIVE)

END
