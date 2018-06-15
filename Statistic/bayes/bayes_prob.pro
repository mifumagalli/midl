;computes the probabilities for a set of candidate counterparts



;calls bayes_lr
;R --> estimate of reliability
;Pnond --> out: probability of non detection 
;Puniq --> out: probability of unique detection

PRO bayes_prob, R, Pnond, Puniq 

;compute product PR(1-R)
prod=PRODUCT(1-R)
coeff=R/(1-R)

;compute S
S=TOTAL(coeff*prod)+prod

;compute Pnond 
Pnond=prod/S

;compute Puniq  
Puniq=coeff*prod/S


END
