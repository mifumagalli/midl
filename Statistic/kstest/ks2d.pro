;+
;
;
;
;
;Stolen from http://www.astro.washington.edu/users/yoachim/code.php
;
;compute the pobability that two 2-d arrays are drawn from the same
;distribution as in numerical recipies.  
;
;
;-


function ks2d, x1,y1,x2,y2
;;x1,y1 points for model 1
;;x2,y2 points for model 2

da1=0
;;loop through the first array as center points

n1=n_elements(x1)
n2=n_elements(x2)
for i=0,n1-1 do begin
  d1=quad_fracs(x1(i),y1(i),x1,y1)
  d2=quad_fracs(x1(i),y1(i),x2,y2)
  da1=max([da1,max(abs(d1-d2))])
endfor

da2=0
;;loop throught the second array
for i=0,n2-1 do begin
  d1=quad_fracs(x2(i),y2(i),x1,y1)
  d2=quad_fracs(x2(i),y2(i),x2,y2)
  da2=max([da2,max(abs(d1-d2))])
endfor

;;average the d's
d=mean([da1,da2])

n=(n1*n2)/(n1+n2)

;;get linear correlation coef for each sample
r1=correlate(x1,y1)
r2=correlate(x2,y2)
rr=sqrt(1.-0.5*(r1^2+r2^2))

;;ESTIMATE probability using raw K-S prob
;;this line was ain error, I coppied teh wrong equation from NR.
;;Caught by Stephane Blondin
;;lambda= ( sqrt(n)*d)/(1.+sqrt(1-rr^2)*(0.25-.75/sqrt(n)))

;;new version of Numerical Recipies says:
lambda=(sqrt(n)*d)/(1.+rr*(0.25-.75/sqrt(n)))

;;stop
raw_prob_ks, lambda, prob

return, prob
end



