;;test my ks2d function

loadct, 39


;wow, that really sucks, lets just try a 1d ks


;ok, one-d works, I should check my raw thing and make sure it's good

n=100

d=findgen(101)/99
p1=fltarr(n)
p2=fltarr(n)

nne=50.

for i=0,n-1 do begin
  prob_ks, d(i), nne, prob
  p1(i)=prob
  raw_prob_ks, (sqrt(nne)+0.12+0.11/sqrt(nne))*d(i),prob
  p2(i)=prob


endfor


;yes, my raw_prob_ks is working fine.




;stop
probs=fltarr(n)

;for i=0,n-1 do begin

x1=randomn(seed,100)
x2=randomn(seed,100)
;;###########do a loop around this that plots average probability v
;;              offset.  Go ahead and make the clouds seperated in
;; both x and y

y1=randomn(seed,100)
y2=randomn(seed,100)


plot, x1,y1, psym=2, xrange=[min([x1,x2]),max([x1,x2])], yrange=[min([y1,y2]),max([y1,y2])]
oplot, x2,y2, psym=2,color=250


;probs(i)=ks2d(x1,y1,x2,y2)
prob=ks2d( x1,y1,x2,y2)
;probs(i)=prob

print, prob

;endfor


;plothist, probs, bin=.01



end
