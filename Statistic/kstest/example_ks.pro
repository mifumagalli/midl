;;example of using ks2d.pro


;try out 2 identical distributions

x1=randomn(seed,500)
y1=randomn(seed,500)

x2=randomn(seed2,2000)
y2=randomn(seed2,2000)

print, 'probability drawn from same parent distribution=',$
     ks2d(x1,y1,x2,y2)
print, '(should be close to 1, or at least not near zero)'

print,''


;make a subtle shift in one of the distributions

print, 'probability drawn from same parent distribution=',$
     ks2d(x1+.05,y1-.05,x2,y2)
print, '(should be a little less than last value)'

print,''


plot,  x2, y2, psym=2

oplot, x1, y1, psym=1, color=250


;now see how they compare when distributions have different spreads
x2=randomn(seed3,2000)*2.
y2=randomn(seed3,2000)*2.

print, 'probability drawn from same parent distribution=',$
     ks2d(x1,y1,x2,y2)
print, '(should be much less than 1)'



end

