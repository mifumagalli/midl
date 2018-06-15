
data1 = readfits('lblue0139.fits')
data2 = readfits('lblue0140.fits')
 star=data1[900:1050,*]-median(data1[700:1100,*])
 uband=data2[900:1050,*]-median(data2[700:1100,*])
 st = total(star,1)
 plot, st
 ub = total(uband,1)
 plot, ub
 st = st - min(st)
 ub = ub - min(ub)
xdata = 3650.84  - (findgen(2048)-972)*2.138
result = where(st gt 0)
temp = ub(result)/st(result)
ub1 = ub*0.0
ub1(result) = temp
r3 = where(xdata lt 3000)
ub1(r3) = ub1(r3)*0.0

plot, xdata, ub/st*100.0, yrange=[0,100.0], xrange=[3000,4000]
plot, xdata, ub1*100.0, yrange=[0,100.0], xrange=[2500,5500], xstyle=1
r2 = where(ub1 gt 0.5*max(ub1[200:1200]))
print, ub1[r2[0]]
print, max(ub1)
plot, ub1
FWHM = xdata[r2[0]] - xdata[r2[323]]
cl = (xdata[r2[0]] + xdata[r2[323]])/2.0
print, FWHM, cl, xdata[r2[0]], xdata[r2[323]]

set_plot, 'ps'
device, filename='u_measured.ps'
plot, xdata, ub1*100.0, yrange=[0,100.0], xrange=[2500,5500],xstyle=1, xtitle='wavelength (Angstroms)', ytitle='Transmission (%)'
xyouts , 4000, 80 , 'HPP' ,charsize=2.0
xyouts , 3500, 73 , string(xdata[r2[0]]) ,charsize=2.0
xyouts , 4200, 73 , string(xdata[r2[323]]) ,charsize=2.0

xyouts , 4000, 65 , 'FWHM' ,charsize=2.0
xyouts , 3800, 58 , FWHM ,charsize=2.0

xyouts , 4000, 50 , 'center wavelength' ,charsize=2.0
xyouts , 3800, 43 , cl ,charsize=2.0
device, /close
set_plot, 'X'

forprint, xdata, ub1, textout="ukeck_filter.dat"

end