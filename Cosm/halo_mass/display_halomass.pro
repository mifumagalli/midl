;+
;
; Display to the screen a ST mass function, sigma m and m star 
; obtained via shethtormen.pro function. All the assumptions made
; there applies here as well
;
;
;
;
;-


;; Do stuff with the halo mass function

 getshhalomass, 0.,model='DEF',  outmass=m0_0,outhalo=n0_0
 getshhalomass, 0.,model='WMAP5',outmass=m1_0,outhalo=n2_0
 getshhalomass, 0.,model='BOLS', outmass=m2_0,outhalo=n1_0

 getshhalomass, 1.,model='DEF',  outmass=m0_1,outhalo=n0_1
 getshhalomass, 1.,model='WMAP5',outmass=m1_1,outhalo=n2_1
 getshhalomass, 1.,model='BOLS', outmass=m2_1,outhalo=n1_1

 getshhalomass, 2.,model='DEF',  outmass=m0_2,outhalo=n0_2
 getshhalomass, 2.,model='WMAP5',outmass=m1_2,outhalo=n2_2
 getshhalomass, 2.,model='BOLS', outmass=m2_2,outhalo=n1_2

 getshhalomass, 3.,model='DEF',  outmass=m0_3,outhalo=n0_3
 getshhalomass, 3.,model='WMAP5',outmass=m1_3,outhalo=n2_3
 getshhalomass, 3.,model='BOLS', outmass=m2_3,outhalo=n1_3

 getshhalomass, 4.,model='DEF',  outmass=m0_4,outhalo=n0_4
 getshhalomass, 4.,model='WMAP5',outmass=m1_4,outhalo=n2_4
 getshhalomass, 4.,model='BOLS', outmass=m2_4,outhalo=n1_4

 getshhalomass, 5.,model='DEF',  outmass=m0_5,outhalo=n0_5
 getshhalomass, 5.,model='WMAP5',outmass=m1_5,outhalo=n2_5
 getshhalomass, 5.,model='BOLS', outmass=m2_5,outhalo=n1_5

;; display 

!x.style=1
!y.style=1

;;Halo mass function
window, 0


h0=0.72
h1=0.705
h2=0.7

;;m_psopen, 'a.eps', /enc, /lan

plot, m0_0+alog10(h0), n0_0/h0^3, xtitle=Textoidl("log h^{-1} M (M_o)"), $
  ytitle=Textoidl(" n (h^3/Mpc^3/dlog10m)"), /ylog, xrange=[9.5,15], yrange=[10^(-4.5),1]
oplot, m1_0+alog10(h1), n1_0/h1^3, line=1
oplot, m2_0+alog10(h2), n2_0/h2^3, line=2

oplot, m0_1+alog10(h0), n0_1/h0^3, color=200
oplot, m1_1+alog10(h1), n1_1/h1^3, line=1, color=200
oplot, m2_1+alog10(h2), n2_1/h2^3, line=2, color=200

oplot, m0_2+alog10(h0), n0_2/h0^3, color=100
oplot, m1_2+alog10(h1), n1_2/h1^3, line=1, color=100
oplot, m2_2+alog10(h2), n2_2/h2^3, line=2, color=100

oplot, m0_3+alog10(h0), n0_3/h0^3, color=50
oplot, m1_3+alog10(h1), n1_3/h1^3, line=1, color=50
oplot, m2_3+alog10(h2), n2_3/h2^3, line=2, color=50

oplot, m0_4+alog10(h0), n0_4/h0^3, color=25
oplot, m1_4+alog10(h1), n1_4/h1^3, line=1, color=25
oplot, m2_4+alog10(h2), n2_4/h2^3, line=2, color=25

oplot, m0_5+alog10(h0), n0_5/h0^3, color=250
oplot, m1_5+alog10(h1), n1_5/h1^3, line=1, color=250
oplot, m2_5+alog10(h2), n2_5/h2^3, line=2, color=250

;;m_psclose


;;Compute and plot mstar at different redshift

;;WMAP5
z1_0=halomass(14,15,0,prec=1,/WMAP5,emsta=ms1_0)
z1_1=halomass(14,15,1,prec=1,/WMAP5,emsta=ms1_1)
z1_2=halomass(14,15,2,prec=1,/WMAP5,emsta=ms1_2)
z1_3=halomass(14,15,3,prec=1,/WMAP5,emsta=ms1_3)
z1_4=halomass(14,15,4,prec=1,/WMAP5,emsta=ms1_4)
z1_5=halomass(14,15,5,prec=1,/WMAP5,emsta=ms1_5)


;;BOLSH
z3_0=halomass(14,15,0,prec=1,/BOLSH,emsta=ms3_0)
z3_1=halomass(14,15,1,prec=1,/BOLSH,emsta=ms3_1)
z3_2=halomass(14,15,2,prec=1,/BOLSH,emsta=ms3_2)
z3_3=halomass(14,15,3,prec=1,/BOLSH,emsta=ms3_3)
z3_4=halomass(14,15,4,prec=1,/BOLSH,emsta=ms3_4)
z3_5=halomass(14,15,5,prec=1,/BOLSH,emsta=ms3_5)


;;DEF
z0_0=halomass(14,15,0,prec=1,/BOLSH,emsta=ms0_0)
z0_1=halomass(14,15,1,prec=1,/BOLSH,emsta=ms0_1)
z0_2=halomass(14,15,2,prec=1,/BOLSH,emsta=ms0_2)
z0_3=halomass(14,15,3,prec=1,/BOLSH,emsta=ms0_3)
z0_4=halomass(14,15,4,prec=1,/BOLSH,emsta=ms0_4)
z0_5=halomass(14,15,5,prec=1,/BOLSH,emsta=ms0_5)

;;m_psopen, 'b.eps', /enc, /lan


window, 1

plot, [5,4,3,2,1,0], alog10(0.705)+[ms1_5,ms1_4,ms1_3,ms1_2,ms1_1,ms1_0], $
  xtitle='Redshift', ytitle=Textoidl("log h^{-1} M_{star}"), yrange=[6,13]

oplot,  [5,4,3,2,1,0], alog10(0.72)+[ms0_5,ms0_4,ms0_3,ms0_2,ms0_1,ms0_0], line=2, color=50
oplot,  [5,4,3,2,1,0], alog10(0.70)+[ms3_5,ms3_4,ms3_3,ms3_2,ms3_1,ms3_0], line=1, color=250

;;m_psclose


end
