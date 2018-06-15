;;;test procedure to look at the random LR stored in structures


common colori

sim3=mrdfits("LR_rand_SB3.fits",1)
obs3=mrdfits("LR_rand_OB3.fits",1)
obs25=mrdfits("LR_rand_OB25.fits",1)
obs2=mrdfits("LR_rand_OB2.fits",1)
ne3=mrdfits("LR_rand_OB0_3.fits",1)
ne25=mrdfits("LR_rand_OB0_25.fits",1)
ne2=mrdfits("LR_rand_OB0_2.fits",1)

help, sim3.lrrand
help, obs3.lrrand
help, obs25.lrrand
help, obs2.lrrand
help, ne3.lrrand
help, ne25.lrrand
help, ne2.lrrand


s=histogram(ALOG10(sim3.lrrand[where(sim3.lrrand GT 0.)]), bin=.5,locations=sbas)
s1=histogram(ALOG10(obs3.lrrand[where(obs3.lrrand GT 0.)]), bin=.5,locations=sbas1)
s2=histogram(ALOG10(obs25.lrrand[where(obs25.lrrand GT 0.)]), bin=.5,locations=sbas2)
s3=histogram(ALOG10(obs2.lrrand[where(obs2.lrrand GT 0.)]), bin=.5,locations=sbas3)
s4=histogram(ALOG10(ne3.lrrand[where(ne3.lrrand GT 0.)]), bin=.5,locations=sbas4)
s5=histogram(ALOG10(ne25.lrrand[where(ne25.lrrand GT 0.)]), bin=.5,locations=sbas5)
s6=histogram(ALOG10(ne2.lrrand[where(ne2.lrrand GT 0.)]), bin=.5,locations=sbas6)



plot,  sbas,   s, psym=10, xrange=[-15,0], yrange=[0,8000], xtitle="Log LR rand", ytitle="Freq"
oplot, sbas1, s1, psym=10, color=rosso
oplot, sbas2, s2, psym=10, color=blu
oplot, sbas3, s3, psym=10, color=verde
oplot, sbas4, s4, psym=10, color=rosso, line=1
oplot, sbas5, s5, psym=10, color=blu  , line=1
oplot, sbas6, s6, psym=10, color=verde, line=1

END
