
   
PRO  man_gain


fits=xmrdfits("lblue0097ff.fits", 0, header, /fscale, /silent)
SXADDPAR, header, "COMMENT", "GAIN APPLIED [1.55,1.55,1.64,1.68]"
fits[3277:4619,*]=1.71*fits[3277:4619,*]
fits[2252:3276,*]=1.62*fits[2252:3276,*]
fits[0:2251,*]=1.55*fits[0:2251,*]
mwrfits, fits,  "lblue0097ff.fits", header, /create

end
