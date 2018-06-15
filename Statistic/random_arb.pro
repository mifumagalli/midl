;+
;PURPOSE
;	draw from an arbitrary random distribution to arbitrarily high accuracy 
;	(accuracy set by fine keyword)
;SYNTAX
;	rands=random_arb(seed, d1,...,d8, pdf=pdf, [fine=fine, cumu_pdf=cumu_pdf])
;INPUTS
;	seed: seed used for random number generation, behaves
;		the same as built-in IDL routines
;	d1,...,d8: dimension of output
;	pdf: the probability density function. This *must*
;		be regularly spaced along the abscissa 
;		(a.k.a. x-axis of the histogram) if xpdf is not set
;	xpdf: the x axis of the probaility density function
;	fine: by what factor do you want to interpolate the given pdf
;		(see note below) [default to 1d3]
;		for faster performance, interpolate your pdf before hand
;		and set fine to 1 (see note below)
;	cumu_pdf: this is meant only to allow an increase in speed. if you 
;		pass this total(pdf)/max(pdf) the code will work properly and will
;		not do this calculation.
;NOTE:
;	this code matches an arbitary distribution only in an approximate
;	sense to within the interval defined by in the PDF. As a result there 
;	is a discreteness in the returned distribution this can be seen
;	with the following command.	
;	    plothist, random_arb(seed, 1d6, pdf=findgen(10), fine=1), bin=0.1
;	compared with
;	    plothist, random_arb(seed, 1d6, pdf=findgen(10), fine=1d3), bin=0.1
;
;	For absolute greatest speed do the following (example is for a straight line pdf)
;	pdf=findgen(10)
;	npdf=n_elements(pdf)
;	fine=1d3
;       pdf1=interpol(double(pdf), findgen(npdf), (npdf-1)*findgen(fine)/fine)
;	npdf1=n_elements(pdf1)
;	cumu_pdf1=total(pdf1, /cumu)
;	cumu_pdf1=cumu_pdf1/cumu_pdf1[npdf1-1]
;       xpdf1=interpol(double(xpdf), findgen(npdf), (npdf-1)*findgen(fine)/fine)
;	rands=random_arb(seed, 1d6, pdf=pdf1, xpdf=xpdf1, fine=1, cumu_pdf=cumu_pdf1)
;	
;	If small discreteness effects matter to you, I suggest double checking
;	that the distribution matches the PDF
;Written by R. da Silva, UCSC, 11-22-2010
;-
FUNCTION random_arb, seed, d1, d2, d3, c4, d5, d6, d7, d8, pdf=pdf, xpdf=xpdf, fine=fine,$
	cumu_pdf=cumu_pdf
;check inputs
if keyword_set(pdf) EQ 0 then begin
   print, 'RANDOM_ARB: PDF NOT SET, you must set PDf'
   doc_library, 'random_arb'
   return, -1
endif
;case statement to have same inputs as randomn and randomu
case n_params() of
1:res1=randomu(seed)
2:res1=randomu(seed, d1)
3:res1=randomu(seed, d1, d2)
4:res1=randomu(seed, d1, d2, d3)
5:res1=randomu(seed, d1, d2, d3, d4)
6:res1=randomu(seed, d1, d2, d3, d4, d5)
7:res1=randomu(seed, d1, d2, d3, d4, d5, d6)
8:res1=randomu(seed, d1, d2, d3, d4, d5, d6, d7)
9:res1=randomu(seed, d1, d2, d3, d4, d5, d6, d7, d8)
endcase

;create the interpolated pdf if needed
npdf=n_elementS(pdf)
if keyword_set(xpdf) EQ 0 then xpdf=lindgen(npdf)
if keyword_set(fine) EQ 0 then fine=1d3
if fine GT 1.0 then begin
pdf1=interpol(double(pdf), findgen(npdf), (npdf-1)*findgen(fine)/fine)
xpdf1=interpol(double(xpdf), findgen(npdf), (npdf-1)*findgen(fine)/fine)
endif else begin
pdf1=pdf
xpdf1=xpdf
endelse
npdf2=n_elementS(pdf1)
;calculate the cumlative pdf
if keyword_set(cumu_pdf) EQ 0 then begin
  cumu_pdf=total(pdf1, /cumu)
  cumu_pdf=cumu_pdf/cumu_pdf[npdf2-1]
endif

;plot, cumu_pdf, xpdf1, ps=4
;xx=findgen(1d6)/1d5
;oplot, xx, interpol(xpdf1, cumu_pdf, xx),ps=3
;fine=1d3
;pdf2=interpol(double(pdf), findgen(npdf), (npdf-1)*findgen(fine)/fine)
;xpdf2=interpol(double(xpdf), findgen(npdf), (npdf-1)*findgen(fine)/fine)
;help, xpdf2
;help, total(pdf2, /cumu)
;oplot, total(pdf2, /cumu)/total(pdf2), xpdf2, ps=3, color=fsc_color('red')

;transform uniform random numbers to match given pdf
;STOP
res_pdf=interpol(xpdf1, cumu_pdf, res1)
heap_gc
return, res_pdf
end
