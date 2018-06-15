;+
;
;
;   Once wcs are in the header, this pro does a quick combine
;   without dealing with the reprojection of images. Works fine
;   for a quick look or for very stable instruments like ESI.
;
;
;
;
;
;-



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;A couple of work functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;take an image and wgt plus header and does a very simple source
;;extractions, returning position in ra and dec
pro img_quickstack_match, image, wgt, header, ra_img, dec_img, mag_img, instr=instr
                      

;;find pixel size
extast, header, ast
pxsize=ast.cd[1,1]*3600.

;;select non zero weight
nozero=where(wgt gt 0)

;;get a sky estimate to select 5 sigma
mmm, image[nozero], skymod, skysig, skysk

;;clip the edges using weight map
mmm, wgt[nozero], wgtmod, wgtsig, wgtsk
setzero=where(wgt lt wgtmod-wgtsig,xzer)

imgnozr=image
if(xzer gt 0) then imgnozr[setzero]=0.

if (instr eq 'LBC') then nsig=8 else nsig=5
find, imgnozr,  x, y, flux, sharp, round, nsig*skysig, 1./pxsize, [-1,1], [0.2,1.0] 


mag_img=-2.5*alog10(flux)

;;compute the wcs 
xy2ad, x, y, ast, ra_img, dec_img

end

;;take an list of ra,dec in two images and plot
;;the comparisons
pro img_quickstack_printcfr, ra1, dec1, ra2, dec2, mag2, stringtitle

;;match catalogues
spherematch, ra1, dec1, ra2, dec2, 1/3600., match1, match2, distance, maxmatch=1

;;find median distance
djs_iterstat, distance*3600., median=meddist

;;fits plot space for match withing 1 arcsec
plot, ra1[match1], dec1[match1], xtitle="RA (deg)", ytitle="DEC (deg)", psym=symcat(6), $
  title=stringtitle, thick=4
oplot, ra2[match2], dec2[match2], psym=symcat(9), symsize=1.5, thick=4, color=250

;;now plot distance error with respect to magnitude
plot, mag2[match2], distance*3600., psym=symcat(6), $
  title="Median WCS offset "+rstring(meddist), $
  xtitle="mag",  ytitle="Distance Offset (arcsec)"
oplot, [-100,100], [meddist,meddist], line=2

;;now plot distance error with respect to ditence from the center
ra_center=median(ra2)
dec_center=median(dec2)
distcenter=sqrt((ra2[match2]-ra_center)^2+(dec2[match2]-dec_center)^2)*60.

plot, distcenter, distance*3600., psym=symcat(6), $
  title="Median WCS offset "+rstring(meddist), $
  xtitle="Distance from center FOV (arcmin)",  ytitle="Distance Offset (arcsec)"
oplot, [-100,100], [meddist,meddist], line=2


end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Main procedure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro img_quickstack, listfile, outname, path=path,  $
                    sdss=sdss, instr=instr, nocfr=nocfr, refer=refer


splog, 'Starting img_checkwcs at ', systime(), filename='splog.checkwcs', /append


;;set default
if ~keyword_set(path) then path='./Raw/swarp/output/'
if ~keyword_set(instr) then instr='LRIS'

;;load the info from the listfile
readcol, listfile, filename, filt, format='A,A'

;;Find out how many objects
n_obj=n_elements(filename)

if ~keyword_set(nocfr) then begin
    ;;open diagnostic plots
    !x.style=1
    !y.style=1
    m_psopen, "check_wcs_"+filename[0]+".ps", /lan
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Assume images are already projected to a common grid, 
;;so a shift and maybe a rotation is all that this does
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

splog, "Register images"


;;set the reference image
if keyword_set (refer) then reference=refer else reference=filename[0]


;;first load reference image
ref_fit=mrdfits(path+'sci_'+reference,1,ref_header,/sil)
ref_wgh=mrdfits(path+'wgt_'+reference,1,/sil)
ref_fit=temporary(ref_fit*ref_wgh)


;;find the position of reference sources
if ~keyword_set(nocfr) then $
  img_quickstack_match, ref_fit, ref_wgh, ref_header, ra_ref, dec_ref, magref, instr=instr

;;now do the stack after alignemt and compare the relative 
;;position after images are registered is done
for i=0, n_obj-1 do begin

    ;;load images
    tmp=mrdfits(path+'sci_'+filename[i],0,main_img_header,/sil)
    sci=mrdfits(path+'sci_'+filename[i],1,img_header,/sil)
    wgt=mrdfits(path+'wgt_'+filename[i],1,/sil)
    
    ;;match in wcs    
    hastrom, sci, img_header, sh_sci,   newhd, ref_header, missing=0.
    hastrom, wgt, img_header, sh_wgt,   newhd, ref_header, missing=0.
    
    ;;update stack
    if(i eq 0) then begin
        totweight=sh_wgt
        stack=sh_sci*sh_wgt
    endif else begin
        totweight+=sh_wgt
        stack+=sh_sci*sh_wgt
    endelse
    
    if ~keyword_set(nocfr) then begin
        ;;find the position of other sources
        img_quickstack_match, sh_sci, sh_wgt, newhd, ra_con, dec_con, mag_con, instr=instr
        
        ;;compare with reference
        img_quickstack_printcfr, ra_ref, dec_ref, ra_con, dec_con, mag_con,$
          "Compare "+reference+" with "+filename[i]
    endif
    
endfor


;;make final image
nozero=where(totweight gt 0)
final=stack-stack
final[nozero]=stack[nozero]/totweight[nozero]

;;write final image
mkhdr, hh,  final
mwrfits, UNDEF, outname, hh, /create
mwrfits, final, outname, ref_header
mwrfits, totweight, outname

if keyword_set(nocfr) then return

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Find all the stars in the white image and make
;;a comparison with USNO or SDSS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;find pixel size
extast, ref_header, ast
pxsize=ast.cd[1,1]*3600.
splog, 'Pixel size ', pxsize

;;find size of field and get object catalogue
nx=fxpar(ref_header,"NAXIS1")
ny=fxpar(ref_header,"NAXIS2")
size = nx > ny 
field_rad = 0.5*size*pxsize/60.
splog, "FOV size arcmin ", 2*field_rad

;;find center 
xy2ad, 0.5*nx, 0.5*ny, ast, racen, decen

;;query
if keyword_set(sdss) then cat='II/294' else cat='USNO-B1'
splog, 'Query to ', cat
info=m_queryvizier(cat,[racen,decen],field_rad)
splog, "Found ", n_elements(info), " objects"

;;extract magnitude, ra and dec
if keyword_set(sdss) then begin
    raj2000=info.raj2000
    dej2000=info.dej2000
    rmag=info.rmag
endif else begin
    rmag=info.R1MAG
    raj2000=info.raj2000
    dej2000=info.dej2000
endelse


;;plot comparison absolute/white
img_quickstack_match, final, totweight, ref_header, ra_whi, dec_whi, mag_whi, instr=instr

img_quickstack_printcfr, ra_whi, dec_whi, raj2000, dej2000, rmag,$
      "Absolute WCS in registered image"


;;close plot
m_psclose



splog, 'All done at ', systime(), /close
    
    
end



