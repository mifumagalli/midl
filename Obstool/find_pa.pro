;+
;compute pa and draw finding chart
; 
;
;ref_pos     --> position in hh:mm:ss and dd:mm:ss or x,y of the reference object 
;tar1_pos    --> position in hh:mm:ss and dd:mm:ss of the main target
;[tar2_pos]  --> position in hh:mm:ss and dd:mm:ss
;imagename   --> file name of the fits file of the field, to produce image finding chart
;name        --> root for the file name to save
;[nsg]       --> number of sigma within which the image is plotted
;             (regolate scale)
;[zoom]      --> The number of pixel in the PA fndchatrt
;ext         --> extension to open 
;
;shifts are computed from reference to target1
;PA is compute from target1 to target2
;
;
;-


pro find_pa, ref_pos, tar1_pos, tar2_pos, imgname=imgname, name=name, nsg=nsg, zoom=zoom, ext=ext, smooth=smooth, sdss=sdss

if ~keyword_set(nsg) then nsg=2.
if ~keyword_set(zoom) then zoom=200.
if ~keyword_set(imgname) then imgname='tmp.fits'
if ~keyword_set(tar2_pos) then begin
    tar2_pos=tar1_pos
    tar1_pos=ref_pos
endif
if ~keyword_set(name) then name='fc_target'
if ~keyword_set(ext) then ext=0.


;;get deg values
x_radec, ref_pos[0], ref_pos[1], ref_rag, ref_deg
x_radec, tar1_pos[0], tar1_pos[1], tar1_rag, tar1_deg
x_radec, tar2_pos[0], tar2_pos[1], tar2_rag, tar2_deg


;;find shift reference to target
observing_night_geometry, tar1_rag, tar1_deg, ref_rag, ref_deg, shift=shift, /silent
splog, "Shift from REF to TARGET (ra positive E, dec positve N): ", shift


;;find PA
observing_night_geometry, tar2_rag, tar2_deg, tar1_rag, tar1_deg, pa=pa, /silent
splog, "PA (N-E): ", pa


;;make finding chart of the object
x_fndchrt, [name+'_fnd',ref_pos[0],ref_pos[1]], imsize=6., /radec,$
  twocirc=[-shift[0]/60.,shift[1]/60.], sdss=sdss


;;plot the slit
if keyword_set(imgname) then begin
   

   
   ;;open fits
   img=mrdfits(imgname,ext,header)
      
   ;;now rotate to have NE orientation
   getrot, header, Rot, CDelt, /SILENT

   adxy, header,  tar2_rag, tar2_deg, tar2_x, tar2_y 
   hrot, img, header, -1, -1, Rot, tar2_x, tar2_y, 0
   
   ;;Extract a subimage from an array and update astrometry in FITS header
   adxy, header,  ref_rag, ref_deg, ref_x, ref_y
   adxy, header,  tar2_rag, tar2_deg, tar2_x, tar2_y 
   adxy, header,  tar1_rag, tar1_deg, tar1_x, tar1_y 
   
   ;;find new corner
   x0=fix(min([tar1_x,tar2_x,ref_x])-zoom) > 0
   x1=fix(zoom+max([tar1_x,tar2_x,ref_x])) < n_elements(img[*,0])-1
   y0=fix(min([tar1_y,tar2_y,ref_y])-zoom) > 0
   y1=fix(zoom+max([tar1_y,tar2_y,ref_y])) < n_elements(img[0,*])-1
   
   hextract, img, header, img, header, x0, x1, y0, y1,/silent
   
   
   ;;get updates sizes
   ysize=n_elements(img[0,*])
   xsize=n_elements(img[*,0])
   
   ;;go from wcs to image (after rotation)
   adxy, header,  ref_rag, ref_deg, ref_x, ref_y
   adxy, header,  tar2_rag, tar2_deg, tar2_x, tar2_y 
   adxy, header,  tar1_rag, tar1_deg, tar1_x, tar1_y 
   
   
   ;;smooth
   if keyword_set(smooth) then img=filter_image(img,smooth=smooth)

   
   ;;psfile
   psfil=strtrim(name,2)+'_pa.ps'
   ;;must use x psopen
   x_psopen, psfil, /portrait
   
   loadct, 0
   
   ;;fix nan
   good = where( finite(img),complement=bed,ncomplement=nbed) 
   if(nbed gt 0) then img[bed]=0.
   med = median(img)
   img=img-med
   mx = max(img, min=mn)
   sig = stddev(img)
   pltmax = (+(nsg*sig)) < mx
   pltmin = (-(nsg*sig)) > mn
   print, med, sig, pltmin, pltmax, mx, mn
   
   imdisp, bytscl(img, min=pltmin, max=pltmax, /nan),$
           pos=[0.108,0.11,0.905,0.778], xstyle=1,ystyle=1, $
           xmargin=[0,0], ymargin=[0,0], charsize=1.5, /NEGATIVE, $
           /axis, XTickformat='(A1)', YTickformat='(A1)', $
           XTICKLEN=1D-5,YTICKLEN=1D-5, title='North', ytitle='East'
   
   
   ;;Label
   loadct, 39
   ;common colori
   
   xyouts, 0.5, 0.95, name, alignment=0.5, charsize=3., $
           color=0, /normal
   
   xsh=strmid(strtrim(shift[0],2),0,8)
   ysh=strmid(strtrim(shift[1],2),0,8)
   
   xyouts, 0.5, 0.9, string("Off(ra): ",xsh,"''  Off(dec): ",ysh,"''"),$
           alignment=0.5, charsize=2.5, color=0, /normal
   xyouts, 0.5, 0.87, string("PA: ", strtrim(pa,2)), $
           alignment=0.5, charsize=2., color=0, /normal
   
   ;;Circle
   x_oplotcirc, 10., x0=ref_x, y0=ref_y, color=fsc_color('green')
   x_oplotcirc, 10., x0=tar1_x, y0=tar1_y, color=fsc_color('blue')
   x_oplotcirc, 10., x0=tar2_x, y0=tar2_y, color=fsc_color('red')
   
   
   ;;oplot 1" slit
   
   ;;get PS
   ps=fxpar(header,'CD2_2')*3600.
   
   ;;draw slit
   m=(tar1_y-tar2_y)/(tar1_x-tar2_x)
   ascissa=mkarr(0,xsize,1.)
   ordinata=m*(ascissa-tar1_x)+tar1_y
   
   ;;oplot slit
   oplot, ascissa, ordinata, thick=1.8, line=2, color=255 
   oplot, ascissa+0.5/PS, ordinata, thick=1.8, color=255 
   oplot, ascissa-0.5/PS, ordinata, thick=1.8, color=255 
   
   
   
   ;;use x psopen
   x_psclose
   
endif

spawn, 'rm -f tmp.fits'


end



