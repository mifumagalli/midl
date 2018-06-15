;+
;
;  Take an image, a flat, info on the gain and exposure time and make
;  a weight image using the inverse variance. 	
;
;  var=cnt * texp + RN^2 + D*t
;  weight=1/var to which apply also mask 
;  ivar=1/var
; 
;-


pro img_makeweight, gimage, flatfits, weight, time=time, instr=instr, side=side,$
                    readnoise=readnoise, date=date, name=name
  
  common ccdinfo, xfull,yfull,nchip,namp,biaslev,satur, goodata

 
  if(instr eq 'ESI') then begin
     readnoise=2.7               ;e-
     darkcurrent=2.              ;e-/pix/hour
  endif
  if(instr eq 'LBC') then begin
     readnoise=10.0              ;e-
     darkcurrent=72.             ;e-/pix/hour
  endif
  if(instr eq 'LRIS') then begin
     if(side eq 'R') then begin
        if(date ge '2013-01-01') then readnoise=4.6 else readnoise=4.0 ;e-
        darkcurrent=0.                                                 ;e-/pix/hour
     endif else begin
        readnoise=3.8           ;e-
        darkcurrent=0.          ;e-/pix/hour
     endelse
  endif
  if(instr eq 'LRISr') then begin
     if(side eq 'R') then begin
        readnoise=6.2  ;e-
        darkcurrent=0. ;e-/pix/hour
     endif else stop 
  endif
  if(instr eq 'IMACS') then begin
     readnoise=3.76            ;e-
     darkcurrent=0.            ;e-/pix/hour
  endif


  ;;construct the variance
  variance=(gimage+readnoise^2.+darkcurrent*time/3600.)
  
  ;;get the inverse variance   
  weight=1./(variance > 1d-99)
  
  ;;zero edges
  zero=where(flatfits lt 0)
  weight[zero]=0.
  
  ;;add bad columns
  if(instr eq 'ESI') then begin
     ;;ESI does not have any before 2012-03-17
     if(strjoin(strsplit(date,'-',/EXTRACT))*1d gt 20120316d) then begin
        
        ;;read bad mask 
        readcol, getenv("MIDL")+"/Imaging/imgredux/utility/esi_badmask.txt", x0, y0, x1, y1
        
        
        for mmm=0, n_elements(x0)-1 do begin
           
           ;;fix border
           if(x0[mmm] lt 0) then x0[mmm]=0
           if(x1[mmm] ge xfull) then x1[mmm]=xfull-1
           if(y0[mmm] lt 0) then y0[mmm]=0
           if(y1[mmm] ge yfull) then y1[mmm]=yfull-1
           
           weight[x0[mmm]:x1[mmm],y0[mmm]:y1[mmm]]=0
           
        endfor
        
     endif 
  endif

  if(instr eq 'LBC') then begin
     
     if(side eq 'R') then bad=mrdfits(getenv("MIDL")+"/Imaging/imgredux/utility/lbc_badmaskr.fits.gz") else $
             bad=mrdfits(getenv("MIDL")+"/Imaging/imgredux/utility/lbc_badmaskb.fits.gz")
     weight=weight*bad
     undefine, bad

  endif
  
  if(instr eq 'LRIS') then begin
     ;;no bad on blue 
     if(side eq 'R') then begin
        
        if(date ge '2013-01-01') then begin
           
           readcol, getenv("MIDL")+"/Imaging/imgredux/utility/lrisr_bad.txt", badcol 
           for mmm=0, n_elements(badcol)-1 do begin
              ;;fix border
              weight[badcol[mmm]:badcol[mmm],*]=0
           endfor
           
        endif else begin
           
           readcol, getenv("MIDL")+"/Imaging/imgredux/utility/lrisr_bad_lbl.txt", x0, y0, x1, y1

           for mmm=0, n_elements(x0)-1 do begin
              
              ;;fix border
              if(x0[mmm] lt 0) then x0[mmm]=0
              if(x1[mmm] ge xfull) then x1[mmm]=xfull-1
              if(y0[mmm] lt 0) then y0[mmm]=0
              if(y1[mmm] ge yfull) then y1[mmm]=yfull-1
              weight[x0[mmm]:x1[mmm],y0[mmm]:y1[mmm]]=0
           endfor

        endelse
     endif
     
  endif

  if(instr eq 'LRISr') then begin
     ;;no bad on blue 
     if(side eq 'R') then begin
        
        readcol, getenv("MIDL")+"/Imaging/imgredux/utility/lrisr_bad_old.txt", x0, y0, x1, y1
        
        for mmm=0, n_elements(x0)-1 do begin
           
           ;;fix border
           if(x0[mmm] lt 0) then x0[mmm]=0
           if(x1[mmm] ge xfull) then x1[mmm]=xfull-1
           if(y0[mmm] lt 0) then y0[mmm]=0
           if(y1[mmm] ge yfull) then y1[mmm]=yfull-1
           weight[x0[mmm]:x1[mmm],y0[mmm]:y1[mmm]]=0
        endfor
        
     endif else stop 
  endif
  if(instr eq 'IMACS') then begin
     ;;read bad mask 
     readcol, getenv("MIDL")+"/Imaging/imgredux/utility/imacs_badmask.txt", x0, y0, x1, y1
     
     for mmm=0, n_elements(x0)-1 do begin
        ;;fix border
        if(x0[mmm] lt 0) then x0[mmm]=0
        if(x1[mmm] ge xfull) then x1[mmm]=xfull-1
        if(y0[mmm] lt 0) then y0[mmm]=0
        if(y1[mmm] ge yfull) then y1[mmm]=yfull-1
        
        weight[x0[mmm]:x1[mmm],y0[mmm]:y1[mmm]]=0
     endfor
  endif

  ;;apply manual masking if needed 
  cckfil=file_test('mask/'+strmid(name,0,strpos(name,'.fits'))+'.reg')
   if(cckfil gt 0) then mask_region, 'mask/'+strmid(name,0,strpos(name,'.fits'))+'.reg', weight
   
end





