;+
;  Make a zoom-in ps to identify all the lines in a spectrum. 
;  spect  --> the data 
;  system --> a file that contains the name, z_abs, type
;  (e.g. LLS,DLA), color to plot.
;
;-


pro line_identify, spect, system=system, inflg=inflg, $
                   start=start, dend=dend, delta=delta, save=save, $
                   path=path, miny=miny


  if ~keyword_set(inflg) then inflg=1
  if ~keyword_set(start) then start=3800.
  if ~keyword_set(dend) then dend=9300.
  if ~keyword_set(delta) then delta=300.
  if ~keyword_set(save) then save=spect+"_line.ps"
  if ~keyword_set(path) then path="./"
  if ~keyword_set(miny) then miny=-1.

  ;;load data
  if(inflg eq 1) then begin
     fil_sig=strmid(spect,0,strpos(spect,"_F.fits"))+"_E.fits"
     flx=x_readspec(path+spect,wav=wav,sig=sig,fil_sig=path+fil_sig) 
  endif else flx=x_readspec(path+spect,wav=wav,sig=sig,inflg=inflg)
  
  ;;load tell
  readcol, getenv("MIDL")+"/Spectra/lists/telluric.txt", tel_id, tel_lam, tel_int, tel_fwhm, $
           format='I,D,F,F'


  ;;load abs
  if keyword_set(system) then readcol, system, abs_name, abs_z, abs_type, abs_color, $
                                       format='A,D,A,A'


  ;;open save
  m_psopen, save, /land, xsize=40, ysize=20

  block=fix((dend-start)/delta)
  
  !x.style=1
  !y.style=1
  

  for bb=0, block-1 do begin

     xrange=[start+bb*delta,start+(bb+1)*delta]
     yrange=[miny,max(flx[where(wav gt xrange[0] and wav lt xrange[1])])*1.1]
     
     plot, wav, flx, psym=10, xrange=xrange, /nodata, yrange=yrange

     ;;oplot lines
     for tt=0, n_elements(tel_id)-1 do begin
        
        if(tel_int[tt] lt 0.75) then $
        oplot, [tel_lam[tt],tel_lam[tt]], [-1d10,1d10], line=1, color=fsc_color("tan4") 
     endfor


     ;;oplot absorbers
     if keyword_set(system) then begin

        ;;loop over absorbers
        for aa=0, n_elements(abs_z)-1 do begin
           
           ;;load lines
           if(abs_type[aa] eq 'lls') then  readcol, $
              getenv("XIDL_DIR")+"/Spec/Lines/Lists/lls.lst", $
              sys_lam, sys_line, sys_lwav, sys_stre, format='D,A,A,F' 
        

           ;;plot
           for tt=0, n_elements(sys_lam)-1 do begin
              
              if(sys_stre[tt] gt 0.001) then begin
                 oplot, (1+abs_z[aa])*[sys_lam[tt],sys_lam[tt]], [-1d10,1d10], line=1,$
                        color=fsc_color(abs_color[aa])

                 xyouts, (1+abs_z[aa])*sys_lam[tt], miny, sys_line[tt]+sys_lwav[tt], $
                         Orie=90, ali=0, charsi=1., color=fsc_color(abs_color[aa])
                 
              endif
           endfor
        endfor
     endif


     ;;replot
     plot, wav, flx, psym=10, xrange=xrange, yrange=yrange, /noera
     oplot, wav, sig, psym=10, color=fsc_color('red')

     
  endfor

  m_psclose

  




end
