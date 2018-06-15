;+
;
; Simple script to write the pixel tables in fits format  
;
;-

pro muse_writepixtab, outname, head=head, str=str, rawunit=rawunit
  ;;load main header
  mwrfits, null, outname, head.headmain, /cre
  mwrfits, str.xpix, outname, head.headxpix
  mwrfits, str.ypix, outname, head.headypix
  mwrfits, str.lambda, outname, head.headlambda
  mwrfits, str.data, outname, head.headdata
  mwrfits, str.dq, outname, head.headdq
  mwrfits, str.stat, outname, head.headstat
  mwrfits, str.origin, outname, head.headorigin


  ;;fix the stupit header again 
  spawn, 'fparkey xpos '+outname+'[1] EXTNAME'
  spawn, 'fparkey ypos '+outname+'[2] EXTNAME'
  spawn, 'fparkey lambda '+outname+'[3] EXTNAME'

  if keyword_set(rawunit) then begin 
   spawn, 'fparkey pix '+outname+'[2] BUNIT'
   spawn, 'fparkey pix '+outname+'[1] BUNIT'
   spawn, 'fparkey Angstrom '+outname+'[3] BUNIT'
  endif else begin
   spawn, 'fparkey rad '+outname+'[2] BUNIT'
   spawn, 'fparkey rad '+outname+'[1] BUNIT'
   spawn, 'fparkey Angstrom '+outname+'[3] BUNIT'
  endelse
  
  spawn, 'fparkey data '+outname+'[4] EXTNAME'
 
  spawn, 'fparkey dq '+outname+'[5] EXTNAME'
  
  spawn, 'fparkey stat '+outname+'[6] EXTNAME'

  spawn, 'fparkey origin '+outname+'[7] EXTNAME'

end
