;+
;
; Simple script to load the pixel tables in fits format  
;
; binary   works with pure binary table 
;
;-

pro muse_loadpixtab, nametab, outhead=outhead, outstr=outstr, binary=binary 

  if ~keyword_Set(binary) then begin
     ;;load main header
     null=mrdfits(nametab,0,headmain)
     xpix=mrdfits(nametab,1,headxpix)
     ypix=mrdfits(nametab,2,headypix)
     lambda=mrdfits(nametab,3,headlambda)
     data=mrdfits(nametab,4,headdata)
     dq=mrdfits(nametab,5,headdq)
     stat=mrdfits(nametab,6,headstat)
     origin=mrdfits(nametab,7,headorigin)

     ;;groupheader
     outhead={headmain:headmain,headxpix:headxpix,headypix:headypix,headlambda:headlambda,$
              headdata:headdata,headdq:headdq,headstat:headstat,headorigin:headorigin}
     
     
  endif else begin
     
     fxbopen, unit, filename, 1
     fxbreadm, unit, ['xpos','ypos','lambda','data','dq','stat','origin'], xpix, $
               ypix,lambda, data, dq, stat, origin
     fxbclose, unit
     outhead={headmain:""}
 
  endelse 
  


  ;;unpack the origin  
  
  ;;slice 
  mask = 32B+16B+8B+4B+2B+1B    ;6 Bit
  offset = 0
  slice = ishft((ulong(origin)), offset)
  slice = slice and mask

  
  ;;ifu
  mask = 16B+8B+4B+2B+1B        ;5 Bit
  offset = -6
  ifu = ishft((ulong(origin)), offset)
  ifu = ifu and mask


  ;;ypix
  ;mask = 4096B+2048B+1024B+512B+256B+128B+64B+32B+16B+8B+4B+2B+1B ;13 Bit
  ;offset = -11
  ;ypix = ishift(byte(ulong(origin)), offset)
  ;ypix = ypix and mask

  ;;xslice
  ;mask = 64B+32B+16B+8B+4B+2B+1B    ;7 Bit
  ;offset = -24
  ;xslice = ishift(byte(ulong(origin)), offset)
  ;xslice = xslice and mask

  ;;groupdata
  outstr={xpix:xpix,ypix:ypix,lambda:lambda,data:data,dq:dq,stat:stat,origin:origin, slice:slice,ifu:ifu}
  
end

