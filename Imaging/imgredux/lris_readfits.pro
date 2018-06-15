;+
; NAME: lris_readfits.pro
;
; PURPOSE: 
;       Given a multi-HDU FITS file, this routine will assemble the
;       various components into a single array based on the header
;       keywords describing the layout.
;
;  NOTE: The blue image is returned as [VidInp1,VidInp2;VidInp3,VidInp4]
;                                      [CCD;CCD1]  
;
;        The red image is returned as [VidInp2,VidInp1;VidInp4,VidInp3]
;                                     [device 1:12; device 1:13]
;  
;
; CALLING SEQUENCE: 
;       array = lris_readfits(filename)
;
; INPUTS:
;	filename = name of the multi-HDU FITS file to read
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       VERBOSE  - if set, give feedback
;       NOTRIM   - if set, data are not trimmed
;       LAYOUT   - set to a named variable that returns new layput
;                  info when notrim is set
;       DISPL    - Open a xatv window                             
;                  AMP: PREDATA[STAR-END]
;                       DATA[STAR-END]
;                       POSTDATA[START-END]
;                       Y[MIN-MAX]
;
;
; OPTIONAL OUTPUT KEYWORD PARAMETERS:
;       HEADER  - retrieve the header from the primary HDU and return
;                 it as a string array
;
; OUTPUTS:
;       This function will return a 2-D floating-point array
;       representing the assembled image.
;
; REQUIREMENTS:
;       - Requires the IDL Astronomy User's Library routines 
;       (http://idlastro.gsfc.nasa.gov/)
;
; EXAMPLES:
;	
; PROCEDURE:
;	- Uses the fits keyword in the header extentions 
;        DETSEC  = '[4096:3073,1:4096]' / NOAO mosaic detector section for ds9
;	to piece the mosiac together.
;
; AUTHOR:
;       Marc Kassis, W. M. Keck Obseravtory
;
; MODIFICATION HISTORY:
;	2009-May-27	MKassis v0.0	Original version
;	2009-Jun-02	GWirth	v0.1	- added secparse routine
;					- now returns type FLOAT array
;					- added optional baseline removal
;	2009-Jun-16	GDW	v0.2	- fix problem with NOTRIM mode
;                                       - use MRDFITS instead of FITSREAD
;       2009-Jun-19     GDW     v0.3    - adapt for binned data
;                                       - write lris_read_amp function
;       2009-Jun-19     GDW     v0.5    - fix bug with PRELINE
;       2009-Jun-24     JMS     v0.6    - fix bug with BZERO
;       2009-Jul-02     GDW     v0.7    - add LINEBIAS option
;                                       - add GAINDATA option
;-      2010-Jan        Mfumagalli      - Extract only a part 
;------------------------------------------------------------------------

;-----------------------------------------------------------------------
function stringify, value, format
;-----------------------------------------------------------------------
;+
; NAME:
;	STRINGIFY
;
; PURPOSE:
;	This function converts a real value into a string based on the 
;	specified format.
;
; CATEGORY:
;	Formatting.
;
; CALLING SEQUENCE:
;	Result = STRINGIFY( Value, Format)
;
; INPUTS:
;	Value	Value which is to be converted.
;
; OPTIONAL INPUTS:
;	Format	Optional format statement
;	
; OUTPUTS:
;	This function returns a string version of the input array with 
;	leading and trailing spaces removed.
;
; RESTRICTIONS:
;	None
;
; EXAMPLE:
;	1) format a number into a string:
;		label = "X=" + stringify(x)
;
;	1) format a number into a string, with formatting:
;		label = "X=" + stringify(x,'(f15.3)')
;
; AUTHOR:
;	Gregory D. Wirth, W. M. Keck Observatory
;
; MODIFICATION HISTORY:
; 	1999-Dec-14	GDW	Original version
;-
;-----------------------------------------------------------------------

; verify input...
np = n_params()
if np lt 1 or np gt 2 then message, 'wrong number of parameters'

; format the value into an appropriate string...
if np gt 1 then begin
   text = string( format=format, value)
endif else begin
   text = string( value)
endelse

; remove leading and trailing spaces...
return, strtrim( text, 2)
end


;------------------------------------------------------------------------
pro secparse, section, x1, x2, y1, y2
;------------------------------------------------------------------------
; Parse a section string of the form [x1:x2,y1:y2], returning
; four values
;------------------------------------------------------------------------
buf = (stregex( section, '\[([0-9]+):([0-9]+),([0-9]+):([0-9]+)\]', $
             /subexp, /extract))[1:4] - 1
x1 = buf[0]
x2 = buf[1]
y1 = buf[2]
y2 = buf[3]
end

;------------------------------------------------------------------------
function lris_read_amp, filename, ext, $
 predata=predata, postdata=postdata, header=header, $
 x1=x1, x2=x2, y1=y1, y2=y2
;------------------------------------------------------------------------
; Read one amp from LRIS mHDU image
;------------------------------------------------------------------------

;; Get the pre and post pix values
;; for LRIS red POSTLINE = 20, POSTPIX = 80, PRELINE = 0, PRECOL = 12
header   = headfits(filename)
precol   = sxpar(header, 'precol')
postpix  = sxpar(header, 'postpix')

;; Deal with binning
BINNING = sxpar(header, 'BINNING')
buf = (stregex( binning, '([0-9]+),([0-9]+)', /subexp, /extract))[1:2]
xbin = buf[0]
ybin = buf[1]
precol = precol/xbin
postpix = postpix/xbin

;; get entire extension...
temp = mrdfits( filename, ext, header, /silent, /unsigned)
tsize = size(temp)
nxt = tsize[1]

;; parse the DETSEC keyword to determine the size of the array.
detsec = sxpar(header, 'DETSEC')
secparse, detsec, x1, x2, y1, y2

;; parse the DATASEC keyword to determine the size of the science region
datasec = sxpar(header, 'DATASEC')
secparse, datasec, xdata1, xdata2, ydata1, ydata2

;; grab the components...
predata  = temp[0:precol-1,*]
data     = temp[xdata1-1:xdata2-1,*]
postdata = temp[nxt-postpix:nxt-1,*]

;; flip in X as needed...
if(x1 gt x2) then begin 
   xt=x2
   x2=x1
   x1=xt
   data = reverse(temporary(data),1)
end

;; flip in Y as needed...
if(y1 gt y2) then begin 
   yt=y2
   y2=y1
   y1=yt
   data  = reverse(temporary(data), 2)
   predata  = reverse(temporary(predata), 2)
   postdata = reverse(temporary(postdata),2)
end

return, data

end

;------------------------------------------------------------------------
function lris_readfits, filename, HEADER=header, $
 NOTRIM=notrim, VERBOSE=verbose, layout=layout, DISPL=displ
;------------------------------------------------------------------------

;; init vars
exten=0                         ; number of extensions
n_ext=0                        ; number of extensions in input image
errmsg=''                       ;  error message set to null
xmax = 0
ymax = 0
xmin = 10000
ymin = 10000

;; verify access to file...
if ~ file_test( filename, /read) then begin
   message, 'Cannot open requested file :' + filename
endif

;; count extensions...
fits_info, filename, /SILENT, N_ext=n_ext

;; allocate array to hold the number of columns in each extension...
xcol = intarr(n_ext)

;; Get the pre and post pix values
;; for LRIS red POSTLINE = 20, POSTPIX = 80, PRELINE = 0, PRECOL = 12
header  = headfits(filename)
PRECOL = sxpar(header, 'PRECOL')
POSTPIX = sxpar(header, 'POSTPIX')
PRELINE = sxpar(header, 'PRELINE')
POSTLINE = sxpar(header, 'POSTLINE')

if keyword_set(verbose) then $
 message, 'PRECOL='     + stringify(PRECOL) $
          + ' POSTPIX='  + stringify(POSTPIX) $
          + ' PRELINE='  + stringify(PRELINE) $
          + ' POSTLINE=' + stringify(POSTLINE), /info


;; get the x and y binning factors...
BINNING = sxpar(header, 'BINNING')
buf = (stregex( binning, '([0-9]+),([0-9]+)', /subexp, /extract))[1:2]
xbin = buf[0]
ybin = buf[1]

;; First read over the header info to determine the size of the output
;; array...
for i=1,n_ext do begin
   theader = headfits(filename, EXTEN=i, ERRMSG=errmsg)
   detsec = sxpar(theader, 'DETSEC')
   if(detsec NE '0') then begin 

      ;;parse the DETSEC keyword to determine the size of the array.
       secparse, detsec, x1, x2, y1, y2

       ;; find the range of detector space occupied by the data 
       ;; [xmin:xmax,ymin:ymax]
       xt = x2 > x1
       xmax = xt > xmax 
       yt = y2 > y1
       ymax = yt > ymax 

       ;; find the min size of the array
       xt = x1 < x2
       xmin = xmin < xt
       yt = y1 < y2
       ymin = ymin < yt

       xcol[i-1] = xt

   endif 
endfor

;; determine the output array size...
nx = xmax - xmin + 1
ny = ymax - ymin + 1

;; change size for binning...
nx = nx / xbin
ny = ny / ybin

;; Update PRECOL and POSTPIX
precol = precol / xbin
postpix = postpix / xbin

;; change size for pre/postscan...
if keyword_set(NOTRIM) then begin
   nx += n_ext*(precol+postpix)  ;; JXP
   ny += preline + postline
endif 

;; allocate output array...
array = fltarr(nx, ny)
if keyword_set(VERBOSE) then $
 message, 'Creating an array of size '+stringify(nx)+ ' by '+stringify(ny), /info

order = sort(xcol)

if keyword_set(NOTRIM) then begin
   ;store layout information
   layout={PREDATA:[0,0],DATA:[0,0],POSTDATA:[0,0],Y:[0,0]}
   layout=replicate(layout,4)
endif

;; insert extensions into master image...
for i=1,n_ext do begin

   ;; grab complete extension...
   data = lris_read_amp( filename, i, $
                         predata=predata, postdata=postdata, $
                         x1=x1, x2=x2, y1=y1, y2=y2)

   ;; insert components into output array...
   if keyword_set(NOTRIM) then begin

 
       ;; insert predata...
       buf = size(predata)
       nxpre = buf[1]
       xs = order[i-1]*PRECOL
       xe = xs + nxpre - 1
     
       ;store info
       layout[i-1].predata[0]=xs
       layout[i-1].predata[1]=xe
       
       if keyword_set(VERBOSE) then begin
           section = '['+stringify(xs)+':'+stringify(xe)+',*]'
           message, 'inserting extension '+stringify(i)+ $
                    ' predata  in '+section, /info
       endif 
       array[xs:xe,*] = predata

       ;; insert data...
       buf = size(data)
       nxdata = buf[1]
       xs = n_ext*precol + (x1-xmin)/xbin
       xe = xs + nxdata - 1
       
       ;store info   
       layout[i-1].data[0]=xs
       layout[i-1].data[1]=xe
  
   
       if keyword_set(VERBOSE) then begin
           section = '['+stringify(xs)+':'+stringify(xe)+',*]'
           message, 'inserting extension '+stringify(i)+ $
                    ' data     in '+section, /info
       endif 
       array[xs:xe,*] = data

       ;; insert postdata...
       buf = size(postdata)
       nxpost = buf[1]
       xs = nx - n_ext*postpix + order[i-1]*postpix
       xe = xs + nxpost - 1

       ;store info
       layout[i-1].postdata[0]=xs
       layout[i-1].postdata[1]=xe
  
 
       
       if keyword_set(VERBOSE) then begin
           section = '['+stringify(xs)+':'+stringify(xe)+',*]'
           message, 'inserting extension '+stringify(i)+ $
                    ' postdata in '+section, /info
       endif 
       array[xs:xe,*] = postdata
       
       
       ;find y info
       buf = size(data)
       nydata = buf[2]
       ys = (y1-ymin)/ybin
       ye = ys + nydata - 1 - postline
       
       layout[i-1].y[0]=ys
       layout[i-1].y[1]=ye
  
   

    endif else begin
       buf = size(data)
       nxdata = buf[1]
       nydata = buf[2]

       xs = (x1-xmin)/xbin
       xe = xs + nxdata - 1
       ys = (y1-ymin)/ybin
       ye = ys + nydata - 1 - postline

       yin1 = PRELINE
       yin2 = nydata - POSTLINE - 1

       if keyword_set(VERBOSE) then begin
           section = '['+stringify(xs)+':'+stringify(xe)+ $
                     ','+stringify(ys)+':'+stringify(ye)+']'
           message, 'inserting extension '+stringify(i)+ $
                    ' data     in '+section, /info
       endif 

       array[xs:xe,ys:ye] = data[*,yin1:yin2]

   endelse
end

;; make sure BZERO is a valid integer for IRAF
OBZERO = sxpar(header, 'BZERO')
sxaddpar, header, 'O_BZERO', obzero
sxaddpar, header, 'BZERO', 32768-obzero

if keyword_set(displ) then xatv, array, /block

return, array

end
