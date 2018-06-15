;+
;PURPOSE
;	in a multiple extension file replace a single extension
;SYNTAX
;	wr_fits_1exten, filename, extension, data [,header=header, out=out]
;INPUTS
;	filename: the name of the file
;	extension: the extenstion to overwrite
;	data: the data to overwrite
;KEYWORDS
;	header: if you want to change the header then set it to this keyword
;	out: if you want to write to a different file then use this
;	clobber: if set then will write over file of the same fileanme
;NOTES
;	writes to a file called 'new_'+filename
;
;Written by R. da Silva, UCSC, 9-24-09
;-

pro wr_fits_1exten, filename, extension, data, header=header, $
	out=out, clobber=clobber
if NOT keyword_set(out) then outfil='new_'+filename else outfil=out

datai=mrdfits(filename, 0, h, /silent)
if extension EQ 0 then begin 
    if keyword_set(header) then mwrfits, data, outfil, header, /create $
    else mwrfits, data, outfil, h, /create
endif else begin
    mwrfits, datai, outfil, h, /create
endelse


exist=1
exten=1
while exist do begin
    datai=mrdfits(filename,exten, h, /silent)
    exist=keyword_set(datai)
    if exist then begin
        if exten NE extension then mwrfits, datai, outfil, h else begin
            if keyword_set(header) then mwrfits, data, outfil, header else $
              mwrfits, data, outfil, h
        endelse
        
    endif
    exten++
endwhile

if keyword_set(clobber) then spawn, 'mv '+'new_'+filename+' '+filename

end
