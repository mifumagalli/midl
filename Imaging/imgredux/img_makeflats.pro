;+
;procedure that make flats from flat frames 
;
;instr  --> the instrument used
;name   --> list of file names with flats
;filter --> list of filter for each image
;root   -->  observation date
;bias   --> the median bias from makebias if you want to subtract the
;           mean bias
;MINMAX --> array of min and max value to esclude faint or saturated flats  
;STATUS --> the status structre from ccdproc 
;
;-


pro img_makeflats, name, filter,  bias=bias, status=status, root=root, minmax=minmax, $
                   side=side, path=path, instr=instr


;;find the filters
flat=find_different(filter)
splog, "For side ", side, " found filters: ", flat 
  
;;make flats
for flfra=0, n_elements(flat)-1 do begin

    ;;check if flat done
    sss=where(strcompress(status.filflat,/remove_all) eq flat[flfra],numbfl)
    
    if(numbfl gt 0) then splog, "Found flat ", flat[flfra] else begin
        splog, "Working on filter ", flat[flfra]
        
        ;;find images for this filter
        thisfilter=where(strcompress(filter,/remove_all) eq flat[flfra], nfilt)
        
        ;;check if everything is ok (nfilt MUST be > 0)
        if(nfilt eq 0) then begin
            splog, "I cannot find images with this filter: ", flat[flfra], " Returning..."
            return
        endif
        
        ;;call loop flats to make final image
        img_loopflats, name[thisfilter], finalflat, bias=bias, minmax=minmax, path=path, $
          instr=instr
        ;;save final output
        mwrfits, finalflat, 'proc/'+side+root+"_"+flat[flfra]+"flat.fits", /create
        
        ;;update status structure
        status.filflat[flfra]=flat[flfra]
        mwrfits, status, side+root+"status.str", /create
    endelse
endfor

splog, "All done with flats for side ",side 


  
end
