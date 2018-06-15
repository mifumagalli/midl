;+
;
;This procedure makes a median bias from different bias frames
;
;
;name          -->  file list of bias to process
;root          --> name to append to names
;side          --> side of the instrument
;path          --> where data are
;biasmedian    --> in output the median bias
;instr         --> the instrument use  
;
;
;
;
;-

PRO img_makebias, name, biasmean, root=root, side=side, path=path, instr=instr
 
common ccdinfo, xfull, yfull, nchip, namp, biaslev, satur, goodata
 
;;speed up if too many
if(n_elements(name) gt 20) then begin
    name=name[0:19]
    splog, "considering only 20 bias frames!"
endif

;;make storage for the mean bias
biasmean=make_array(xfull,yfull,/float)
nbia=n_elements(name)

;;loop over bias and make them
for pos=0, nbia-1 do begin
    ;;open files 
    splog, "Working on bias ", name[pos]
    ;;load chip (telescope specific)
    ;;layout is a structure that contain specific information on the
    ;;data section in each amplifier
    if(instr eq 'LRIS') then begin
        head=headfits(path+name[pos],exten=0)
        date=fxpar(head,"DATE")
        if(date gt '2009-04-01') then chip=lris_readfits(path+name[pos],/notrim,layout=layout)
        if(date lt '2009-04-01') then chip=lris_readold(path+name[pos],layout=layout)
    endif
    if(instr eq 'LRISr') then chip=lris_readold(path+name[pos],layout=layout)
    if(instr eq 'LBC') then  lbtc_readfits,path+name[pos],mosaic=chip,/notrim,layout=layout,/nodisp
    if(instr eq 'ESI') then chip=esi_readfits(path+name[pos],layout=layout)
    if(instr eq 'FEINC') then chip=feinc_readfits(path+name[pos],layout=layout)
    if(instr eq 'IMACS') then chip=imacs_readfits(path+name[pos],layout=layout)
    
    ;;go for oscan subtraction
    img_oscan, chip, biasimage, layout=layout, instr=instr
    undefine, chip
    ;;do median bias for CR
    if(pos eq 0) then begin
       size=size(biasimage)
       biasmean=fltarr(nbia,size[1],size[2]) 
       biasmean[0,*,*]=biasimage
    endif else biasmean[pos,*,*]=biasimage
 endfor

;;combine the final bias
splog, "Creating mean bias... "  
;;biasmean=biasmean/nbia
biasmean=djs_median(biasmean,1)

;;save median bias
mwrfits, biasmean, 'proc/'+side+root+"medbias.fits", /create
splog, "All done with meanbias ", side+root+"medbias.fits"

  
end
