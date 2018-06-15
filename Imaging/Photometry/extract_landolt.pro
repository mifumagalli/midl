;procedure that extract from the Landolt table
;a group of stars and add info from sdss
;Use list updated to 2009.

;group --> the lable name to identify the field
;write  --> if set, write a structure


pro extract_landolt, group, write=write


readcol, 'Landolt2009.txt', skipline=32, SGroup, Star, RAh, RAm, RAs, DEd, $
  DEm, DEs, Vmag, BV, UB, VR, RI, VI, n, m, e_Vmag, e_BV, e_UB, e_VR, e_RI, e_VI, $
  format='a,a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f'



;identify
index=where(strtrim(group,2) eq strtrim(SGroup,2), match)


if(match lt 1) then begin
    splog, 'I cannot find the group...'
    return
endif else begin

    ;extract
    Star    = 	Star[index]
    RAh     = 	RAh[index]
    RAm     = 	RAm[index]
    RAs     = 	RAs[index]
    DEd     = 	DEd[index]
    DEm     = 	DEm[index]
    DEs     = 	DEs[index]
    Vmag    = 	Vmag[index]
    BV	    = 	BV[index]
    UB	    = 	UB[index]
    VR	    = 	VR[index]
    RI	    = 	RI[index]
    VI	    = 	VI[index]
    n	    = 	n[index]
    m	    = 	m[index]
    e_Vmag  = 	e_Vmag[index]
    e_BV    = 	e_BV[index]
    e_UB    = 	e_UB[index]
    e_VR    = 	e_VR[index]
    e_RI    = 	e_RI[index]
    e_VI    = 	e_VI[index]

    ;get sdss
    
    sdss_name=strarr(match)
    u_sdss=fltarr(match)
    g_sdss=fltarr(match)
    r_sdss=fltarr(match)
    i_sdss=fltarr(match)
    z_sdss=fltarr(match)
    warning=fltarr(match)

        
    for i=0, match-1 do begin
        
        sta_sdss=queryvizier('II/294',[ten(RAh[i],RAm[i],RAs[i])*15.,ten(DEd[i],DEm[i],DEs[i])],0.03)
        
            
          ;check if found 
        fnd=size(sta_sdss,/type)
        if(fnd eq 8 ) then begin                
             ;if found store 
            if(n_elements(sta_sdss) gt 1) then begin
                splog, "WARNING: Found multiple SDSS stars for ", star[i], $
                " Picked ", sta_sdss[0].SDSS
         ;grab value (and set warning)
                warning[i]=1.
                sdss_name[i]=sta_sdss[0].SDSS
                u_sdss[i]=sta_sdss[0].umag
                g_sdss[i]=sta_sdss[0].gmag
                r_sdss[i]=sta_sdss[0].rmag
                i_sdss[i]=sta_sdss[0].imag
                z_sdss[i]=sta_sdss[0].zmag
                
            endif else begin
                
                sdss_name[i]=sta_sdss.SDSS
                u_sdss[i]=sta_sdss.umag
                g_sdss[i]=sta_sdss.gmag
                r_sdss[i]=sta_sdss.rmag
                i_sdss[i]=sta_sdss.imag
                z_sdss[i]=sta_sdss.zmag
                
            endelse
            
        endif else begin
            splog, "WARNING: I didn't find any SDSS star for ",  star[i]
        endelse
    endfor
    

    if keyword_set(write) then begin
    ;write a structure
        str_write={Star:Star, RAh:RAh, RAm:RAm,  RAs:RAs, DEd:DEd, DEm:DEm, DEs:DEs, Vmag:Vmag, $
                   BV:BV,  UB:UB,  VR:VR, RI:RI, VI:VI, n:n, m:m, e_Vmag:e_Vmag, e_BV:e_BV, $
                   e_UB:e_UB,  e_VR:e_VR, e_RI:e_RI, e_VI:e_VI, sdss_name:sdss_name, warning:warning,$
                   u_sdss:u_sdss, g_sdss:g_sdss, r_sdss:r_sdss, i_sdss:i_sdss, z_sdss:z_sdss}
        
        mwrfits, str_write, group+'_landolt.fits', /create
    endif

    ;print to the screan some basic stuff
    Bmag=BV+Vmag
    Umag=UB+Bmag
    Rmag=Vmag-VR
    Imag=Vmag-VI
    
    print, 'Mag in Johnson-Kron-Cousins system'
    print, 'Star ', 'Umag ', 'Bmag ', 'Vmag ', 'Rmag ', 'Imag '
    print, '----------------------------------------------------'
    forprint, star, Umag, Bmag, Vmag, Rmag, Imag, textout=1
    

    print, 'SDSS Mag (AB system)'
    print, 'Star', 'u ', 'g ', 'r ', 'i ', 'z '
    print, '------------------------------------'
    forprint, star, u_sdss, g_sdss, r_sdss, i_sdss, z_sdss, warning, textout=1

    

endelse




end
