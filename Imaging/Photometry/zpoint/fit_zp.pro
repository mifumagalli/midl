;+
;program that reads a table of mag info
;and fit a zp
;
;table_file has to be formatted in the following way
;star_name, mag_instr, err_magintr, mag_true, color_instr, air_mass
;
;
;nclr       --> If set, do not fit for color term
;am_coeff   --> If set to a number, the air_mass term is not fit.
;fix_color  --> If set to a number, the color term will be subtracted
;           --> before fitting the zero point. (Supported by fixed AM only)
;ab_mag     --> if set to a filter (U,B,V,R,I), return ZP in AB mag.          
;apply      --> If set to image, it writes the solution in the header
;
;-

pro fit_zp, table_file, nclr=nclr, am_coeff=am_coeff, ab_mag=ab_mag, apply=apply, $
            fix_color=fix_color


readcol, table_file, star_name, mag_instr, err_magintr, mag_true, $
  color_instr, air_mass, format='a,f,f,f,f,f'


;if set AB, define shift

if keyword_set(ab_mag) then begin

    type_system=' AB mag'
    ;offsets are from  Cooke at al. 2005 
    case ab_mag of 
        'u': zp_ab=0.70
        'b': zp_ab=-0.15
        'v': zp_ab=-0.01
        'r': zp_ab=0.18
        'i': zp_ab=0.43
        
        else: print, 'No offset for this filter!'
    endcase
endif else begin
    type_system=' Johnson-Kron-Cousins mag'
    zp_ab=0.
endelse

splog, 'System: ', type_system

;open image and find stuff
if keyword_set(apply) then begin
    fits=mrdfits(apply,0,header)
endif






;calibrate if no color, fix AM
if keyword_set(nclr) and keyword_set(am_coeff) then begin

    ;if fixed color
    if keyword_set(fix_color) then corr=-fix_color*color_instr else corr=0.
    
    
    ;get ZP
    x_photsol1, mag_instr+corr, mag_true, err_magintr, air_mass, am_coeff, $
      coeffs, sigma_coeff, CHISQ=chi 
    
    splog, 'ZP: ', coeffs+zp_ab
    splog, 'eZP: ', sigma_coeff
    splog, 'AM: ', abs(am_coeff)
    splog, 'eAM: ', 0.
    if keyword_set(fix_color) then begin
        splog, 'CLR: ', fix_color
        splog, 'eCLR: ', 0.
    endif
    splog, 'Chi^2: ', chi 

    ;plot each star
    m_calib=mag_instr+coeffs-abs(air_mass*am_coeff)
    
    plot, mag_true, m_calib-mag_true, psym=6, $
      xtitle='Star Mag', ytitle='Residual', /ynozero
    
    if(n_elements(m_calib) gt 1) then splog, 'Mean ABS(Residual) ', $
      mean(abs(m_calib-mag_true))

    ;apply 
    if keyword_set(apply) then begin
        fxaddpar, header, 'P_COMM1', 'Photometry done with fit_zp on '+systime()
        fxaddpar, header, 'P_COMM2', 'm_true=m_inst+P_ZP-P_AM*am'
        fxaddpar, header, 'P_ZP', coeffs+zp_ab, 'Photometry on '+type_system
        fxaddpar, header, 'P_ERR_ZP', sigma_coeff
        fxaddpar, header, 'P_AM', abs(am_coeff), 'Fixed AM'
        if keyword_set(fix_color) then begin
           ;prompt for color info
            head_color=''
            read, 'Specify color: ', head_color
            fxaddpar, header, 'P_CLR', fix_color, 'Fixed term. Color '+head_color
        endif
        fxaddpar, header, 'P_RESID',  mean(abs(m_calib-mag_true))
    endif
    
    
endif




;calibrate if only fix AM
if ~keyword_set(nclr) and keyword_set(am_coeff) then begin


    ;get ZP
    x_photsol2, mag_instr, mag_true, err_magintr, air_mass, $
      coeffs, sigma_coeff, setam=am_coeff, color=color_instr,  CHISQ=chi 

    splog, 'ZP: ', coeffs[0]+zp_ab
    splog, 'eZP: ', sigma_coeff[0]
    splog, 'AM: ', am_coeff
    splog, 'eAM: ', 0.
    splog, 'CLR: ', coeffs[1]
    splog, 'eCLR: ', sigma_coeff[1]
    splog, 'Chi^2: ', chi 

    ;plot each star (nb x_photsol2 already includes Am term in mag_instr)
    m_calib=mag_instr+coeffs[0]-coeffs[1]*color_instr
    
    plot, mag_true, m_calib-mag_true, psym=6, $
      xtitle='Star Mag', ytitle='Residual', /ynozero
    
    if(n_elements(m_calib) gt 1) then splog, 'Mean ABS(Residual) ', $
      mean(abs(m_calib-mag_true))


    ;apply 
    if keyword_set(apply) then begin
    
        ;prompt for color info
        head_color=''
        read, 'Specify color: ', head_color
        fxaddpar, header, 'P_COMM1', 'Photometry done with fit_zp on '+systime()
        fxaddpar, header, 'P_COMM2', 'm_true=m_inst+P_ZP-P_AM*am-P_CLR*clr'
        fxaddpar, header, 'P_ZP', coeffs[0]+zp_ab, 'Photometry on '+type_system
        fxaddpar, header, 'P_ERR_ZP', sigma_coeff[0]
        fxaddpar, header, 'P_AM', am_coeff, 'Fixed AM'
        fxaddpar, header, 'P_CLR', coeffs[1], 'color '+head_color
        fxaddpar, header, 'P_ERR_CLR', sigma_coeff[1]
        fxaddpar, header, 'P_RESID',  mean(abs(m_calib-mag_true))
    endif
    

endif

 
;calibrate with everything
if ~keyword_set(nclr) and ~keyword_set(am_coeff) then begin
    

  ;get ZP
    x_photsol3, mag_instr, mag_true, err_magintr, air_mass, color_instr, $
      coeffs, sigma_coeff,  CHISQ=chi 
    
    splog, 'ZP: ', coeffs[0]+zp_ab
    splog, 'eZP: ', sigma_coeff[0]
    splog, 'AM: ', coeffs[1]
    splog, 'eAM: ', sigma_coeff[1]
    splog, 'CLR: ', coeffs[2]
    splog, 'eCLR: ', sigma_coeff[2]
    splog, 'Chi^2: ', chi 

    ;plot each star
    m_calib=mag_instr+coeffs[0]-coeffs[1]*air_mass-coeffs[2]*color_instr
    
    plot, mag_true, m_calib-mag_true, psym=6, $
      xtitle='Star Mag', ytitle='Residual', /ynozero
    
    if(n_elements(m_calib) gt 1) then splog, 'Mean ABS(Residual) ', $
      mean(abs(m_calib-mag_true))
                                
     ;apply
    if keyword_set(apply) then begin
        head_color=''
        read, 'Specify color: ', head_color
        
        fxaddpar, header, 'P_COMM1', 'Photometry done with fit_zp on '+systime()
        fxaddpar, header, 'P_COMM2', 'm_true=m_inst+P_ZP-P_AM*am-P_CLR*clr'
        fxaddpar, header, 'P_ZP', coeffs[0]+zp_ab, 'Photometry on '+type_system
        fxaddpar, header, 'P_ERR_ZP', sigma_coeff[0]
        fxaddpar, header, 'P_AM', coeffs[1]
        fxaddpar, header, 'P_ERR_AM', sigma_coeff[1]
        fxaddpar, header, 'P_CLR', coeffs[2], 'color '+head_color
        fxaddpar, header, 'P_ERR_CLR', sigma_coeff[2]
        fxaddpar, header, 'P_RESID',  mean(abs(m_calib-mag_true))
    endif

endif


splog, 'Calibrate with m_true=m_inst+ZP-AM*am[-CLR*clr]'

;write header
if keyword_set(apply) then mwrfits, fits, apply, header, /create


end

