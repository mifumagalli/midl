pro sim_physunits, phys

phys        = create_struct('Grav',6.673d-8, $           ; Newton's gravity constant [cgs]
                            'c',2.99792458d10,    $        ; speed of lights [cm/s]
                            'm_H',1.67262158d-24, $        ; proton mass; [g]
                            'msun', 1.989d33, $            ; solar mass [g]
                            'lsun', 3.862d33, $            ; solar luminosity [erg/s]
                            'yr',3.1558d7, $               ; yr [s]
                            'Mpc',3.0856776d+24, $         ; Mpc in cm
                            'sigma0',6.3462964d-18, $      ; HI photo-ionisation cross section [cm^2]
                            'sigma_Lya',4.4813564d-18, $   ; Lya cross section in [cm^2]
                            'f',0.4164, $                  ; f-value
                            'lya',1215.67d-8, $            ; cm
                            'sigmaT',6.65250d-25, $        ; Thomson cross section [cm^2]
                            'e',4.80320425d-10, $          ; electron charge [esu]
                            'me',9.10938291d-28, $         ; electron mass [g]
                            'pc',3.0856776d+18, $          ; pc in cm
                            'kpc',3.0856776d+21,$          ; kpc in cm
                            'km',1d5,$                    ; km in cm
                            'ccms',2.99792458d10)         ; cm/s

end
