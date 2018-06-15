
function bspline_radial_valu, r, theta,  sset, mode=mode

; WRITTEN: S. Burles, MIT
     
      nr = n_elements(r)
      nt = n_elements(theta)
      if (nr NE nt) OR (nt EQ 0) OR (nr EQ 0) then return, 0

      rs = sort(r)
      bf1 = bspline_action(r[rs], sset, lower=loweraction, upper=upperaction)
      yfit = r * 0.

      nord = sset.nord
      npoly = n_elements(sset.ntheta)
      filler = replicate(1,sset.nord)

      tempset = create_bsplineset(sset.fullbkpt, sset.nord)

      if n_elements(mode) EQ 1 then begin
        if (mode GE 0 AND mode LT npoly) then begin
         
          if sset.ntheta[mode] LT 0 then $
                   ft =sin(-1.0*sset.ntheta[mode]*theta[rs]) $
          else ft = cos(sset.ntheta[mode]*theta[rs])
           
          action = bf1 * (ft # filler)
          tempset.coeff = sset.coeff[mode,*]
          yfit[rs] = bspline_valu(r[rs], tempset, action=action, $
                                  upper=upperaction, lower=loweraction)
        endif
      endif else begin
        for j=0, npoly-1 do begin
          if sset.ntheta[j] EQ 0 then action = bf1 $
          else begin
            if sset.ntheta[j] LT 0 then ft =sin(-1.0*sset.ntheta[j]*theta[rs]) $
            else ft = cos(sset.ntheta[j]*theta[rs])
            action = bf1 * (ft # filler)
          endelse

          tempset.coeff = sset.coeff[j,*]
          yfit[rs] = yfit[rs] + bspline_valu(r[rs], tempset, action=action, $
                                  upper=upperaction, lower=loweraction)
           
        endfor
      endelse

   return, yfit
end

