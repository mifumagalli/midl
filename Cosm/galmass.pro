;these abudnance matching fits determined/compiled by Kyle Stewart.
;Also available as a smartphone app on the Android Market as "GalMass"
;Documentation and further information available at 
;http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1109.3207

;print out help screen
pro galmasshelp
  print, "GalMass uses abundance matching to convert between central galaxy"
  print, "stellar mass (Mstar) and halo virial mass (Mvir)."
  print, " "
  print, "GalMass then uses fits derived from galaxy gas fraction observations"
  print, "at z=0 (McGaugh 2005)and z=2 estimates (Erb et al. 2006) to estimate"
  print, "Mgas as a function of Mstar or vice versa."
  print, " "
  print, "The 3 abubundance matching models (Mstar<->Mvir) are:"
  print, "  "
  print, "Conroy & Wechsler 2009:"
  print, "  called by x=MhaloToMstarCW(Mhalo,z)"
  print, "     and by y=MstarToMhaloCW(Mstar,z)"
  print, "  "
  print, "Moster et al. 2010:"
  print, "  called by x=MhaloToMstarM(Mhalo,z)"
  print, "     and by y=MstarToMhaloM(Mstar,z)"
  print, "  "
  print, "Behroozi et al. 2010:"
  print, "  called by x=MhaloToMstarB(Mhalo,z)"
  print, "     and by y=MstarToMhaloB(Mstar,z)"
  print, "  "
  print, "The galaxy gas model (Mstar<->Mgas) is given by Stewart et al. 2009: "
  print, "  called by x=MstarToMgas(Mstar,z,sig)"
  print, "     and by y=MgasToMstar(Mgas ,z,sig)"
  print, "  "
  print, "  where the sig variable governs the scatter in the relation:"
  print, "  sig=0 looks at only the average relation"
  print, "  sig=[some value] gives a random log-scatter with stdev=[some value]"
  print, "  sig=-1 gives a random z-dependent log-scatter according to Stewart et al. 2009"
  print, "  "
  print, "Note that the various relations are typically valid for z<2 or so (and most robust for z<1)."
end


;this function uses an analytic function derived from Conroy & Wechsler 09
;only valid for z<=2
function MhaloToMstarCW, Mhalo, z
  M1=Mhalo
  Zgas = z < 2 ;can't assign correctly beyond z=2 -- no stellar masses there.
  Mcs  = 10.0^(0.056*Zgas^2 + .068*Zgas + 9.5)
  Mch  = 10.0^(0.320*Zgas^2 + .018*Zgas +11.2)
  alpha = .021*Zgas^(4.86) + 3.39
  beta  = .085*Zgas + .36
  M2 = alog10(Mcs/(Mch^beta *2.0^(beta-alpha))) + (beta-alpha)*alog10(M1+Mch) + alpha*alog10(M1)
  M2=10.0^(M2)
  return, M2
end


;this function uses an analytic function derived from Conroy & Wechsler 09
;it is only valid for z<=2
function MstarToMhaloCW, Mstar, z
  M1=Mstar
  Zgas = z < 2 ;can't assign correctly beyond z=2 -- no stellar masses there.
  Mcs   = 10.0^(0.3517*Zgas^2 -0.2847*Zgas + 12.359)
  Mch   = 10.0^(0.1018*Zgas^2 -0.0885*Zgas +  10.7395)
  alpha = -0.0183*Zgas^2 + 0.01726*Zgas + 0.3238
  beta  = 0.1214*Zgas^2 -0.9053*Zgas + 3.170
  M2 = alog10(Mcs/(Mch^beta *2.0^(beta-alpha))) + (beta-alpha)*alog10(M1+Mch) + alpha*alog10(M1)
  M2=10.0^(M2)
  return, M2
end


;this function uses an analytic function derived from data points
;created by the fitting function to Moster et al 2010
;use same functional form ast the Conroy & Wechsler fits, for convenience
function MstarToMhaloM, Mstar, z
  M1=Mstar
  Zgas = z < 3 ;can't assign correctly beyond z=3 -- no stellar masses there.
  Mcs   = 10.0^(0.0222*Zgas^2 -0.0185*Zgas + 12.544)
  Mch   = 10.0^(0.03205*Zgas^2 -0.13775*Zgas + 10.846)
  alpha = 0.4733*Zgas - 0.03296
  beta  = -0.588*Zgas^0.48886 + 2.64
  ;print, Zgas, Mcs, Mch, alpha, beta
  M2 = alog10(Mcs/(Mch^beta *2.0^(beta-alpha))) + (beta-alpha)*alog10(M1+Mch) + alpha*alog10(M1)
  M2=10.0^(M2)
  return, M2
end

;this function uses an analytic function derived from abundance matching,
;equation given in from Moster, Somerville et al 2010
;claims to be valid to z<=3
function MhaloToMstarM, Mhalo, z
  Zgas = z < 3 ;can't assign correctly beyond z=3 -- no stellar masses there.
  ;four parameters, mstar/Mhalo, B, gamma, M1
  Mhalo_new=Mhalo ;/0.7 ;get into Msun units, without h factor
  mM= 0.0282 * (1+Zgas)^(-0.72)
  B= 1.06 + 0.17*Zgas
  gamma= 0.556 * (1+Zgas)^(-0.26)
  M1= 10.0^( 11.884 * (1+Zgas)^(0.019) )
  M2= 2 * Mhalo_new * mM * (  (Mhalo_new/M1)^(-1*B) + (Mhalo_new/M1)^(gamma) )^(-1.0)
  return, M2
end

;this function uses an analytic function from Behroozi et al 2010
function MstarToMhaloB, Mstar, z
  ms=Mstar
  zs=z < 2  ;can't assign correctly beyond z=2 -- no stellar masses there.
  M0 = 10.0^( 10.72   - (zs/(1+zs))* 0.59 )
  M1 = 10.0^( 12.35   - (zs/(1+zs))* 0.30 )
  beta     = ( 0.43   - (zs/(1+zs))* 0.18  )
  delta    = ( 0.56   - (zs/(1+zs))* 0.18  )
  gamma    = ( 1.54   - (zs/(1+zs))* 2.52  )
  mv =  alog10(M1) + beta*alog10(ms/M0) - 0.5 + ((ms/M0)^delta) / ((1+(ms/M0)^(-1.0*gamma)) ) 
  mv=10.0^(mv)
  return, mv
end

  
;this function uses an analytic function derived from data points
;created by the fitting function to Behroozi et al 2010
;use same functional form as the Conroy & Wechsler fits, for convenience
function MhaloToMstarB, Mhalo, z
  M1=Mhalo
  Zgas = z < 2 ;can't assign correctly beyond z=2
  Mcs   = 10.0^(0.03495*Zgas^2 -0.1919*Zgas + 10.1994)
  Mch   = 10.0^(0.0050869*Zgas^2 +0.0029872*Zgas +  11.82366)
  alpha = -0.207588*Zgas^2 + 0.751719*Zgas + 2.42311
  beta  = 0.120483*Zgas^2 -0.099373*Zgas + 0.20647
  M2 = alog10(Mcs/(Mch^beta *2.0^(beta-alpha))) + (beta-alpha)*alog10(M1+Mch) + alpha*alog10(M1)
  M2=10.0^(M2)
  return, M2
end


;function returns gas mass for given stellar mass by equation for Mgas/Mstar vs. Mstar
;as described in Stewart et al. 2009 (gas-rich mergers)
;only valid for z<=2
;sig -- impose a scatter in the average relation
;       sig=0  -> only return the average value
;       sig=-1 -> add random z-dependent scatter as outlined in Stewart et al.
function MstarToMgas, Mstellar, z, sig
  M1=Mstellar
  Zgas = z < 2 ;can't assign correctly beyond z=2 -- no stellar masses there.
  A=[.04, .587, .45]
  mu = A[1]*(1.0+Zgas)^A[2]
  Rz = A[0] * (M1/4.5e11)^( -1*mu )
  Rz=alog10( Rz )  ;work in log space, for now
  if sig lt 0 then sig=(0.34) - .19*alog10(1.0+Zgas)
  Rz = Rz + randomn(seed, size(Rz,/N_ELEMENTS))*sig  ;deviate from the "average"
  Mg=(10.0^Rz) * M1  
  return, Mg
end

;function returns stellar mass for given galaxy gas mass by equation for Mgas/Mstar vs. Mstar
;as described in Stewart et al. 2009 (gas-rich mergers)
;only valid for z<=2
;sig -- impose a scatter in the average relation
;       sig=0  -> only return the average value
;       sig=-1 -> add random z-dependent scatter as outlined in Stewart et al.
function MgasToMstar, Mgas, z, sig
  M1=Mgas
  Zgas = z < 2 ;can't assign correctly beyond z=2 -- no stellar masses there.
  A=[25, .587, .45]
  mu = A[1]*(1.0+Zgas)^A[2]
  Rz = A[0] * (M1/1.8e10)^( mu/(1-mu) )
  Rz=alog10( Rz )  ;work in log space, for now
  if sig lt 0 then sig=(0.34) - .19*alog10(1.0+Zgas)
  Rz = Rz + randomn(seed, size(Rz,/N_ELEMENTS))*sig  ;deviate from the "average"
  Ms=(10.0^Rz) * M1  
  return, Ms
end

