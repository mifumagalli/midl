;+
;
;
;  Take a ray output file from runhalo and make an 
;  estimate of the lofx in that box.
;
;
;-



pro find_lofx_rays, datafile, scale=scale


  ;;open data file 
  str=mrdfits(datafile,1)

  
  ;;turn scale factor into redshift 
  redshift=1./scale-1

  ;;approximate the redshift separation with finite delta
  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, $
          cosm_L, cosm_r, sigma_8
  m_cosm_common, /wmap7
  delta_z=0.001*str.boxsize*sqrt(cosm_dm*(1+redshift)^3+cosm_L)*(1+redshift)*(cosm_h/2.99d5)
  
  ;;find limits 
  limits=where(str.halo_raylimit gt 0,numlim)
  ;;total path 
  path=str.numrays*delta_z
  
  dndz=1.*numlim/path
  
  ;;get dzdx
  hz=m_cosm_hubble(redshift,/wmap7)
  dxdz=cosm_h/hz*(1+redshift)^2
  dndx=dndz/dxdz
  
  splog, "Redshift, Path, dxdz, dndz, dndx"
  splog, redshift, path, dxdz, dndz, dndx
  
end


