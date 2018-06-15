;This procedure calls genmf to produce a grid 
;of halo mass functions.
;for different cosmological parameters. 
;Be aware that the input power spectrum is the one that comes 
;with the code (based on WMAP3), hence sigma_8 is tuned onto that.


;type     --> the type of halo you want
;
; 
;             opt = Option for mass function 
;             opt = 0 -Reed et al 2006, with n_eff dependence
;             opt = 1 -Reed et al 2006, without n_eff dependence
;             opt = 2 -Sheth-Tormen
;             opt = 3 -Jenkins et al 2001
;             opt = 4 -Warren et al 2005
;             opt = 5 -Reed et al 2003
;             opt = 6 -Press-Schechter

;nameout  --> the name of models
;zrange   --> the interval of redshift
;zbin     --> the step in redshift 
;_EXTRA   --> the cosmo parameters supported by m_cosm_commom



pro make_hm_grid, type, nameout, zrange=zrange, zbin=zbin, _EXTRA=extra, silent=silent


common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, cosm_L, cosm_r, sigma_8

;initialize common block
m_cosm_common, SILENT=silent, _EXTRA=extra


;set default
if ~keyword_set(zrange) then zrange=[0.,5.]
if ~keyword_set(zbin) then zbin=0.01


;make grid
z_grid=mkarr(zrange[0],zrange[1],zbin)



;get models
for z=0, n_elements(z_grid)-1 do begin

    splog, 'Compute redhsift ', z_grid[z]
    strz=strmid(rstring(z_grid[z]),0,6)
    name=nameout+strz
    
    spawn, ['./genmf',string(cosm_dm),rstring(cosm_L),rstring(sigma_8),rstring(z_grid[z]),name,rstring(type)], /noshell
    
   
endfor






end
