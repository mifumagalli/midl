;+
;
; NAME: lac3d
;
; VERSION: 1.4, 25 Oct 11, Ric Davies
;
; PURPOSE: 
;  detect and correct bad pixels, which are identified using a 3D
;  Laplacian Edge Detection method based on van Dokkum's LACosmic.
;  The scheme is as follows:
;   1. remove continuum from spectrally extended sources (e.g. stars, AGN)
;   2. remove spatially extended sources (e.g. OH lines)
;   3. remove remaining fine structure (similar to LACosmic)
;   4. derive noise cube 
;   5. calculate bad pixel signal (i.e. Laplacian) in 3D
;   6. find bad pixels based on ratio of bad pixel signal to noise &
;      flux
;   7. repeat for negative bad pixels
;   8. interpolate bad pixels & add back in removed structure
;   9. update noise cube and repeat
;
; NOTES: 
;  this is designed to work with datacubes for which the spectral
;  axis is assumed to be 3 (but it doesn't really matter much).
;
;  On my machine it takes ~2mins per million pixels processed, and
;  most of this effort is in deriving the noise cube
;
;  While adjusting s_lim, set /wnoise the first time, and then
;  afterwards use noise='name_n' to read it back in. This saves a lot
;  of time.
;
; CALLING SEQUENCE:
;  lac3d, objcube, $
;    noise=noise, s_lim=s_lim, f_lim=f_lim, max_iter=max_iter, $
;    /nocont, /nospatial, /noneg, /alt, /wnoise, /wbadpix, $
;    /test, /nofine, msize_z=msize_z, msize_xy=msize_xy, $
;    msize_f=msize_f, msize_n=msize_n
;
; INPUT:
;
;  objcube = name of fits cube in which bad pixels should be detected
;            and corrected
;
; OPTIONAL INPUT & KEYWORDS:
;
; noise     name of fits noise cube associated with objcube. Generating
;           the noise cube is what takes the longest, so if you have
;           one, use it.
;
; s_lim     detection threshold as a ratio of the signal in the
;           Laplacian image to the estimated noise (with factor 2
;           scaling). Default value is 3. This parameter is likely to
;           need adjusting.
;
; f_lim     detection threshold as a ratio of the signal in the
;           Laplacian image to the signal in the structure image
;           (which contains point-source like things rather than
;           emission which is spectrally or spatially extended).
;           Default value = 2
;
; max_iter  maximum number of iterations for detecting bad pixels.
;           Default value = 5 (unless /test is set, then max_iter=1)
;
; /nocont   set keyword if continuum sources should not be traced and
;           removed before detecting bad pixels 
;
; /nospatial set if sources which are spatially (but not
;           necessarily spectrally) very extended should not be
;           traced and removed before detecting bad pixels 
;
; /noneg    do not look for or correct any negative bad pixels
;
; /nofine   do not subtract fine structure from cube (although it is
;           still calculated)
;
; /alt      use alternative convolution kernal which uses all 26
;           neighbours (including diagonals) rather than only 6
;           direct neighbours. Tests indicate this should be better.
;           Note it is also the Laplacian, just with a different sigma
;           for the Gaussian form which it is derived.
;
; /wnoise   force write out of derived noise cube (name is same as
;           input object cube, with _n appended), overwriting if
;           filename already exists.
;
; /wbadpix  write out cube indicating positions of detected bad
;           pixels (name is same as input cube, with _b
;           appended). This cube contains '1' at position of each
;           bad pixel and '0' everywhere else.
;
; more advance parameters:
;
; /test     write out extra files useful for debugging
;
; msize_z  is median filter size for continuum source estimation
;          default is 49, which uses a [1,1,49] filter.
;
; msize_xy is median filter size for spatial source estimation
;          default is 7, which uses a [7,7,1] filter
;
; msize_f  is median filter size for fine structure estimation
;          default is 3, which uses a [3,3,3] filter
;
; msize_n  is median filter size for noise estimation
;          default is 5, which uses a [5,5,5] filter
;
;----------------------------------------------------
;
;updates since version 0.5
;19.02.07 - added option of using modified convolution kernal
;         \ and keywords for writing out cubes for noise & bad pixels
;26.02.07 - used boundary extend for spatial structure filter
;         \ revised intro material
;         \ swapped order of functions/procedures 
;08.05.07 - corrected error in noise scaling for detection threshold
;14.11.07 - split cube into blocks (to keep memory usage manageable)
;16.11.07 - iteratively update noise & fine structure
;         \ corrected a small error in output cube (oops!)
;19.11.07 - modified usage of /wnoise keyword; and also s_lim_default
;21.11.07 - corrected error in noise estimation (oops!)
;         \ added option to not correct negative bad pixels
;22.11.07 - modified slightly pclip estimation for noise
;04.12.07 - modified pclip thresholds
;15.07.08 - added keywords for fine structure & median filter sizes
;25.10.11 - ignore edge planes filled with at least half NaN values
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function boundary_extend, d_raw
;wrap 3-pixel boundaries to extend each axis by 3 in each direction

dsize = size(d_raw)
xsize = dsize[1]
ysize = dsize[2]
zsize = dsize[3]

d = fltarr(xsize+6,ysize+6,zsize+6)
d[3:xsize+2,3:ysize+2,3:zsize+2] = d_raw
;add front end
d[3:xsize+2,3:ysize+2,0:2] = d_raw[*,*,zsize-3:zsize-1]
;add back end
d[3:xsize+2,3:ysize+2,zsize+3:zsize+5] = d_raw[*,*,0:2]
;add left side
d[0:2,3:ysize+2,*] = d[xsize+0:xsize+2,3:ysize+2,*]
;add right side
d[xsize+3:xsize+5,3:ysize+2,*] = d[3:5,3:ysize+2,*]
;add top
d[*,ysize+3:ysize+5,*] = d[*,3:5,*]
;add bottom
d[*,0:2,*] = d[*,ysize+0:ysize+2,*]

return,d
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lac3dn, cube, noise=noise, s_lim=s_lim, f_lim = f_lim, max_iter=max_iter, nocont=nocont, nospatial=nospatial, noneg=noneg, alt=alt, wnoise=wnoise, wbadpix=wbadpix, test=test, msize_z=msize_z, msize_xy=msize_xy,msize_f=msize_f,msize_n=msize_n,nofine=nofine

resolve_all,/continue_on_error,/quiet

t_start = systime(1)

;----------------------------
;default parameters

s_lim_default = 3
f_lim_default = 2

max_iter_default = 5

msize_z_default  = 49       ; box size for continuum source estimation
msize_xy_default =  7       ; box size for spatial source estimation
msize_f_default  =  3       ; box size for fine structure estimation
msize_n_default  =  5       ; box size for noise estimation


;rel_scale = 0.0

;----------------------------

;define laplacian

lap = [[[0.,0,0],[0,-1,0],[0,0,0]],$
       [[0,-1,0],[-1,6,-1],[0,-1,0]],$
       [[0,0,0],[0,-1,0],[0,0,0]]]

if keyword_set(alt) then $
lap = [[[-1.,-1,-1],[-1,-1,-1],[-1,-1,-1]],$
       [[-1,-1,-1],[-1,26,-1],[-1,-1,-1]],$
       [[-1,-1,-1],[-1,-1,-1],[-1,-1,-1]]]

lap = 1./max(lap) * lap

;check parameters
if (strmid(cube,4,5,/reverse_offset) ne '.fits') then cube=cube+'.fits'

;mk_noise = 1
;if (keyword_set(noise)) then mk_noise = 0

if (not keyword_set(s_lim)) then s_lim = s_lim_default
if (not keyword_set(f_lim)) then f_lim = f_lim_default
if (not keyword_set(max_iter)) then max_iter = max_iter_default
if (not keyword_set(msize_z)) then msize_z = msize_z_default
if (not keyword_set(msize_xy)) then msize_xy = msize_xy_default
if (not keyword_set(msize_f)) then msize_f = msize_f_default
if (not keyword_set(msize_n)) then msize_n = msize_n_default
;if (keyword_set(test)) then max_iter = 1

;read object data
d_raw0 = readfits(cube,hdr)
dsize = size(d_raw0)
xsize0 = dsize[1]
ysize0 = dsize[2]
zsize0 = dsize[3]

;trim of excess NaNs (but put back later)
;left edge
fp = 0
n1 = 0
while (fp lt 1) do begin
    slice = d_raw0[n1,*,*]
    npix = n_elements(slice)
    junk = where(finite(slice),nfin)
    if (nfin ge npix/2.) then fp = 1 else n1 = n1 + 1
    ;print,fp,n1,nfin,npix
endwhile
;right edge
fp = 0
n2 = xsize0-1
while (fp lt 1) do begin
    slice = d_raw0[n2,*,*]
    npix = n_elements(slice)
    junk = where(finite(slice),nfin)
    if (nfin ge npix/2.) then fp = 1 else n2 = n2 - 1
    ;print,fp,n2,nfin,npix
endwhile
;bottom edge
fp = 0
n3 = 0
while (fp lt 1) do begin
    slice = d_raw0[*,n3,*]
    npix = n_elements(slice)
    junk = where(finite(slice),nfin)
    if (nfin ge npix/2.) then fp = 1 else n3 = n3 + 1
    ;print,fp,n3,nfin,npix
endwhile
;top edge
fp = 0
n4 = ysize0-1
while (fp lt 1) do begin
    slice = d_raw0[*,n4,*]
    npix = n_elements(slice)
    junk = where(finite(slice),nfin)
    if (nfin ge npix/2.) then fp = 1 else n4 = n4 - 1
    ;print,fp,n2,nfin,npix
endwhile
;front edge
fp = 0
n5 = 0
while (fp lt 1) do begin
    slice = d_raw0[*,*,n5]
    npix = n_elements(slice)
    junk = where(finite(slice),nfin)
    if (nfin ge npix/2.) then fp = 1 else n5 = n5 + 1
    ;print,fp,n5,nfin,npix
endwhile
;back edge
fp = 0
n6 = zsize0-1
while (fp lt 1) do begin
    slice = d_raw0[*,*,n6]
    npix = n_elements(slice)
    junk = where(finite(slice),nfin)
    if (nfin ge npix/2.) then fp = 1 else n6 = n6 - 1
    ;print,fp,n6,nfin,npix
endwhile
print,format='(a,6i)','ignoring edges to:',n1,n2,n3,n4,n5,n6
d_raw =  d_raw0[n1:n2,n3:n4,n5:n6]
dsize = size(d_raw)
xsize = dsize[1]
ysize = dsize[2]
zsize = dsize[3]

;stop

;replace nan with interpolated values or zeros
nan_r = where(finite(d_raw,/nan),nan_i)
if (nan_i gt 0) then begin
    print,'Replacing '+strtrim(string(nan_i),2)+' non-finite pixels'
    d_ext = boundary_extend(d_raw)
    repval = fltarr(nan_i)
    pos = array_indices(d_raw,nan_r)
    xpos = 3+reform(pos[0,*])
    ypos = 3+reform(pos[1,*])
    zpos = 3+reform(pos[2,*])
    xmin = xpos-1
    xmax = xpos+1
    ymin = ypos-1
    ymax = ypos+1
    zmin = zpos-1
    zmax = zpos+1
    for i=long(0),nan_i-1 do begin
        sec = d_ext[xmin[i]:xmax[i],ymin[i]:ymax[i],zmin[i]:zmax[i]]
        secfr = where(finite(sec),secfi)
        if (secfi gt 0) then repval[i] = median(sec[secfr])
    endfor
    d_raw[nan_r] = repval
endif

;there should now be only finite pixels

;unless /nocont is set then remove continuum sources & subtract from raw cube
dz49 = fltarr(xsize,ysize,zsize)
if not(keyword_set(nocont)) then begin
    print,'Tracing continuum sources'
    halfm = (msize_z-1)/2
    for iz=0,zsize-1 do begin
        zmin = max([0,iz-halfm])
        zmax = min([zsize-1,iz+halfm])
        for iy=0,ysize-1 do begin
            for ix=0,xsize-1 do begin
            dz49[ix,iy,iz] = median(d_raw[ix,iy,zmin:zmax])
        endfor
    endfor
endfor  
d_raw = temporary(d_raw) - dz49
endif

;unless /nospatial is set then remove sources which are spatial only,
; & subtract from raw cube
dxy7 = fltarr(xsize,ysize,zsize)
if not(keyword_set(nospatial)) then begin
    print,'Tracing spatial sources'
    halfm = (msize_xy-1)/2
    d_ext = boundary_extend(d_raw)
    for iy=0,ysize-1 do begin
        ymin = iy+3-halfm
        ymax = iy+3+halfm
        for ix=0,xsize-1 do begin
            xmin = ix+3-halfm
            xmax = ix+3+halfm
            for iz=0,zsize-1 do begin
                dxy7[ix,iy,iz] = median(d_ext[xmin:xmax,ymin:ymax,iz+3])
        endfor
    endfor
endfor
d_raw = temporary(d_raw) - dxy7
endif

d_out = d_raw

;calculate true fine structure;
; subtract from raw cube
dm3 = fltarr(xsize,ysize,zsize)
print,'Generating fine structure cube'
halfm = (msize_f-1)/2
d_ext = boundary_extend(d_out)
for iz=0,zsize-1 do begin
    zmin = iz+3-halfm
    zmax = iz+3+halfm
    for iy=0,ysize-1 do begin
        ymin = iy+3-halfm
        ymax = iy+3+halfm
        for ix=0,xsize-1 do begin
            xmin = ix+3-halfm
            xmax = ix+3+halfm
            dm3[ix,iy,iz] = median(d_ext[xmin:xmax,ymin:ymax,zmin:zmax])
        endfor
    endfor
endfor  
f = abs(dm3)
;f = abs(dm3 + rel_scale*(dxy7 + dz49)) ; this fraction is to take account 
                                ; of the fact that the signal in the 
                                ; highly smoothed cubes is more certain

if not(keyword_set(nofine)) then d = d_out - dm3 else d = d_out

;create noise cube if it is not supplied
if (keyword_set(noise)) then begin
    if (strmid(noise,4,5,/reverse_offset) ne '.fits') then noise=noise+'.fits'
    nc = readfits(noise)
endif else begin

; this is done by looking at stddev of each 5x5x5 pixels after rejecting
; deviant ones using pclip
    print,'Estimating noise map'
    halfm = (msize_n-1)/2
    nc = fltarr(xsize,ysize,zsize)
    d_ext = boundary_extend(d)
    for iz=0,zsize-1 do begin
        zmin = iz+3-halfm
        zmax = iz+3+halfm
        for iy=0,ysize-1 do begin
            ymin = iy+3-halfm
            ymax = iy+3+halfm
            for ix=0,xsize-1 do begin
                xmin = ix+3-halfm
                xmax = ix+3+halfm
                sec = d_ext[xmin:xmax,ymin:ymax,zmin:zmax]
                secfr = where(finite(sec) and sec ne 0.00,secf_i)
                if (secf_i gt 52) then begin
;;;;; Here's the rejection method
                    secf = sec[secfr]
                    medval = median(secf,/even)
                    abssecf = abs(secf-medval)
                    med_o = abssecf[sort(abssecf)]
                    n_med = n_elements(secf)-1
                    pclip = med_o[fix(0.8*n_med)]
                    ppix = where(abssecf le pclip*5)
                    pvec = secf[ppix]
                    pmed = median(pvec,/even)
                    psdv = stddev(pvec)
                    goodpix = where(abs(pvec-pmed) le psdv*3)
                    nc[ix,iy,iz] = stddev(pvec[goodpix])
;;;; that was it
                endif
            endfor
        endfor
    endfor
endelse

;finished creating noise map: nc

;now start the iterations on the Laplacian & bad pixel detection
r_bp_tot = 0
loop = 1
new_bad = 11
badpix = fltarr(xsize,ysize,zsize)
badpixn = fltarr(xsize,ysize,zsize)
lp = fltarr(xsize,ysize,zsize)
ln = fltarr(xsize,ysize,zsize)
while (loop le max_iter and new_bad gt 10) do begin
print,'iteration '+strtrim(string(loop),2)

if (loop gt 1) then begin

;re-estimate fine structure (around bad pixels)
halfm = (msize_f-1)/2
r_bp_cube = fltarr(xsize,ysize,zsize)*0.
r_tmp = [r_bpp,r_bpn]
r_bp_pn = r_tmp[where(r_tmp ge 0)]
r_bp = array_indices(r_bp_cube,r_bp_pn)
for i=long(0),n_elements(r_bp_pn)-1 do begin
    xmin = max([r_bp[0,i]-halfm,0])
    xmax = min([r_bp[0,i]+halfm,xsize-1])
    ymin = max([r_bp[1,i]-halfm,0])
    ymax = min([r_bp[1,i]+halfm,ysize-1])
    zmin = max([r_bp[2,i]-halfm,0])
    zmax = min([r_bp[2,i]+halfm,zsize-1])
    r_bp_cube[xmin:xmax,ymin:ymax,zmin:zmax] = 1
endfor
r_bp = array_indices(r_bp_cube,where(r_bp_cube eq 1,i_bp))
print,'updating '+strtrim(string(n_elements(r_bp)/3),2)+$
  ' pixels in fine structure cube'
d_ext = boundary_extend(d_out)
for ibp = long(0),i_bp-1 do begin
    ix = reform(r_bp[0,ibp])
    iy = reform(r_bp[1,ibp])
    iz = reform(r_bp[2,ibp])
    zmin = iz+3-halfm
    zmax = iz+3+halfm
    ymin = iy+3-halfm
    ymax = iy+3+halfm
    xmin = ix+3-halfm
    xmax = ix+3+halfm
    dm3[ix,iy,iz] = median(d_ext[xmin:xmax,ymin:ymax,zmin:zmax])
endfor
f = abs(dm3)
;f = abs(dm3 + rel_scale*(dxy7 + dz49))
if not(keyword_set(nofine)) then d = d_out - dm3

;re-estimate noise (around bad pixels)
if (not (keyword_set(noise))) then begin
halfm = (msize_n-1)/2
r_bp_cube = fltarr(xsize,ysize,zsize)*0.
r_tmp = [r_bpp,r_bpn]
r_bp_pn = r_tmp[where(r_tmp ge 0)]
r_bp = array_indices(r_bp_cube,r_bp_pn)
for i=long(0),n_elements(r_bp_pn)-1 do begin
    xmin = max([r_bp[0,i]-halfm,0])
    xmax = min([r_bp[0,i]+halfm,xsize-1])
    ymin = max([r_bp[1,i]-halfm,0])
    ymax = min([r_bp[1,i]+halfm,ysize-1])
    zmin = max([r_bp[2,i]-halfm,0])
    zmax = min([r_bp[2,i]+halfm,zsize-1])
    r_bp_cube[xmin:xmax,ymin:ymax,zmin:zmax] = 1
endfor
r_bp = array_indices(r_bp_cube,where(r_bp_cube eq 1,i_bp))
print,'updating '+strtrim(string(n_elements(r_bp)/3),2)+$
  ' pixels in noise cube'
d_ext = boundary_extend(d)
for ibp = long(0),i_bp-1 do begin
    ix = reform(r_bp[0,ibp])
    iy = reform(r_bp[1,ibp])
    iz = reform(r_bp[2,ibp])
    zmin = iz+3-halfm
    zmax = iz+3+halfm
    ymin = iy+3-halfm
    ymax = iy+3+halfm
    xmin = ix+3-halfm
    xmax = ix+3+halfm
    sec = d_ext[xmin:xmax,ymin:ymax,zmin:zmax]
    secfr = where(finite(sec) and sec ne 0.00,secf_i)
    if (secf_i gt 52) then begin
                    secf = sec[secfr]
                    medval = median(secf,/even)
                    abssecf = abs(secf-medval)
                    med_o = abssecf[sort(abssecf)]
                    n_med = n_elements(secf)-1
                    pclip = med_o[fix(0.8*n_med)]
                    ppix = where(abssecf le pclip*10)
        pvec = secf[ppix]
        pmed = median(pvec,/even)
        psdv = stddev(pvec)
        goodpix = where(abs(pvec-pmed) le psdv*5)
        nc[ix,iy,iz] = stddev(pvec[goodpix])
    endif
endfor
endif

endif ; loop gt 1

;start on the Laplacian
print,'Creating Laplacian'

;do the following in sections
zmin=-101
zmax=1
while(zmax lt zsize) do begin
 zmin = zmin+100
 zmax = zmax+100
 tzmin = max([0,zmin])
 tzmax = min([zmax,zsize-1])
 tzsize = tzmax-tzmin+1  
 ;print,zmin,zmax,tzmin,tzmax,tzsize,zsize
 ;subsample by factor 2
 d2 = rebin(d[*,*,tzmin:tzmax],xsize*2,ysize*2,tzsize*2,/sample)
 
 ;convolve with laplacian (detect +ve residuals)
 l2p = convol(d2,lap,/center,/edge_wrap)
 l2n = -l2p

 ;set negative values to zero
 r_l2 = where(l2p lt 0,i_l2)
 if (i_l2 gt 0) then l2p[r_l2] = 0

 ;resample to original size
 tlp = rebin(l2p,xsize,ysize,tzsize)
 lp[*,*,tzmin+1:tzmax-1] = tlp[*,*,1:tzsize-2]

 if (not keyword_set(noneg)) then begin
 ;set negative values to zero
     r_l2 = where(l2n lt 0,i_l2)
     if (i_l2 gt 0) then l2n[r_l2] = 0

 ;resample to original size
     tln = rebin(l2n,xsize,ysize,tzsize) 
     ln[*,*,tzmin+1:tzmax-1] = tln[*,*,1:tzsize-2]
 endif
 
endwhile ; end of sectioning

;finished making +ve & -ve Laplacians

;calculate ratio S
spp = lp / (0.5*nc)

;calculate ratio L+/f
lp_f = lp/f

;identify & replace bad pixels
; +ve bad pixels
r_bpp = where(lp_f gt f_lim and spp gt s_lim,i_bpp)
if (i_bpp gt 0) then begin
    badpix[r_bpp] = badpix[r_bpp]+d_out[r_bpp]-dm3[r_bpp]
    d_out[r_bpp] = dm3[r_bpp]
endif
print,'corrected '+strtrim(string(i_bpp),2)+' +ve pixels'

r_bpn=-1
i_bpn=0
 if (not keyword_set(noneg)) then begin
 ;calculate ratio S
     snp = ln / (0.5*nc)

 ;calculate ratio L+/f
     ln_f = ln/f

 ;identify & replace bad pixels
 ; -ve bad pixels
     r_bpn = where(ln_f gt f_lim and snp gt s_lim,i_bpn)
     if (i_bpn gt 0) then begin
         badpix[r_bpn] = badpix[r_bpn]+d_out[r_bpn]-dm3[r_bpn]
         d_out[r_bpn] = dm3[r_bpn]
     endif
     print,'corrected '+strtrim(string(i_bpn),2)+' -ve pixels'

 endif

new_bad = i_bpp + i_bpn
r_bp_tot = r_bp_tot + new_bad
d = d_out
loop = loop + 1

print,'corrected a total of ',+strtrim(string(r_bp_tot),2)+' pixels'

endwhile
;finished the loop

frac_str = string(format='(f5.2)',100.*r_bp_tot/(xsize*ysize*zsize))
print,'corrected '+frac_str+'% of all pixels in cube'

badpixn[where(badpix ne 0.)] = 1

; replace structure into cube
d_out = temporary(d_out) + dxy7 + dz49

; replace NaN pixels in output cube
if (nan_i gt 0) then d_out[nan_r] = !values.f_nan

;replace excised slices
tmp = fltarr(xsize0,ysize0,zsize0) + !values.f_nan
tmp[n1:n2,n3:n4,n5:n6] = d_out
d_out = tmp
writefits,strmid(cube,0,strlen(cube)-5)+'_x.fits',d_out,hdr

;finish
print,'took '+strtrim(string(fix(0.5+systime(1)-t_start)),2)+' sec'

;write out extra cubes as required
if keyword_set(wbadpix) then begin
    tmp[n1:n2,n3:n4,n5:n6] = badpixn
    badpixn = tmp
    writefits,strmid(cube,0,strlen(cube)-5)+'_b.fits',badpixn,hdr
endif

;;MF this crashes...
;nc_name = strmid(cube,0,strlen(cube)-5)+'_n.fits'
;nc_write = file_test(nc_name)
;if (nc_write eq 1 and keyword_set(wnoise)) then begin
;    print,'over-writing noise cube '+nc_name
;    tmp[n1:n2,n3:n4,n5:n6] = nc
;    nc = tmp
;    writefits,nc_name,nc,hdr
;endif
;if (nc_write eq 0) then begin
;    print,'writing noise cube '+nc_name
;    tmp[n1:n2,n3:n4,n5:n6] = nc
;    nc = tmp
;    writefits,nc_name,nc,hdr
;endif
;if (nc_write eq 1 and not(keyword_set(wnoise))) then begin
;    print,'will not over-write noise cube '+nc_name
;endif

if keyword_set(test) then begin
print,'Writing out extra files'
    writefits,'testlp.fits',lp,hdr
    writefits,'testln.fits',ln,hdr
    writefits,'testnc.fits',nc,hdr
    writefits,'testspp.fits',spp,hdr
    writefits,'testsnp.fits',snp,hdr
    writefits,'testf.fits',f,hdr
    writefits,'testlpf.fits',lp_f,hdr
    writefits,'testlnf.fits',ln_f,hdr
    writefits,'testbp.fits',badpix,hdr
    writefits,'testbpn.fits',badpixn,hdr
    writefits,'testdz49.fits',dz49,hdr
    writefits,'testdxy7.fits',dxy7,hdr
    writefits,'testd_raw.fits',d_raw,hdr

openw,1,'bp_'+strmid(cube,0,strlen(cube)-5)+'.txt'
if (i_bpp gt 0) then begin
    pos = array_indices(badpix,r_bpp)
    for i=long(0),n_elements(r_bpp)-1 do printf,1,pos[0,i]+1,pos[1,i]+1,pos[2,i]+1,badpix[r_bpp[i]]
endif
if (i_bpn gt 0) then begin
    pos = array_indices(badpix,r_bpp)
    for i=long(0),n_elements(r_bpp)-1 do printf,1,pos[0,i]+1,pos[1,i]+1,pos[2,i]+1,badpix[r_bpp[i]]
endif
close,1

endif

;stop
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
