;+
;
;
;
;  Parse the ray output 
;
;
;
;
;-


pro parse_rayoutput, filename


  ;;open file for reading
  openr, lun, filename, /get_lun
  

  ;;toss header
  junk=""
  readf, lun, junk


  ;;read number of rays and box size
  nrays=0.
  boxsize=0.
  readf, lun, nrays, boxsize

  ;;parse rays 
  ray_pos1=fltarr(nrays)
  ray_pos2=fltarr(nrays)
  ray_plane=fltarr(nrays)

  buf=0.
  buf1=0.
  buf2=0.
  
  for rr=0, nrays-1 do begin

     readf, lun, buf, buf1, buf2
     ray_pos1[rr]=buf
     ray_pos2[rr]=buf1
     ray_plane[rr]=buf2
     
  endfor

  
  ;;toss second header 
  readf, lun, junk


  idb=0 
  raymemb=0 
  raylim=0 
  pos0=0 
  pos1=0 
  pos2=0 
  rvirb=0 
  mvirb=0 
  idpar=0
  line=0
  
  while (not EOF(lun)) do begin
     
     readf, lun, idb, raymemb, raylim, pos0, pos1, pos2, rvirb, mvirb, idpar
     
     if(line eq 0) then begin 
    
        id = idb
        ray_member = raymemb
        ray_limit = raylim
        position0 = pos0
        position1 = pos1
        position2 = pos2
        rvir = rvirb
        mvir = mvirb
        id_parent = idpar

     endif else begin

        id = [id,idb]
        ray_member = [ray_member,raymemb]
        ray_limit = [ray_limit,raylim]
        position0 = [position0,pos0]
        position1 = [position1,pos1]
        position2 = [position2,pos2]
        rvir = [rvir,rvirb]
        mvir = [mvir,mvirb]
        id_parent = [id_parent,idpar]

     endelse 
     
     line++
     
  endwhile
  
  close, lun
  free_lun, lun


  ;;create structure
  str={numrays:nrays,boxsize:boxsize,ray_pos1:ray_pos1,ray_pos2:ray_pos2,$
       ray_plane:ray_plane, halo_id:id, halo_raymember:ray_member,$
       halo_raylimit:ray_limit, halo_pos0:position0, halo_pos1:position1, $
       halo_pos2:position2, halo_rvir:rvir, halo_mvir:mvir, halo_idparent:id_parent}
   mwrfits, str, strmid(filename,0,strpos(filename,".txt"))+'.fits', /crea
 
end
