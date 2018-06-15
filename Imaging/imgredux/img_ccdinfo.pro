;+
;
;  Contains the common information that are inizialized by the
;  pipeline 
;
;
;-

pro img_ccdinfo, instr, path, name, side

;;initialise a common with information on the size of the ccd
common ccdinfo, xfull,yfull,nchip,namp,biaslev,satur, goodata 
;;xfull,yfull the full size of the mosaic
;;nchip the number of chip that makes the mosaic
;;namp is the number of amplifier
;;biaslev is the typical counts in the bias (use to select good flats)
;;satur  is the saturation level (use to reject saturated flats and CR)
;;goodata is array of x0,x1,y0,y1 which specify the region of good data in case of vignetting (to normalise flat,wgt image etc..)
if(instr eq 'LRIS') then begin
   head=headfits(path+name[0],exten=0)
   date=fxpar(head,"DATE")
   xfull=4096
   yfull=4096
   nchip=2
   namp=4
   biaslev=1000 ;;blue and red side
   satur=65535
   ;;set good region
   if(date ge '2013-01-01') then begin 
  
      if(side[0] eq 'B') then begin
         goodata=[220,3640,430,2830] 
      endif else goodata=[430,3780,640,3085]
  
   endif 
   if(date lt '2013-01-01' and date ge '2009-04-01') then begin
      
      if(side[0] eq 'B') then begin
         goodata=[350,3750,750,3154] 
      endif else goodata=[450,3780,665,3085]
      
   endif
   if(date lt '2009-04-01') then begin
 
      if(side[0] eq 'B') then begin
         goodata=[335,3740,950,3350] 
      endif else stop           ;use LRISr
      
   endif
   
endif 
;;make this an if with date and remove LRISr
if(instr eq 'LRISr') then begin
    ;;old red side on lris - single chip
    xfull=2048
    yfull=2048
    nchip=1
    namp=2
    biaslev=850
    satur=65535
    goodata=[190,1650,0,2006]
endif 
if(instr eq 'LBC') then begin
    xfull=6292
    yfull=6731
    nchip=4
    namp=4
    biaslev=700                    ;;blue and red 750,200. Use 700 (conservative flat rejection)
    satur=65535
    goodata=[0,xfull-1,0,yfull-1]  ;;take full mosaic
endif 
if(instr eq 'ESI') then begin
    xfull=1100
    yfull=1550
    nchip=1
    namp=2
    biaslev=1000 
    satur=65535
    goodata=[200,860,100,1380]
endif 
if(instr eq 'FEINC') then begin   ;;binned 2x2
    xfull=1024
    yfull=1024
    nchip=1
    namp=1        ;;treat all amps together
    biaslev=1000 
    satur=65535
    goodata=[1,1000,1,1023]
endif
if(instr eq 'IMACS') then begin
    xfull=8192
    yfull=8192
    nchip=8
    namp=8
    biaslev=2000  ;;mean value 
    satur=65000
    goodata=[0,xfull-1,0,yfull-1]  ;;take full mosaic
endif 


if(instr ne 'LRIS' and instr ne 'LBC' and instr ne 'ESI' and $
   instr ne 'LRISr' and instr ne 'FEINC' and instr ne 'IMACS') then begin
    splog, 'Instrument ', instr, ' not supported!'
    return
endif


end
