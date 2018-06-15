;+
;
;  Check the velocity range in the .clm file
;
;
;-


pro check_velocity, clmfile, path=path, list=list, xrange=xrange

  if ~keyword_set(path) then path='./'
  if ~keyword_set(list) then list='dla.lst'


  ;;open file to read
  openr, lun, path+clmfile, /get_lun

  ;;get name
  objname=""
  readf, lun, objname
  
  ;;get flag
  flag=0.
  readf, lun, flag

  ;;get spect file
  specfile=""
  readf, lun, specfile

  ;;get z
  redsh=0d
  readf, lun, redsh

  ;;get tables
  tab=""
  readf, lun, tab

  ;;get colum
  colum=""
  readf, lun, colum
  
  ;;get zero
  zero=0
  readf, lun, zero

  ;;open spec
  spec=x_readspec(specfile,0,INFLG=1,WAV=wav)


  ;;open atom list
  readcol, getenv("XIDL_DIR")+"/Spec/Lines/Lists/"+list, $
          lambda, name, lshort, format="D,A,I"

  ;;open ps
  !x.style=1
  !y.style=1
  m_psopen, "check_"+clmfile+".ps", /land

  
  while not EOF(lun) do begin

      ;;identify current transition

      flag=0
      line=""

      readf, lun, flag
      readf, lun, line
      
      values=1d*strsplit(line,",", /extr) ;;[0] lambda, [1,2] +/- vel

      ;;match lambda
      thisl=where(abs(lambda - values[0]) lt 0.001,nthis)
      
      ;;label 
      if(nthis eq 1) then title=name[thisl]+" "+string(lshort[thisl]) else $
              title=string(values[0])


      
      ;;flag
      if(flag eq 0 or flag eq 2 or flag eq 4) then title=title+" non use"
      if(flag eq 1 or flag eq 3 or flag eq 5) then title=title+" use"

      ;;limits
      if(flag eq 2 or flag eq 3) then title=title+" >"
      if(flag eq 4 or flag eq 5) then title=title+" <"


      ;;find observed range
      obsl=values[0]*(1+redsh)
      velocity=(wav-obsl)*2.99792d5/obsl

      
      ;;limit
      ;if ~keyword_set(xrange) then begin
         if(values[1] lt 0 and values[2] gt 0) then xrange=3.*[values[1],values[2]]
         if(values[1] gt 0) then xrange=[0.3*values[1],3*values[2]]
         if(values[2] lt 0) then xrange=[3*values[1],0.3*values[2]]
      ;endif

      ;;plot
      plot, velocity, spec, xrange=xrange, psym=10,$
              title=title, yrange=[-0.2,1.2]
      oplot, [-1D5,1D5], [1.,1.], color=fsc_color("blue"), line=2
      oplot, [-1D5,1D5], [0.,0.], color=fsc_color("blue"), line=2
      oplot, [0,0], [-10,10], line=1
      oplot, [values[1],values[1]], [-10,10], line=3, color=fsc_color("red")
      oplot, [values[2],values[2]], [-10,10], line=3, color=fsc_color("red")
      
      
  endwhile
  
  
  ;;close
  free_lun, lun
  m_psclose

end
