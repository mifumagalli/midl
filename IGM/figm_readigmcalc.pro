;Procedure that reads the  binary file from figm_calculator
;and computes the IGM transmission.
;Returns observed lambda (AA) and IGM transmission for a given
;redshift. 


;filename --> name of binary file output of figm_calculator
;lam_out  --> output: observed lambda (AA)
;trans_out -> output: IGM transmission  
;XRANGE    -> keyword: xrange in the plot   
;xlog      -> keyword: set xlog for the plot   
;noplot    -> keyword: disable dispaly option  

;Written by M.F. on Oct 2009
;Add boundary values on Apr 2010
;Report bugs to mfumagalli@ucolick.org




PRO figm_readigmcalc, filename, lam_out, trans_out, XRANGE=xrange, xlog=xlog, noplot=noplot
  
  
;open bianray to read
  openr,lun,filename,/GET_LUN
  
  
;get num line
  p=STRPOS(filename,".dat")
  p2=STRPOS(filename,"_c")
  num=LONG(STRMID(filename,p2+2,p-p2-2))
  
  lambda=DBLARR(num)
  trans=DBLARR(num)
  
  for i=0L,num-1 do begin 
     lambda[i]=read_binary(lun,DATA_TYPE=5,DATA_DIMS=1)
     trans[i]=read_binary(lun,DATA_TYPE=5,DATA_DIMS=1)
  endfor
  


  !x.style=1


  if ~keyword_set(XRANGE) then xrange=[min(lambda),max(lambda)]
  
  indx=where(trans NE 0 and trans NE 1, nind)
  if(nind GT 0) then trans[indx]=exp(-trans[indx])
  
  lam_out=lambda
  trans_out=trans
  
  ;add a bunch of zero and 1 left and right
  lam_out=[lambda[0]-10.,lambda[0]-5.,lam_out,lambda[num-1]+5.,lambda[num-1]+10.]
  trans_out=[0.,0.,trans,1.,1.]

  if ~keyword_set(noplot) then begin
     if keyword_set(XLOG) then plot, lambda, trans, xrange=xrange, /xlog $
     else plot, lambda, trans, xrange=xrange
  endif
  
  
  free_lun, lun

END
