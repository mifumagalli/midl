;+
;
;
;  Deal with some calibration for the 30mt telescope at Nov 2011.
;
;  http://www.iram.es/IRAMES/mainWiki/EmirforAstronomers?action=AttachFile&do=view&target=EmirCommissioningReportVers1.1.pdf
;
;
;
;   For an input frequency in GHz, return a S/Ta in Jy/K
;
;
;-


pro iramcalib, infreq, outcalib


  ;;frequency GHz 
  freq=[86.,145.,210.,260.,330.]
  ;;hpbw arcsec
  hpbw=[29.,16.,11.,9.,7.]
  ;;calibration S/Ta [Jy/K]
  calib=[5.9,6.4,7.5,8.4,12.0]

  ;;produce interpolation function 
  outcalib=interpol(calib,freq,infreq)
  ;print, outcalib 

end
