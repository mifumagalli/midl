;+
;
;
; Compute the weighted mean, standard deviation and error on the mean
;   
;
; INPUT:
; value     array of the distribution
; weight    array of the weight
;
; OUTPUT:
; opt_mean      mean
; opt_std       standard deviation   
; opt_meanerr   error on the mean
;
;
;
;-


pro weight_stats, value, weight, opt_mean, opt_std, opt_meanerr

  
  ;;optimal mean
  opt_mean=total(value*weight)/total(weight)
  
  ;;optimal STD
  opt_std=sqrt(total(weight*(value-opt_mean)^2)/total(weight))
  
  ;;mean error 
  opt_meanerr=opt_std*sqrt(total((weight/total(weight))^2))
     
end
