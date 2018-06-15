
;+
;
;
;   Return the mean and std of a histogram that represents 
;   the realizations of a discrete distribution
;
;
;-


function discrete_stats, histo

  ;;find elements
  disc=find_different(histo,neach=neach)
  
  ;;build probability 
  prob_each=neach/n_elements(histo)

  ;;find mean 
  mean=total(disc*prob_each)

  ;;find std 
  std=sqrt(total((disc-mean)^2*prob_each))

  return, [mean, std]

end
