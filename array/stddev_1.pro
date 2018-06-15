;performs stddev of 2d array to return
; along the axis

function stddev_1, imgh
n_imgh=(size(imgh, /dim))[1]
djs_means=total(imgh, 2)/n_imgh
sqrs=total(imgh^2, 2)/n_imgh
return, sqrt(sqrs-djs_means^2)
end
