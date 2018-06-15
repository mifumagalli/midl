;Adapted From the muse pipeline muse_wcs.c

pro muse_applywcs, Xpos, Ypos, aRA, aDEC, RApos, DECpos
    
    CPL_MATH_DEG_RAD = 57.29577951308232087679815
    CPL_MATH_PI_2 = 1.57079632679489661923132
    
    dp = aDEC / CPL_MATH_DEG_RAD
    phi = temporary(Xpos)
    theta = temporary(Ypos) + CPL_MATH_PI_2
    RApos = (atan(cos(theta) * sin(phi), sin(theta) * cos(dp) + cos(theta) * sin(dp) * cos(phi)) * CPL_MATH_DEG_RAD) +aRA
    DECpos = (asin(sin(theta) * sin(dp) - cos(theta) * cos(dp) * cos(phi)) * CPL_MATH_DEG_RAD)

end
