;This procedure computes the frequentist probability that one object
;is an interloper according to number counts statistics.

;magnitude       --> magnitude of the object in the same 
;                    system and filter provided in  cumulat_file 
;bimpact         --> distance of the object in same units as in the cumulative distribution.
;cumulat_file    --> File with the cumulative distribution
;                    of the background sources two colums: mag and
;                    counts (e.g. cnt/arcsec^2). 


;default file is the number counts in the U band from Grazian et al., 2009




PRO prob_frequen, magnitude, bimpact, Prob, cumulat_file


;define path file
if ~keyword_set(cumulat_file)  then $ 
  cumulat_file=getenv("MIDL")+'/Statistic/frequentist/grazian/mag_cum_grazian.dat'


;read data and interpolation
readcol, cumulat_file, mag_data, cnt_data
cnt_interpol=interpol(cnt_data,mag_data,magnitude)

;warning if extrapola
ext1=where(magnitude gt max(mag_data), next1)
ext2=where(magnitude lt min(mag_data), next2)

if(next1 gt 0) then splog, 'WARNING!: Extrapolation for ', magnitude[ext1]
if(next2 gt 0) then splog, 'WARNING!: Extrapolation for ', magnitude[ext2]

;compute absolute number
number=!PI*bimpact^2*cnt_interpol

;compute probability
Prob=1-Exp(-Number)


END
