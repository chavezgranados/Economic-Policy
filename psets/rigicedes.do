*********** RIGIDECES EN PRECIOS **************
Luis Chávez, 2025

*Fuente de datos: BCRP, PBI 2007 y PBI: 1990-2024

ed
gen d= y_n/ y_r
ed
drop if t<1990
ed
gen ly_r=ln(y_r)
gen ly_n=ln(y_n)
regress y_r t
predict y_rhat, xb
twoway (scatter t y_rhat)
twoway (scatter t y_rhat) (scatter t y_r)
twoway (line t y_rhat) (line t y_r)
twoway (line y_rhat t) (line y_r t)
ed
ed
drop if t>2024
ed
compress
ed
gen dly_n= ly_n- y_rhat
gen dly_r= ly_r- y_rhat
regress dly_n dly_r
regress d dly_n
regress dly_r dly_n
ed
regress d dly_n if t<2007
regress d dly_n if t>=2007
h regress postestimation
estat ic
regress d dly_n if t<2007
estat ic
regress d dly_n if t<=2000
regress d dly_n if t<=2019
regress d dly_n if t>=2000
twoway (line y_rhat t) (line y_r t)
doedit