/* Pseudo-code for pooling study-specific non-linear realathinships
 using fractional ploynomial: a step-by-step guide.
 We use two studies in this code; study 1 and study 2*/

/*Clear all*/
clear

/*Confounder covariate model*/
/*The varaibles to use in the confounder model: age,
smoking, no of children, education*/
regress bmi i.smoking i.child i.education age if study==1
predict m1 if study==1
regress bmi i.smoking i.child i.education age if study==2
predict m2 if study==2


/*Combine m1 and m2 to m*/
if m1 != .{
gen m = m1
}
replace m = m2 if m2 != .
drop m1-m2

/*Confounder covariate index is to adjusted for in our main model*/
gen CVMI = bmi - m

/*Find the model choice for study 1 - Step 4 in the paper*/
*Please note the codes below need to be run together, otherwise it won't work, so study 1 code as one go;
fp <bmi>, center scale dimension(4) replace: regress menopauseage <bmi> CVMI if study==1
display "`e(fp_center_mean)'"
display "`e(fp_scale_a)'"
display "`e(fp_scale_b)'"

*Based on the model comprison output FP1 with power 3 seems to be the best. We plot all, liner, FP1 to FP4 to see the difference;
regress menopauseage bmi CVMI if study==1
predict fit1
label variable fit1 "Linear"
fp <bmi>, center scale fp(3) replace: regress menopauseage <bmi> CVMI if study==1
matrix list e(fp_fp)
predict fit2
label variable fit2 "FP1(3)"
fp <bmi>, center scale fp(-2 2) replace: regress menopauseage <bmi> CVMI if study==1
matrix list e(fp_fp)
predict fit3
label variable fit3 "FP2(-2 2)"
fp <bmi>, center scale fp(-2 -2 3) replace: regress menopauseage <bmi> CVMI if study==1
matrix list e(fp_fp)
predict fit4
label variable fit4 "FP3(-2 -2 3)"
fp <bmi>, center scale fp(2 3 3 3) replace: regress menopauseage <bmi> CVMI if study==1
matrix list e(fp_fp)
predict fit5
label variable fit5 "FP4(2 3 3 3)"

scatter menopauseage fit1 fit2 fit3 fit4 fit5 bmi if study==1, c(. l l l l l) m(o i i i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit1 fit2 fit3 fit4 bmi if study==1, c(. l l l l) m(o i i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit2 fit3 fit4 bmi if study==1, c(. l l l) m(o i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit2 fit3 bmi if study==1, c(. l l) m(o i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit1 bmi if study==1, c(. l) m(o i) msize(small) lpattern(. -_.) ytitle("Menopauseage") xtitle("BMI") sort
line fit3 bmi, sort

*The model choice for study 1 is FP2(-2,2);
drop fit1-fit5

/*Find the model choice for study 2 - Step 4 in the paper Darssan et al*/
*Please note the codes below need to be run together, otherwise it won't work, so study 2 code as one go;
fp <bmi>, center scale dimension(4) replace: regress menopauseage <bmi> CVMI if study==2
display "`e(fp_center_mean)'"
display "`e(fp_scale_a)'"
display "`e(fp_scale_b)'"

*Based on the model comprison output FP1(0) seems to be the best. We plot all, FP1 to FP4 to see the difference;
fp <bmi>, center scale fp(0) replace: regress menopauseage <bmi> CVMI if study==2
matrix list e(fp_fp)
predict fit1
label variable fit1 "FP1(0)"
fp <bmi>, center scale fp(1 3) replace: regress menopauseage <bmi> CVMI if study==2
matrix list e(fp_fp)
predict fit2
label variable fit2 "FP2(1 3)"
fp <bmi>, center scale fp(-2 3 3) replace: regress menopauseage <bmi> CVMI if study==2
matrix list e(fp_fp)
predict fit3
label variable fit3 "FP3(-2 3 3)"
fp <bmi>, center scale fp(3 3 3 3) replace: regress menopauseage <bmi> if study==2
matrix list e(fp_fp)
predict fit4
label variable fit4 "FP4(3 3 3 3)"
scatter menopauseage fit1 fit2 fit3 fit4 bmi if study==2, c(. l l l l) m(o i i i i) msize(small) ytitle("Menopauseage") xtitle("BMI")  sort
scatter menopauseage fit1 fit2 fit3 bmi if study==2, c(. l l l) m(o i i i) msize(small) ytitle("Menopauseage") xtitle("BMI")  sort
line fit3 fit2 bmi, sort
*The model choice for study 2 is FP2(1,3);
drop fit1-fit4

/*Estimate the functional form - Step 5 in the paper Darssan et al*/
*Study 1;
local eqxb xb 
marksample touse, novarlist
frac_wgt `"`exp'"' `touse' `"`weight'"'
local wgt `r(wgt)'
fp <bmi>, center scale fp(-2 2) replace: regress menopauseage <bmi> CVMI `wgt' if study==1 & `touse'==1
matrix list e(fp_fp)
display "`e(fp_center_mean)'"
display "`e(fp_scale_a)'"
display "`e(fp_scale_b)'"
local v `e(fp_terms)'
drop `e(fp_terms)'
`=e(fp_gen_cmdline)'
display "`e(fp_gen_cmdline)'"
display "`e(fp_terms)'"
display "`v'"
regress menopauseage `v' CVMI `wgt' if study==1 & `touse'==1
*xpredict me1, with(`v') double `eqxb'
*xpredict me1_con, with(`v') cons double `eqxb'
xpredict me1_CVMI, with(`v' CVMI) cons double `eqxb'
/*Standard error and variance of the estimated effect - prediction*/
xpredict seme1, with(`v' CVMI) cons stdp
gen vme1 = seme1^2
/*Fixed effect weights*/
gen inv_vme1 = cond(seme1==0,.,1/vme1)

/*
/*Just plot to see how things work*/
line me1 CVMI, sort
line me1_con bmi, sort
line me1_CVMI bmi, sort
*/
//drop me1 me1_con me1_CVMI seme1 vme1 inv_vme1

/*Study 2 - The model choice for study 2 is FP2(1,3)*/
fp <bmi>, center scale fp(1 3) replace: regress menopauseage <bmi> CVMI if study==2
matrix list e(fp_fp)
display "`e(fp_center_mean)'"
display "`e(fp_scale_a)'"
display "`e(fp_scale_b)'"
local v `e(fp_terms)'
drop `e(fp_terms)'
`=e(fp_gen_cmdline)'
display "`e(fp_gen_cmdline)'"
display "`e(fp_terms)'"
display "`v'"
regress menopauseage `v' CVMI if study==2
*xpredict me2, with(`v')
*xpredict me2_con, with(`v') cons
xpredict me2_CVMI, with(`v' CVMI) cons
/*Standard error and variance of the estimated effect partial prediction*/
xpredict seme2, with(`v' CVMI) cons stdp
gen vme2 = seme2^2
/*Fixed effect weights*/
gen inv_vme2 = cond(seme2==0,.,1/vme2)
/*
/*Just plot to see how things work*/
line me2 bmi, sort
line me2_con bmi, sort
line me2_CVMI bmi, sort
*/
//drop me2 me2_con me2_CVMI seme2 vme2 inv_vme2


/*Fixed effect weights*/

/*Sum of the fixed effect weights*/
gen sumfixw = inv_vme1 + inv_vme2

/*Standardised fixed effect weights*/
gen sw1 = inv_vme1/sumfixw
gen sw2 = inv_vme2/sumfixw


/*Overall fixed-effect estimate and the variance*/
gen fixme =  sw1*me1_CVMI + sw2*me2_CVMI
gen varfixme = cond(sumfixw==0,.,1/sumfixw)

/*
/*Just plot to see how things work*/
line fixme bmi, sort
line varfixme bmi, sort
*/

/*Random effect weights*/

/*Q*/
gen Q = inv_vme1*(me1_CVMI - fixme)^2 + inv_vme2*(me2_CVMI - fixme)^2

/*Sum of inverse variance squared*/
gen suminv_vsq = inv_vme1^2 + inv_vme2^2

/*Tau squared*/
gen tausq = cond(sumfixw==0,0,max(0,(Q-2+1)/(sumfixw - (suminv_vsq/sumfixw))))

/*sum of random-effect weights*/
gen sumrandw = (1/(inv_vme1 + tausq)) + (1/(inv_vme2 + tausq))

*Standardised random effect weight;
gen srw1 = ((inv_vme1 + tausq)^(-1))/sumrandw
gen srw2 = ((inv_vme2 + tausq)^(-1))/sumrandw

/*Overall random-effect estimate and the variance*/
gen ranme = srw1*me1_CVMI + srw2*me2_CVMI 
gen varrandme = cond(sumrandw==0,.,1/sumrandw)

/*Just plot to see how things work*/
line ranme bmi, sort
line varrandme bmi, sort


/*A plot of confounder covariate index by study ID*/
line CVMI bmi if study==1, sort















































