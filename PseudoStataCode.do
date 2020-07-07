/* Pseudo-code for pooling study-specific non-linear relationships
 using fractional polynomial: a step-by-step guide.
 We use three studies in this code; study 1, 2, and 3*/

/*Clear all*/
clear

/*Change working directory*/
cd "C:\xxx\yyy"

/*Read the data*/
cd "C:\Users\uqddarss\OneDrive\1.SPH\Research\MetaAnalysisMethods\Output"
import delimited hypodat.csv, clear 

/*Step 1: Exploratory analysis*/
*Description of data;
ds
describe
tabulate study
summarize bmi menopauseage age
tab smoking 
tab nochild 
tab education
*histograms;
histogram(bmi), by(study)
histogram(menopauseage), by(study)
*density;
*run next two lines together;
kdensity bmi if study==1, addplot(kdensity bmi if study==2 || ///
kdensity bmi if study==3) legend(ring(0) pos(2) label(1 "Study 1") label(2 "Study 2") label(3 "Study 3"))
*run next two lines together;
kdensity menopauseage if study==1, addplot(kdensity menopauseage if study==2 || ///
kdensity menopauseage if study==3) legend(ring(0) pos(2) label(1 "Study 1") label(2 "Study 2") label(3 "Study 3"))

/*Step 2: Assessing the linearity of the relationship*/
*scatter;
scatter menopauseage bmi, by(study)
scatter menopauseage bmi if study==1

/*Step 3: Modelling the effect of confounding variables on the explanatory variable*/
/*Confounder  model*/
/*The variables to use in the confounder model: age,
smoking, no of children, education*/
regress bmi i.smoking i.nochild i.education age if study==1
predict m1 if study==1
regress bmi i.smoking i.nochild i.education age if study==2
predict m2 if study==2
regress bmi i.smoking i.nochild i.education age if study==3
predict m3 if study==3
/*Combine m1, m2 and m3 to m*/
if m3 != .{
gen m = m3
}
replace m = m2 if m2 != .
replace m = m1 if m1 != .
drop m1-m3
/*Confounder variable model index, i.e. the residuals, is used to adjust for confounders in our main model*/
gen CVMI = bmi - m

/*Step 4: Selecting the best study-specific fractional polynomials for the association between the outcome and exposure variables*/
*Please note the codes below need to be run together; otherwise they won't work, so study 1 code as one go;

*Study 1;
fp <bmi>, center scale dimension(4) replace: regress menopauseage <bmi> CVMI if study==1
display "`e(fp_center_mean)'"
display "`e(fp_scale_a)'"
display "`e(fp_scale_b)'"

/*Based on the model comparison output (Deviance and p-values) select the most suitable model. 
In addition plot all, linear, FP1 to FP4 to see the difference*/
regress menopauseage bmi CVMI if study==1
predict fit1
label variable fit1 "Linear"
fp <bmi>, center scale fp(3) replace: regress menopauseage <bmi> CVMI if study==1
matrix list e(fp_fp)
predict fit2
label variable fit2 "FP1(3)"
fp <bmi>, center scale fp(3 3) replace: regress menopauseage <bmi> CVMI if study==1
matrix list e(fp_fp)
predict fit3
label variable fit3 "FP2(3 3)"
fp <bmi>, center scale fp(3 3 3) replace: regress menopauseage <bmi> CVMI if study==1
matrix list e(fp_fp)
predict fit4
label variable fit4 "FP3(3 3 3)"
fp <bmi>, center scale fp(2 2 3 3) replace: regress menopauseage <bmi> CVMI if study==1
matrix list e(fp_fp)
predict fit5
label variable fit5 "FP4(2 2 3 3)"

scatter menopauseage fit1 fit2 fit3 fit4 fit5 bmi if study==1, c(. l l l l l) m(o i i i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit1 fit2 fit3 fit4 bmi if study==1, c(. l l l l) m(o i i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit2 fit3 fit4 bmi if study==1, c(. l l l) m(o i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit2 fit3 bmi if study==1, c(. l l) m(o i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit1 bmi if study==1, c(. l) m(o i) msize(small) lpattern(. -_.) ytitle("Menopauseage") xtitle("BMI") sort
line fit2 bmi, sort
*For illustration purpose of this code the model choice for study 1 is FP2(3 3);
drop fit1-fit5

*Study 2;
fp <bmi>, center scale dimension(4) replace: regress menopauseage <bmi> CVMI if study==2
display "`e(fp_center_mean)'"
display "`e(fp_scale_a)'"
display "`e(fp_scale_b)'"

/*Based on the model comparison output (Deviance and p-values) select the most suitable model. 
In addition plot all, linear, FP1 to FP4 to see the difference*/
regress menopauseage bmi CVMI if study==2
predict fit1
label variable fit1 "Linear"
fp <bmi>, center scale fp(3) replace: regress menopauseage <bmi> CVMI if study==2
matrix list e(fp_fp)
predict fit2
label variable fit2 "FP1(3)"
fp <bmi>, center scale fp(-2 -1) replace: regress menopauseage <bmi> CVMI if study==2
matrix list e(fp_fp)
predict fit3
label variable fit3 "FP2(-2 -1)"
fp <bmi>, center scale fp(0 0.5 0.5) replace: regress menopauseage <bmi> CVMI if study==2
matrix list e(fp_fp)
predict fit4
label variable fit4 "FP3(0 0.5 0.5)"
fp <bmi>, center scale fp(3 3 3 3) replace: regress menopauseage <bmi> CVMI if study==2
matrix list e(fp_fp)
predict fit5
label variable fit5 "FP4(3 3 3 3)"

scatter menopauseage fit1 fit2 fit3 fit4 fit5 bmi if study==2, c(. l l l l l) m(o i i i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit1 fit2 fit3 fit4 bmi if study==2, c(. l l l l) m(o i i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit2 fit3 fit4 bmi if study==2, c(. l l l) m(o i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit2 fit3 bmi if study==2, c(. l l) m(o i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit2 bmi if study==2, c(. l) m(o i) msize(small) lpattern(. -_.) ytitle("Menopauseage") xtitle("BMI") sort
line fit2 bmi, sort
*For illustration purpose of this code the model choice for study 2 is FP1(3);
drop fit1-fit5

*Study 3;
fp <bmi>, center scale dimension(4) replace: regress menopauseage <bmi> CVMI if study==3
display "`e(fp_center_mean)'"
display "`e(fp_scale_a)'"
display "`e(fp_scale_b)'"

/*Based on the model comparison output (Deviance and p-values) select the most suitable model. 
In addition plot all, linear, FP1 to FP4 to see the difference*/
regress menopauseage bmi CVMI if study==3
predict fit1
label variable fit1 "Linear"
fp <bmi>, center scale fp(2) replace: regress menopauseage <bmi> CVMI if study==3
matrix list e(fp_fp)
predict fit2
label variable fit2 "FP1(2)"
fp <bmi>, center scale fp(0.5 2) replace: regress menopauseage <bmi> CVMI if study==3
matrix list e(fp_fp)
predict fit3
label variable fit3 "FP2(0.5 2)"
fp <bmi>, center scale fp(0.5 3 3) replace: regress menopauseage <bmi> CVMI if study==3
matrix list e(fp_fp)
predict fit4
label variable fit4 "FP3(0.5 3 3)"
fp <bmi>, center scale fp(3 3 3 3) replace: regress menopauseage <bmi> CVMI if study==3
matrix list e(fp_fp)
predict fit5
label variable fit5 "FP4(3 3 3 3)"

scatter menopauseage fit1 fit2 fit3 fit4 fit5 bmi if study==3, c(. l l l l l) m(o i i i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit1 fit2 fit3 fit4 bmi if study==3, c(. l l l l) m(o i i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit2 fit3 fit4 bmi if study==3, c(. l l l) m(o i i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit2 fit3 bmi if study==3, c(. l l) m(o i i) msize(small) ytitle("Menopauseage") xtitle("BMI") sort
scatter menopauseage fit2 bmi if study==3, c(. l) m(o i) msize(small) lpattern(. -_.) ytitle("Menopauseage") xtitle("BMI") sort
line fit3 bmi, sort
*For illustration purpose of this code the model choice for study 2 is FP3(0.5 3 3);
drop fit1-fit5


/*Step 5: Estimating the functional forms for the association between the outcome and exposure variables for each study, adjusting for confounders*/
*Please note the codes below need to be run together, otherwise they won't work, so study 1 code as one go;
/***Please install xpredict package from stata****/

/*Install the stata package xpredict*/

/*Study 1 - The model choice for study 2 is FP2(3,3)*/
fp <bmi>, center scale fp(3 3) replace: regress menopauseage <bmi> CVMI if study==1
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
regress menopauseage `v' CVMI if study==1
*xpredict me1, with(`v')
*xpredict me1_con, with(`v') cons
xpredict me1_CVMI, with(`v' CVMI) cons
/*Standard error and variance of the estimated effect partial prediction*/
xpredict seme1, with(`v' CVMI) cons stdp
gen vme1 = seme1^2

/*Study 2 - The model choice for study 2 is FP1(3)*/
fp <bmi>, center scale fp(3) replace: regress menopauseage <bmi> CVMI if study==2
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

/*Study 3 - The model choice for study 3 is FP3(0.5 3 3)*/
fp <bmi>, center scale fp(0.5 3 3) replace: regress menopauseage <bmi> CVMI if study==3
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
regress menopauseage `v' CVMI if study==3
*xpredict me3, with(`v')
*xpredict me3_con, with(`v') cons
xpredict me3_CVMI, with(`v' CVMI) cons
/*Standard error and variance of the estimated effect partial prediction*/
xpredict seme3, with(`v' CVMI) cons stdp
gen vme3 = seme3^2


/*Step 6: Pooling the functional forms across studies */

/*Fixed effect weights*/
gen inv_vme1 = cond(seme1==0,.,1/vme1)
gen inv_vme2 = cond(seme2==0,.,1/vme2)
gen inv_vme3 = cond(seme3==0,.,1/vme3)

/*Sum of the fixed effect weights*/
gen sumfixw = inv_vme1 + inv_vme2 + inv_vme3

/*Standardised fixed effect weights*/
gen sw1 = inv_vme1/sumfixw
gen sw2 = inv_vme2/sumfixw
gen sw3 = inv_vme3/sumfixw


/*Overall fixed-effect estimate and the variance*/
gen fixme =  sw1*me1_CVMI + sw2*me2_CVMI + sw3*me3_CVMI
gen varfixme = cond(sumfixw==0,.,1/sumfixw)

/*Random effect weights*/

gen Q = inv_vme1*(me1_CVMI - fixme)^2 + inv_vme2*(me2_CVMI - fixme)^2 + inv_vme3*(me3_CVMI - fixme)^2

/*Sum of inverse variance squared*/
gen suminv_vsq = inv_vme1^2 + inv_vme2^2 + inv_vme3^2

/*s squared*/
gen tausq = cond(sumfixw==0,0,max(0,(Q-3+1)/(sumfixw - (suminv_vsq/sumfixw))))

/*sum of random-effect weights*/
gen sumrandw = (1/(inv_vme1 + tausq)) + (1/(inv_vme2 + tausq)) + (1/(inv_vme3 + tausq))

*Standardised random effect weight;
gen srw1 = ((inv_vme1 + tausq)^(-1))/sumrandw
gen srw2 = ((inv_vme2 + tausq)^(-1))/sumrandw
gen srw3 = ((inv_vme3 + tausq)^(-1))/sumrandw

/*Overall random-effect estimate and the variance*/
gen ranme = srw1*me1_CVMI + srw2*me2_CVMI + srw3*me3_CVMI
gen varrandme = cond(sumrandw==0,.,1/sumrandw)

/*Plots*/
twoway (scatter menopauseage bmi) ///
       (lowess me1_CVMI bmi, legend(off)), /// 
       saving(study1, replace) ytitle(Menopause Age)
	   
twoway (scatter menopauseage bmi) ///
       (lowess me2_CVMI bmi, legend(off)), /// 
       saving(study2, replace) ytitle(Menopause Age)
	   
twoway (scatter menopauseage bmi) ///
       (lowess me3_CVMI bmi, legend(off)), /// 
       saving(study3, replace) ytitle(Menopause Age)

gr combine study1.gph study2.gph study3.gph, col(2) iscale(0.5)

*pooled curve on scatter	   
twoway (scatter menopauseage bmi) ///
       (lowess ranme bmi, legend(off)), /// 
       ytitle(Menopause Age)


/*Just plot to see how things work*/
line ranme bmi, sort
line varrandme bmi, sort
















































