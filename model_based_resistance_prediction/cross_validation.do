* Cross validation script * 5 fold cross validation *

* To begin with we need to following data, firstly the agar dilution phenotypes (will be stratified on to balance cross validation)
* We also need the full model output

import delimited "agar_dilution_phenotypes.csv", encoding(ISO-8859-1) clear
rename isolate guuid
merge 1:m guuid using "model_output_for_crossval.dta"
drop if _merge != 3 /* drops the non-agar dilution models */

* Just as a reminder, this is the final model



local final "test_type c_blatem_buj2b_bin e_blatem_buj2b_bin c_blatem_buj2b_log e_blatem_buj2b_log c_blactxm_buj2be_bin e_blactxm_buj2be_bin c_blaoxa_buj2d_bin"
local final "`final' e_blaoxa_buj2d_bin c_blashv_buj2b_bin e_blashv_buj2b_bin c_blashv_buj2b_log e_blashv_buj2b_log c_blaother_bla e_blaother_bla c_temprother_sigtem"
local final "`final'  e_temprother_sigtem c_ampcprother_sigampcpr e_ampcprother_sigampcpr c_ompnon_functioning e_ompnon_functioning blatemXblactxm blactxmXblaoxa"



* mixed mic `final'  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:clsi_test eucast_test , noconstant stddev mle residuals(independent, by(test_type)) covariance(unstructured)




* Next we carve the dataset into 4, but as part of cross validation we want to obtain rougly similar groups phenotype combinations in each cross validation fold
* Therefore we partition the samples into 
* EUCAST resistant/CLSI resistant
* EUCAST Resistant/ CLSI intermediate
* EUCAST Resistant/ CLSI sensitive
* EUCAST Sensitive/CLSI sensitive


gen strata = 0
replace strata = 1 if eucast_sir == "S" & clsi_sir == "S"
replace strata = 2 if eucast_sir == "R" & clsi_sir == "S"
replace strata = 3 if eucast_sir == "R" & clsi_sir == "I"
replace strata = 4 if eucast_sir == "R" & clsi_sir == "R"

bys guuid: assert strata[1] == strata[_N]


* Just double checking strata constant across sample

* We then split each of these strata randomly into 5 folds using runiform
* Sort samples by this
* Then split the dataset into 5 (4 sets of 52 samples, 1 set of 53 samples)

/* seed set for numbers in paper, please unset/rerun if you wish to re-randomise groups, note number picked using numpy random selecting an integer between 0 and 1000000 */
set seed 348872

egen orig_tag = tag(guuid)
egen tag_samp = tag(guuid)
gen partition = runiform() if tag_samp == 1


sort guuid tag_samp
bys guuid: replace partition = partition[_N]
sort strata tag_samp partition
bys strata: replace tag_samp = sum(tag_samp)
sort strata partition guuid tag_samp
bys strata partition guuid: replace tag_samp = tag_samp[_N]

* So with this code, we generated a partition variable to randomly order elements of each strata
* we then applied this so it was there for all entries for each guuid.
* We have then numbered/sorted samples within each by this

* Next we then need to partition each strata into 5 using these values


gen sample_1 = 0
gen sample_2 = 0
gen sample_3 = 0
gen sample_4 = 0
gen sample_5 = 0

bys strata: gen strata_tot = tag_samp[_N]

* Carving up the strata
bys strata: replace sample_1 = 1 if tag_samp/strata_tot <0.2
bys strata: replace sample_2 = 1 if tag_samp/strata_tot <0.4 & tag_samp/strata_tot >= 0.2
bys strata: replace sample_3 = 1 if tag_samp/strata_tot <0.6 & tag_samp/strata_tot >=0.4
bys strata: replace sample_4 = 1 if tag_samp/strata_tot <0.8 & tag_samp/strata_tot >=0.6
bys strata: replace sample_5 = 1 if tag_samp/strata_tot <=1.0 & tag_samp/strata_tot >=0.8



* Just checking the sample sizes
tab sample_1 if orig_tag == 1
tab sample_2 if orig_tag == 1
tab sample_3 if orig_tag == 1
tab sample_4 if orig_tag == 1
tab sample_5 if orig_tag == 1



gen fold_group = 0
replace fold_group = 1 if sample_1 == 1
replace fold_group = 2 if sample_2 == 1
replace fold_group = 3 if sample_3 == 1
replace fold_group = 4 if sample_4 == 1
replace fold_group = 5 if sample_5 == 1

by strata: tab fold_group if orig_tag 

* Now actually fitting the model on data for 4/5 groups and then predicting for the fifth

mixed mic `final'  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:clsi_test eucast_test if sample_1 == 0, noconstant stddev mle residuals(independent, by(test_type)) covariance(unstructured)
predict x_1
mixed mic `final'  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:clsi_test eucast_test if sample_2 == 0, noconstant stddev mle residuals(independent, by(test_type)) covariance(unstructured)
predict x_2
mixed mic `final'  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:clsi_test eucast_test if sample_3 == 0, noconstant stddev mle residuals(independent, by(test_type)) covariance(unstructured)
predict x_3
mixed mic `final'  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:clsi_test eucast_test if sample_4 == 0, noconstant stddev mle residuals(independent, by(test_type)) covariance(unstructured)
predict x_4
mixed mic `final'  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:clsi_test eucast_test if sample_5 == 0, noconstant stddev mle residuals(independent, by(test_type)) covariance(unstructured)
predict x_5




* Getting the EUCAST only predictions out from this, and then transforming to predicted MIC

* First recoding into whole number values, note here we always!!! round up as MIC of intermediate means the first whole number MIC sufficient to inhibit growth will be above

recode x_1 (10/11=11) (9/10=10) (8/9=9) (7/8=8) (6/7=7) (5/6=6) (4/5=5) (3/4=4) (2/3=3) (1/2=2) (0/1=1)
recode x_2 (10/11=11) (9/10=10) (8/9=9) (7/8=8) (6/7=7) (5/6=6) (4/5=5) (3/4=4) (2/3=3) (1/2=2) (0/1=1)
recode x_3 (10/11=11) (9/10=10) (8/9=9) (7/8=8) (6/7=7) (5/6=6) (4/5=5) (3/4=4) (2/3=3) (1/2=2) (0/1=1)
recode x_4 (10/11=11) (9/10=10) (8/9=9) (7/8=8) (6/7=7) (5/6=6) (4/5=5) (3/4=4) (2/3=3) (1/2=2) (0/1=1)
recode x_5 (10/11=11) (9/10=10) (8/9=9) (7/8=8) (6/7=7) (5/6=6) (4/5=5) (3/4=4) (2/3=3) (1/2=2) (0/1=1)

* next getting the EUCAST bit  i.e. genetating test specific predictions *

gen eucast_1 = 0
replace eucast_1 = x_1 if test_type == 1
sort guuid eucast_1
bys guuid: replace eucast_1  = eucast_1[_N]

gen eucast_2 = 0
replace eucast_2 = x_2 if test_type == 1
sort guuid eucast_2
bys guuid: replace eucast_2  = eucast_2[_N]

gen eucast_3 = 0
replace eucast_3 = x_3 if test_type == 1
sort guuid eucast_3
bys guuid: replace eucast_3  = eucast_3[_N]

gen eucast_4 = 0
replace eucast_4 = x_4 if test_type == 1
sort guuid eucast_4
bys guuid: replace eucast_4 = eucast_4[_N]

gen eucast_5 = 0
replace eucast_5 = x_5 if test_type == 1
sort guuid eucast_5
bys guuid: replace eucast_5  = eucast_5[_N]


* next getting the clsi bit *

gen clsi_1 = 0
replace clsi_1 = x_1 if test_type == 0
sort guuid clsi_1
bys guuid: replace clsi_1  = clsi_1[_N]

gen clsi_2 = 0
replace clsi_2 = x_2 if test_type == 0
sort guuid clsi_2
bys guuid: replace clsi_2  = clsi_2[_N]

gen clsi_3 = 0
replace clsi_3 = x_3 if test_type == 0
sort guuid clsi_3
bys guuid: replace clsi_3  = clsi_3[_N]

gen clsi_4 = 0
replace clsi_4 = x_4 if test_type == 0
sort guuid clsi_4
bys guuid: replace clsi_4 = clsi_4[_N]

gen clsi_5 = 0
replace clsi_5 = x_5 if test_type == 0
sort guuid clsi_5
bys guuid: replace clsi_5  = clsi_5[_N]

* Then censoring similar to MICs as the external prediction ,  note however there is an issue here
* With fewer sensitive MICs the wild type distribution in particular is less well characterised, the model sometimes can only predict a minimum EUCAST MIC  8
* Two possible approaches, use rounded MICs to maintain the same range 
* (this does not agree with the biological approach we took for external prediction where we felt MICs between two values would only be observed at the value above*) 
* or round up and truncate the range to be consistent across the whole set (i.e. predicting <= 8 and then higher MICs
*  Note we also took the decision to censor data at the top end at >32/2 to match the external prediction. This criteria can be lifted and makes little difference to overall predictions



* Predictions can only go to <= 4/4 =2
* observed can only go to >32/2 = 6

gen observed_eucastmic = ln(eucast_mic)/ln(2)
replace observed_eucastmic = 3 if observed_eucastmic <= 3
*managing upper censoring
replace observed_eucastmic = 6 if observed_eucastmic >= 6
replace eucast_1 = 6 if eucast_1 >= 6
replace eucast_2 = 6 if eucast_2 >= 6
replace eucast_3 = 6 if eucast_3 >= 6
replace eucast_4 = 6 if eucast_4 >= 6
replace eucast_5 = 6 if eucast_5 >= 6


* Managing the lower censoring issue
replace eucast_1 = 3 if eucast_1 <= 3
replace eucast_2 = 3 if eucast_2 <= 3
replace eucast_3 = 3 if eucast_3 <= 3
replace eucast_4 = 3 if eucast_4 <= 3
replace eucast_5 = 3 if eucast_5 <= 3

gen observed_clsimic = ln(clsi_mic)/ln(2)
replace observed_clsimic = 2 if observed_clsimic <= 2
* managing upper censoring
replace observed_clsimic = 6 if observed_clsimic >= 6

replace clsi_1 = 6 if clsi_1 >= 6
replace clsi_2 = 6 if clsi_2 >= 6
replace clsi_3 = 6 if clsi_3 >= 6
replace clsi_4 = 6 if clsi_4 >= 6
replace clsi_5 = 6 if clsi_5 >= 6

* No lower censoring issue for CLSI



* Then we will measure performance using MIC agreement, essential agreement and categorical agreement.
gen mic_eucastxval1 = abs(eucast_1 - observed_eucastmic) if sample_1 == 1
gen mic_eucastxval2 = abs(eucast_2 - observed_eucastmic) if sample_2 == 1
gen mic_eucastxval3 = abs(eucast_3 - observed_eucastmic) if sample_3 == 1
gen mic_eucastxval4 = abs(eucast_4 - observed_eucastmic) if sample_4 == 1
gen mic_eucastxval5 = abs(eucast_5 - observed_eucastmic) if sample_5 == 1


gen mic_clsixval1 = abs(clsi_1 - observed_clsimic) if sample_1 == 1
gen mic_clsixval2 = abs(clsi_2 - observed_clsimic) if sample_2 == 1
gen mic_clsixval3 = abs(clsi_3 - observed_clsimic) if sample_3 == 1
gen mic_clsixval4 = abs(clsi_4 - observed_clsimic) if sample_4 == 1
gen mic_clsixval5 = abs(clsi_5 - observed_clsimic) if sample_5 == 1

keep if orig_tag == 1

* So now each of these variables give us the prediction performance for the left out fold in each case, these are then summed for figures used in supplementary.

tab mic_eucastxval1 
tab mic_eucastxval2
tab mic_eucastxval3 
tab mic_eucastxval4 
tab mic_eucastxval5 

tab mic_clsixval1 
tab mic_clsixval2
tab mic_clsixval3 
tab mic_clsixval4 
tab mic_clsixval5 



gen observed_eucastmic = ln(eucast_mic)/ln(2)
replace observed_eucastmic = 3 if observed_eucastmic <= 3


* Finally looking at agreement with BD Phoenix

gen bd_crossval = bd_phoenix_transformed
replace bd_crossval = 3 if bd_crossval <= 3

gen bd_xval1 = abs(eucast_1 - bd_crossval) if sample_1 == 1
gen bd_xval2  = abs(eucast_2 - bd_crossval) if sample_2 == 1
gen bd_xval3  = abs(eucast_3 - bd_crossval) if sample_3 == 1
gen bd_xval4 = abs(eucast_4 - bd_crossval) if sample_4 == 1
gen bd_xval5  = abs(eucast_5 - bd_crossval) if sample_5 == 1


tab bd_xval1
tab bd_xval2
tab bd_xval3
tab bd_xval4
tab bd_xval5
