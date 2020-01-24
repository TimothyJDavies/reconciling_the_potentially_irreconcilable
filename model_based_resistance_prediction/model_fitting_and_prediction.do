

qui cap log close
log using model_fitting.log, replace
* Note the stata output of this model is in the log model_fitting.log
* Also not all this assumes STATAs working directory is the same one as the model_base.csv file

import delimited "model_base.csv",  encoding(ISO-8859-1) clear


cap gen out_sample = 0 
replace out_sample = 1 if complete_agardil == 0
replace test_type = 1 if complete_agardil == 0

egen pickone_batch = tag(batch_no)
egen pickone_guuid = tag(guuid)
rename coa_mic bd_phoenix_mic

* Variance components

* The way I currently see this, without using any cap genetic data / epi data, you intrinsically have 
* a 2 level hierarchical or a 2 level cross classified model. note the 2 level cross classified
* or 2 level cross classified with interaction between the 2 classifications
estimates clear
* Models

mixed mic || guuid: 
est store est_1
mixed mic || batch_no: 
est store est_2
mixed mic || _all:R.batch_no || guuid: 
est store est_3
egen batchXguuid = group(batch_no guuid)
mixed mic || _all:R.batch_no || guuid: || batchXguuid:
est store est_4

* The best from these models is the est_3, so visualising est 3
est restore est_3


* Random effects vizualisation
predict u0 , reffects relevel(guuid)
predict v0 , reffects releve(_all)
predict u0se , reses relevel(guuid)
predict v0se , reses relevel(_all)
predict residuals , residuals
cap egen v0rank =  rank(v0) if pickone_batch == 1
cap egen u0rank =  rank(u0) if pickone_guuid == 1
sort v0rank
list batch_no v0 v0se if pickone_batch ==1
qnorm v0 if pickone_batch == 1 , ylab(-.4 (.2) .6) aspectratio(1)
summ u0
qnorm u0 if pickone_guuid == 1, aspectratio(1)
cap gen labheight = v0 + 1.96*v0se +0.05
serrbar v0 v0se v0rank if pickone_batch == 1, scale(1.96) yline(0) addplot(scatter labheight v0rank, msymbol(non) mlabel(batch_no)  mlabposition(1) mlabangle(vertical) mlabcolor(navy))
serrbar u0 u0se u0rank if pickone_guuid == 1 , scale(1.96) yline(0)
qnorm residuals , aspectratio(1)
drop u0 u0se residuals u0rank v0 v0se v0rank labheight

estimates clear
drop batchXguuid


/* Now adding the fixed components
Note for this we have several a_priori variables
test_type => As we can clearly see from our graphical analysis
for each of the "resistance associated effects" we have found, our model will include a representation of these
Given all thinking to this point is presence absence and binary resistance classification , we include each
of these as a binary variable to begin with
so this then means our base model

mic = test_type + blactxm_buj2be_bin + blaother_bla_bin + blaoxa_buj2b_bin + blashv_buj2b_bin + blatem_buj2b_bin + tempromter_sigtem 
+ ompnon_functioning + ampcprother_sigampcpr + ampcprg134a + RandomEffects

Note the tempromoter is only significant if there is a significant mutation and the ompnon_functioning is only significant if there is another blm

We will study interactions/modelling effects in the order of blms, promoters omps, then on commonality of evidence
i.e.
test_type
blatem
blactxm
blaoxa
blashv
blaother
temprother_sigtem 
ampcprother_sigampcpr
ompnon_functioning


note 
blaother
temprother_sigtem 
ampcprg134a 
ampcprother_sigampcpr
ompnon_functioning
all binary

* Each of the "transmissible/multiple copy" elements we will model as
1) Binary
2) 2 level binary + linear
3) 2 level binary + log

for the continuous modelling we will do the following to our copyno stat
we will recentre at 0 , so any elements with copy no between 0 and 1 we will set to 1
we will then minus 1 from the copy number and model that variable
note this var will be valled var_cont

Likewise for the log 
but this time no minus 1 step , as log(1) = 0 
it will be called var_log

Thresholds for inclusion <  0.05
Thresholds for interaction terms, it needs to be present in at least 5 samples
*/
*/ 

local original_nonbinblms "blatem_buj2b blactxm_buj2be blaoxa_buj2d blashv_buj2b"
local non_blms "blaother_bla temprother_sigtem ampcprother_sigampcpr ompnon_functioning"

* Dealikng with outliers
replace blashv_buj2b = 16.26761  if blashv_buj2b > 20

* Seeting up the BLMs
cap program drop make_bin
program make_bin
	local x `1'
	cap gen `x'_bin = 0
	replace `x'_bin = 1 if `x' != 0
end

cap program drop make_cont
program make_cont
	local x `1'
	cap gen `x'_cont = `x'
	replace `x'_cont = 1 if `x' < 1 & `x' > 0 
	replace `x'_cont = `x'_cont - 1 if `x'_cont > 0
end

cap program drop make_log
program make_log
	local x `1'
	cap gen `x'_log = `x'
	replace `x'_log = 1 if `x' < 1 & `x' > 0 
	replace `x'_log = log(`x'_log)/log(2) if `x'_log > 0
end


recode blaother_bla 1/99=1

foreach var in `original_nonbinblms' {
	make_bin `var'
	make_cont `var'
	make_log `var'
}

* baseline model test_type blatem_buj2b_bin blactxm_buj2be_bin blaoxa_buj2d_bin blashv_buj2b_bin `non_blms'

global baseline "test_type blatem_buj2b_bin blactxm_buj2be_bin blaoxa_buj2d_bin blashv_buj2b_bin `non_blms'"

* before checking how each of these effects model, I am going to check if there is any model improvement in adding superimposed mlst strain effects

mixed mic  ${baseline} || _all:R.batch_no || guuid: , stddev mle
est store est_1
mixed mic mlst_131 mlst_69 mlst_73 ${baseline} || _all:R.batch_no || guuid: , stddev mle
est store est_2
lrtest est_1 est_2
mixed mic mlst_131  ${baseline} || _all:R.batch_no || guuid: , stddev mle
est store est_3
lrtest est_1 est_3

*Even at this early stage strain effects don't appear to add much to the model (note the basline is the mlst "OTHER" cat
*Therefore strain not included
* Not looking at modelling each of the fixed effects, chosing the best fit for each one

cap program drop fixed_effect_model
program fixed_effect_model
	estimates clear
	local x `1'
	local subin_`x' "${baseline}"
	local not "`x'_bin"
	local subin_`x' : list subin_`x'- not
	mixed mic `subin_`x'' `x'_bin || _all:R.batch_no || guuid: , stddev mle
	est store bin_
	mixed mic `subin_`x'' `x'_bin `x'_cont || _all:R.batch_no || guuid: , stddev mle
	est store bin_cont_
	mixed mic `subin_`x'' `x'_bin `x'_log || _all:R.batch_no || guuid: , stddev mle
	est store bin_log_
end

cap gen blatem_buj2b_binXtest_type = blatem_buj2b_bin*test_type
cap gen  blatem_buj2b_logXtest_type =  blatem_buj2b_log*test_type
cap gen blactxm_buj2be_binXtest_type = blactxm_buj2be_bin*test_type
cap gen blaoxa_buj2d_binXtest_type = blaoxa_buj2d_bin*test_type
cap gen blashv_buj2b_logXtest_type = blashv_buj2b_log*test_type
cap gen blashv_buj2b_binXtest_type = blashv_buj2b_bin*test_type
cap gen blaother_blaXtest_type = blaother_bla*test_type

cap gen temprother_sigtemXtest_type = temprother_sigtem*test_type
cap gen ampcprother_sigampcprXtest_type = ampcprother_sigampcpr*test_type
cap gen ompnon_functioningXtest_type = ompnon_functioning*test_type
cap gen blatemXblaoxa = blatem_buj2b_bin*blaoxa_buj2d_bin
cap gen blatemXblactxm = blatem_buj2b_bin*blactxm_buj2be_bin

cap gen blatemXblashv = blatem_buj2b_bin*blashv_buj2b_bin
cap gen blatemXblaother = blatem_buj2b_bin*blaother_bla
cap gen blatemXampcprother_sigampcpr = blatem_buj2b_bin*ampcprother_sigampcpr
cap gen blactxmXblaoxa = blactxm_buj2be_bin*blaoxa_buj2d_bin
cap gen blactxmXblashv = blactxm_buj2be_bin*blashv_buj2b_bin
cap gen blactxmXblaother = blactxm_buj2be_bin*blaother_bla
cap gen blactxmXampcprother_sigampcpr = blactxm_buj2be_bin*ampcprother_sigampcpr

cap gen blaoxaXblashv = blaoxa_buj2d_bin*blashv_buj2b_bin
cap gen blaoxaXblaother = blaoxa_buj2d_bin*blaother_bla
cap gen blaoxaXampcprother_sigampcpr = blaoxa_buj2d_bin*ampcprother_sigampcpr
cap gen blashvXblaother  = blashv_buj2b_bin*blaother_bla
cap gen blashvXampcprother_sigampcpr = blashv_buj2b_bin*ampcprother_sigampcpr


* TEM
fixed_effect_model blatem_buj2b
est stat *
* model as bin_log_ given AIC, this is explained in the paper
global baseline "test_type blatem_buj2b_bin blatem_buj2b_log blactxm_buj2be_bin blaoxa_buj2d_bin blashv_buj2b_bin `non_blms'"

* CTXM
fixed_effect_model blactxm_buj2be
est stat *
* Here now binary has the lowest AIC, but note node of these modellings particularly change the
* AIC, so we stick with binary here

*OXA
fixed_effect_model blaoxa_buj2d
est stat *
* Here now binary has the lowest AIC,
* AIC, so we stick with binary here

*SHV
fixed_effect_model blashv_buj2b
est stat *

* Bin log the best on AIC although the difference is smaller her

global baseline "test_type blatem_buj2b_bin blatem_buj2b_log blactxm_buj2be_bin blaoxa_buj2d_bin blashv_buj2b_bin blashv_buj2b_log `non_blms'"

* Now we have out main effects for the fixed part of the model
* The next step is to look at interaction effects in this fixed part

* Test type interactions . Given one of our questions was to explain the difference between the two testing methodologies, we checked for test type interactions first

estimates clear


mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A
mixed mic $baseline  blatem_buj2b_binXtest_type blatem_buj2b_logXtest_type   || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B
est stat *


* Significant improvement in the models by introducing this interaction
* redoing the baseline 
global baseline "test_type blatem_buj2b_bin blatem_buj2b_log blatem_buj2b_binXtest_type blatem_buj2b_logXtest_type blactxm_buj2be_bin blaoxa_buj2d_bin blashv_buj2b_bin blashv_buj2b_log `non_blms' "


mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A
mixed mic $baseline  blactxm_buj2be_binXtest_type || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B
est stat *


* Significant improvement in the models by introducing this interaction
* redoing the baseline 
global baseline "test_type blatem_buj2b_bin blatem_buj2b_log blatem_buj2b_binXtest_type blatem_buj2b_logXtest_type blactxm_buj2be_bin blactxm_buj2be_binXtest_type blaoxa_buj2d_bin blashv_buj2b_bin blashv_buj2b_log `non_blms' "


mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A
mixed mic $baseline blaoxa_buj2d_binXtest_type || _all:R.batch_no || guuid: , stddev mle
est store B

lrtest A B
est stat *


* Significant improvement in the models by introducing this interaction
* redoing the baseline
global baseline "test_type blatem_buj2b_bin blatem_buj2b_log blatem_buj2b_binXtest_type blatem_buj2b_logXtest_type blactxm_buj2be_bin blactxm_buj2be_binXtest_type blaoxa_buj2d_bin blaoxa_buj2d_binXtest_type blashv_buj2b_bin blashv_buj2b_log `non_blms' "


mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A
mixed mic $baseline blashv_buj2b_logXtest_type blashv_buj2b_binXtest_type || _all:R.batch_no || guuid: , stddev mle
est store B

lrtest A B
est stat *


* Significant improvement in the models by introducing this interaction
* redoing the baseline
* local non_blms "blaother_bla temprother_sigtem ampcprg134a ampcprother_sigampcpr ompnon_functioning"

global baseline "test_type blatem_buj2b_bin blatem_buj2b_log blatem_buj2b_binXtest_type blatem_buj2b_logXtest_type blactxm_buj2be_bin blactxm_buj2be_binXtest_type blaoxa_buj2d_bin blaoxa_buj2d_binXtest_type blashv_buj2b_bin blashv_buj2b_log blashv_buj2b_logXtest_type blashv_buj2b_binXtest_type `non_blms' "


mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A
mixed mic $baseline blaother_blaXtest_type || _all:R.batch_no || guuid: , stddev mle
est store B

lrtest A B
est stat *



* Significant improvement in the models by introducing this interaction
* redoing the baseline
* local non_blms has now been encoded local _2

local _1 "test_type blatem_buj2b_bin blatem_buj2b_log blatem_buj2b_binXtest_type blatem_buj2b_logXtest_type blactxm_buj2be_bin blactxm_buj2be_binXtest_type blaoxa_buj2d_bin blaoxa_buj2d_binXtest_type blashv_buj2b_bin blashv_buj2b_log blashv_buj2b_logXtest_type blashv_buj2b_binXtest_type"
local _2 "blaother_bla blaother_blaXtest_type temprother_sigtem ampcprother_sigampcpr ompnon_functioning"
global baseline "`_1' `_2'"


mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A
mixed mic $baseline temprother_sigtemXtest_type || _all:R.batch_no || guuid: , stddev mle
est store B

lrtest A B
est stat *


* Significant improvement in the models by introducing this interaction
* redoing the baseline

local _2 "blaother_bla blaother_blaXtest_type temprother_sigtem temprother_sigtemXtest_type  ampcprother_sigampcpr ompnon_functioning"
global baseline "`_1' `_2'"


mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A
mixed mic $baseline ampcprother_sigampcprXtest_type || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B
est stat *


* Significant improvement in the models by introducing this interaction
* redoing the baseline

local _2 "blaother_bla blaother_blaXtest_type temprother_sigtem temprother_sigtemXtest_type ampcprother_sigampcpr ampcprother_sigampcprXtest_type ompnon_functioning"
global baseline "`_1' `_2'"


mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A
mixed mic $baseline ompnon_functioningXtest_type || _all:R.batch_no || guuid: , stddev mle
est store B

est stat *

* Significant improvement in the models by introducing this interaction
* Given this is true for almost all elements tested, 
* redoing the baseline

* Now I will check for first order cap gene cap gene interactions only


**** BLATEM ******
local _2 "blaother_bla blaother_blaXtest_type temprother_sigtem temprother_sigtemXtest_type ampcprother_sigampcpr ampcprother_sigampcprXtest_type ompnon_functioning ompnon_functioningXtest_type"
global baseline "`_1' `_2'"

tab blatemXblactxm   if pickone_guuid == 1
mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A
mixed mic $baseline blatemXblactxm || _all:R.batch_no || guuid: , stddev mle
est store B

lrtest A B
est stat *


* Significant improvement in the models by introducing this interaction
* redoing the baseline
mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store N

local _3 "blatemXblactxm"
global baseline "`_1' `_2' `_3'"

tab blatemXblaoxa   if pickone_guuid == 1
mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A
mixed mic $baseline blatemXblaoxa|| _all:R.batch_no || guuid: , stddev mle
est store B
local _3 "blatemXblaoxa"
global baseline "`_1' `_2' `_3'"
mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store C


* So models in the above are
* A BlaCTXM blaTEM
* B BLACTXMblaTEM and blaOXA
* C blaOXA blaTEM
* N is no interaction terms
lrtest A N
lrtest C N
lrtest B A
lrtest B C
* Including each of blaOXA and CTXM appears interactions initially important,
*  however, when then removing each there is no significant difference (p = 0.11 on LR test) the models , 
* and hence as bla OXA the stronger effect only this one is kept.

est stat *

* OXA only term appears to be best here
* Significant improvement in the models by introducing this interaction
* redoing the baseline
local _3 "blatemXblaoxa"
global baseline "`_1' `_2' `_3'"

tab blatemXampcprother_sigampcpr  if pickone_guuid == 1
* =3 so omitted


**** BLACTXM ******
local _3 "blatemXblaoxa"
global baseline "`_1' `_2' `_3'"

tab blactxmXblaoxa if pickone_guuid == 1
mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A
mixed mic $baseline blactxmXblaoxa || _all:R.batch_no || guuid: , stddev mle
est store B

est stat *
* Significant improvement in the models by introducing this interaction
* redoing the baseline
local _3 "blatemXblaoxa blactxmXblaoxa"
global baseline "`_1' `_2' `_3'"


* blactxmXblashv = 0 so omitted


tab blactxmXblaother if pickone_guuid == 1
* = 1 so omitted



tab blactxmXampcprother_sigampcpr  if pickone_guuid == 1
* = 0 so omitted

**** BLAOXA ***

tab blaoxaXblashv  if pickone_guuid == 1
* =0 so omitted


tab blaoxaXblaother if pickone_guuid == 1
* = 0



tab  blaoxaXampcprother_sigampcpr if pickone_guuid == 1
* = 0 so omitted

*** BLASHV ****


tab blashvXblaother if pickone_guuid == 1
* = 1 so omitted


tab blashvXampcprother_sigampcpr if pickone_guuid == 1
* = 0 so omitted
*/
**** SO THAT IS ALL THE INTERACTIONS. NOTE I HAVEN'T Tested the TEM promoter / omp non functionality given how these variables work, + also havent tested more than
*** Binary effects

**** So Now testing if any of it can be removed (in the same order as above) 



***ORIGINAL ****
local _1 "test_type blatem_buj2b_bin blatem_buj2b_log blatem_buj2b_binXtest_type blatem_buj2b_logXtest_type blactxm_buj2be_bin blactxm_buj2be_binXtest_type blaoxa_buj2d_bin blaoxa_buj2d_binXtest_type blashv_buj2b_bin blashv_buj2b_log blashv_buj2b_logXtest_type blashv_buj2b_binXtest_type"
local _2 "blaother_bla blaother_blaXtest_type temprother_sigtem temprother_sigtemXtest_type ampcprother_sigampcpr ampcprother_sigampcprXtest_type ompnon_functioning ompnon_functioningXtest_type"
local _3  "blatemXblaoxa blactxmXblaoxa"

local a_1 `_1'
local a_2 `_2'
local a_3 `_3'
global baseline "`a_1' `a_2' `a_3'"

mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store A

local b_1 `_1'
local b_2 `_2'
local b_3 " blactxmXblaoxa"
global baseline "`b_1' `b_2' `b_3'"
mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store B

lrtest A B
* Now testing each one to see its importance


* TEM TEST TYPE = > Stays
local b_1 "test_type blatem_buj2b_bin  blatem_buj2b_log  blactxm_buj2be_bin blactxm_buj2be_binXtest_type blaoxa_buj2d_bin blaoxa_buj2d_binXtest_type blashv_buj2b_bin blashv_buj2b_binXtest_type blashv_buj2b_log blashv_buj2b_logXtest_type "
local b_2 `_2'
local b_3 `_3'
global new "`b_1' `b_2' `b_3'"
mixed mic $new  || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B

* CTXM Test Type => Doesn't appear to add significantly to the model
local b_1 "test_type blatem_buj2b_bin blatem_buj2b_binXtest_type blatem_buj2b_log  blatem_buj2b_logXtest_type blactxm_buj2be_bin blaoxa_buj2d_bin blaoxa_buj2d_binXtest_type blashv_buj2b_bin blashv_buj2b_binXtest_type blashv_buj2b_log blashv_buj2b_logXtest_type "
local b_2 `_2'
local b_3 `_3'
global new "`b_1' `b_2' `b_3'"
mixed mic $new  || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B

*OXA test type => Stays
local b_1 "test_type blatem_buj2b_bin blatem_buj2b_binXtest_type blatem_buj2b_log  blatem_buj2b_logXtest_type blactxm_buj2be_bin blactxm_buj2be_binXtest_type blaoxa_buj2d_bin blashv_buj2b_bin blashv_buj2b_binXtest_type blashv_buj2b_log blashv_buj2b_logXtest_type "
local b_2 `_2'
local b_3 `_3'
global new "`b_1' `b_2' `b_3'"
mixed mic $new  || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B

*SHV test_type => Stays
local b_1 "test_type blatem_buj2b_bin blatem_buj2b_binXtest_type blatem_buj2b_log  blatem_buj2b_logXtest_type blactxm_buj2be_bin blactxm_buj2be_binXtest_type blaoxa_buj2d_bin blaoxa_buj2d_binXtest_type blashv_buj2b_bin  blashv_buj2b_log "
local b_2 `_2'
local b_3 `_3'
global new "`b_1' `b_2' `b_3'"
mixed mic $new  || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B

* bla_other test_type => Stays
local b_1 `_1'
local b_2 "blaother_bla  temprother_sigtem temprother_sigtemXtest_type ampcprother_sigampcpr ampcprother_sigampcprXtest_type ompnon_functioning ompnon_functioningXtest_type"
local b_3 `_3'
global new "`b_1' `b_2' `b_3'"
mixed mic $new  || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B

* TEMPROM => stays
local b_1 `_1'
local b_2 "blaother_bla blaother_blaXtest_type  temprother_sigtem ampcprother_sigampcpr ampcprother_sigampcprXtest_type ompnon_functioning ompnon_functioningXtest_type"
local b_3 `_3'
global new "`b_1' `b_2' `b_3'"
mixed mic $new  || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B

* ampC other => stays
local b_1 `_1'
local b_2 "blaother_bla blaother_blaXtest_type  temprother_sigtem temprother_sigtemXtest_type  ampcprother_sigampcpr  ompnon_functioning ompnon_functioningXtest_type"
local b_3 `_3'
global new "`b_1' `b_2' `b_3'"
mixed mic $new  || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B

* omptest_type => Stays but weaker evidence from LRT
local b_1 `_1'
local b_2 "blaother_bla blaother_blaXtest_type  temprother_sigtem temprother_sigtemXtest_type  ampcprother_sigampcpr   ampcprother_sigampcprXtest_type ompnon_functioning "
local b_3 `_3'
global new "`b_1' `b_2' `b_3'"
mixed mic $new  || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B

* tem oxa ==> Keeping the interaction improves model
local b_1 `_1'
local b_2 `_2'
local b_3 " blactxmXblaoxa"
global baseline "`b_1' `b_2' `b_3'"
mixed mic $baseline  || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B

* CTXM_oxa oxa => Keeping the interaction improves model
local b_1 `_1'
local b_2 `_2'
local b_3 "blatemXblaoxa"
global new "`b_1' `b_2' `b_3'"
global new "`b_1' `b_2' `b_3'"
mixed mic $new  || _all:R.batch_no || guuid: , stddev mle
est store B
lrtest A B

* Everything apart from CTXM test type interaction improves the model. 


**** SO Ones with little evidence at top level => CTXM test type interaction , NOTE HERE WE HAVE NO SIGNIFICANTLY HIGH COPY NUMBER SAMPLES
***** Hence little resistance appears to be caused by this. so maybe more in keeping with Low copy number unable to cause resistance
*** ALSO OF NOTE THERE IS SOME SATURATION EFFECT , gene interactions tend to be negative ****

**** SO now given the importance of the test type effects , and to see if we can include them throughout the model
*** Note HOWEVER GIVEN  the importance of test type effects, we want to keep it provided it doesn't make any other estimates unstable
**** I Will check what happens to point estimates if you add in the extra test type separated terms



* THIS SECTION IS ALL ABOUT CHECKING STABILITY OF POINT ESTIMATES

local _1 "test_type blatem_buj2b_bin blatem_buj2b_log blatem_buj2b_binXtest_type blatem_buj2b_logXtest_type blactxm_buj2be_bin blactxm_buj2be_binXtest_type blaoxa_buj2d_bin blaoxa_buj2d_binXtest_type blashv_buj2b_bin blashv_buj2b_log blashv_buj2b_logXtest_type blashv_buj2b_binXtest_type"
local _2 "blaother_bla blaother_blaXtest_type temprother_sigtem temprother_sigtemXtest_type ampcprother_sigampcpr ampcprother_sigampcprXtest_type ompnon_functioning ompnon_functioningXtest_type"
local _3  "blatemXblaoxa blactxmXblaoxa"

local a_1 `_1'
local a_2 `_2'
local a_3 `_3'
di "`a_1"'
global A  "`a_1' `a_2' `a_3'"
mixed mic $A  || _all:R.batch_no || guuid: , stddev mle
est store A


local b_1 "test_type blatem_buj2b_bin blatem_buj2b_binXtest_type blatem_buj2b_log  blatem_buj2b_logXtest_type blactxm_buj2be_bin  blaoxa_buj2d_bin blaoxa_buj2d_binXtest_type blashv_buj2b_bin blashv_buj2b_binXtest_type blashv_buj2b_log blashv_buj2b_logXtest_type "
local b_2 `_2'
local b_3 `_3'
global B "`b_1' `b_2' `b_3'"
mixed mic $B  || _all:R.batch_no || guuid: , stddev mle
est store B
* Removing term

lrtest A B
* B is margically the best model , (AIC 3698 , 3699, 3700).
* However given how close this, and the greater interpretability of separating each by test type we will keep this model for interpretatbility
* Note given the small numbers I will only look for interactions between random effect and test type. (this is the only possible given the large numbers of interaction terms each method could produce)


local a_1 `_1'
local a_2 `_2'
local a_3 `_3'
global $B "`a_1' `a_2' `a_3'"

**** Model B *****
estimates clear
mixed mic $B  || _all:R.batch_no || guuid: , stddev mle
est store est_1 
mixed mic $B  || _all:R.batch_no || guuid: , stddev mle residuals(independent, by(test_type))
est store est_2
mixed mic $B  || _all:R.batch_no || guuid:test_type , stddev mle
est store est_3
mixed mic $B  || _all:R.batch_no || guuid:test_type , stddev mle covariance(unstructured)
est store est_4
mixed mic $B || _all:R.guuid || batch_no:test_type, stddev mle 
est store est_5
*mixed mic $B || _all:R.guuid || batch_no:test_type, stddev mle covariance(unstructured)
*est store est_6
** Estimate 3 is the best here

tab batch_no, gen(id1_)
foreach var of varlist id1_* {
	cap gen `var'Xtest_type = `var'*test_type
}

mixed mic $B  || _all:R.batch_no || guuid:test_type , stddev mle residuals(independent, by(test_type)) covariance(unstructured)
est store est_7
mixed mic $B  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:test_type , stddev mle covariance(unstructured)
est store est_8
mixed mic $B  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:test_type , stddev mle residuals(independent, by(test_type)) covariance(unstructured)
est store est_9

est stat *

***** So adding in each of the other levels appers to improve fit. Note this model is now really really complicated in terms of variance components
**** Checking the effect of the extra mechanisms


estimates clear
mixed mic $A  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:test_type ,  stddev mle residuals(independent, by(test_type)) covariance(unstructured)
est store est_1
mixed mic $B  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:test_type , stddev mle residuals(independent, by(test_type)) covariance(unstructured)
est store est_2

est stat *
***** Estimates appear stable to each of these 3 fixed components. THerefore to aid interpretability estimate C Chosen.

* REPARAMETERISING FOR THE GRAPH
* INITIAL
local _1 "test_type blatem_buj2b_bin blatem_buj2b_log blatem_buj2b_binXtest_type blatem_buj2b_logXtest_type blactxm_buj2be_bin blactxm_buj2be_binXtest_type blaoxa_buj2d_bin blaoxa_buj2d_binXtest_type blashv_buj2b_bin blashv_buj2b_log blashv_buj2b_logXtest_type blashv_buj2b_binXtest_type"
local _2 "blaother_bla blaother_blaXtest_type temprother_sigtem temprother_sigtemXtest_type ampcprother_sigampcpr ampcprother_sigampcprXtest_type ompnon_functioning ompnon_functioningXtest_type"
local _3  "blatemXblaoxa blactxmXblaoxa"

local a_1 `_1'
local a_2 `_2'
local a_3 `_3'
di "`a_1"'
global A  "`a_1' `a_2' `a_3'"
mixed mic $A  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:test_type , mle residuals(independent, by(test_type)) covariance(unstructured)
*NEW
cap gen clsi_test = 0
replace clsi_test=1 if test_type == 0
cap gen eucast_test = test_type	


local elements "blatem_buj2b_bin blatem_buj2b_log blactxm_buj2be_bin blaoxa_buj2d_bin blashv_buj2b_bin blashv_buj2b_log blaother_bla temprother_sigtem ampcprother_sigampcpr ompnon_functioning"
local clsi ""
local eucast  ""

foreach var in `elements' {
	cap gen c_`var' = `var'*(1-test_type)
	cap gen e_`var' = `var'*test_type
	local clsi "`clsi' c_`var'"
	local eucast "`eucast'  e_`var'"
}

local final "test_type"
local m  : word count `clsi'
forvalues i = 1/`m' {
local k1 : word `i' of `clsi'
local k2 : word `i' of `eucast'
local final "`final' `k1' `k2'"
}

local final "`final' blatemXblactxm blactxmXblaoxa"

di "`final'"


mixed mic `final'  || _all:R.batch_no || _all:id1_*Xtest_type, noconstant cov(identity)|| guuid:clsi_test eucast_test , noconstant stddev mle residuals(independent, by(test_type)) covariance(unstructured)

* Generating predictions
predict orig_fitted /* this is an untransformed original linear prediction */
*We now predict an actual MIC. To to this we take linear predictions and then round upwards (as MICs within ranges e.g. 2.4 would only result in an observed MIC of 3)
predict pred_mic 

replace pred_mic = 11 if pred_mic <=11 & pred_mic > 10
replace pred_mic = 10 if pred_mic <=10 & pred_mic > 9
replace pred_mic = 9 if pred_mic <=9 & pred_mic > 8
replace pred_mic = 8 if pred_mic <=8 & pred_mic > 7
replace pred_mic = 7 if pred_mic <=7 & pred_mic > 6
replace pred_mic = 6 if pred_mic <=6 & pred_mic > 5
replace pred_mic = 5 if pred_mic <=5 & pred_mic > 4
replace pred_mic = 4 if pred_mic <=4 & pred_mic > 3
replace pred_mic = 3 if pred_mic <=3 & pred_mic > 2
replace pred_mic = 2 if pred_mic <=2 & pred_mic > 1
replace pred_mic = 1 if pred_mic <=1 & pred_mic > 0

* First then putting the BD Phoenix MICs into a computer readable format

gen bd_phoenix_transformed = 0

replace bd_phoenix_transformed = 1 if bd_phoenix_mic == "<=2/2"
replace bd_phoenix_transformed = 2 if bd_phoenix_mic == "4/2"
replace bd_phoenix_transformed = 3 if bd_phoenix_mic == "8/2"
replace bd_phoenix_transformed = 4 if bd_phoenix_mic == "16/2"
replace bd_phoenix_transformed = 5 if bd_phoenix_mic == "32/2"
replace bd_phoenix_transformed = 6 if bd_phoenix_mic == ">32/2"

* Next they predict on the different ranges, so censoring so they predict on the same range,
* The model cannot predict < 4/2 (estimated mean of the "wild type population"
* BD Phoenix only goes up to > 32/2
replace pred_mic = 6 if pred_mic >= 6
replace bd_phoenix_transformed = 2 if bd_phoenix_transformed <=2

*****************************************
* COMPARING PREDICTIONS ON THE TESTING SET
*****************************************

gen pred_mic_diff = abs(pred_mic - bd_phoenix_transformed)
* First removing exclusions
tab bd_phoenix_transformed pred_mic if prediction_exclusion == 0 & inagar == "False"
tab pred_mic_diff if prediction_exclusion == 0 & inagar == "False"

*Then retaining 11 exclusions
tab bd_phoenix_transformed pred_mic if inagar == "False"
tab pred_mic_diff if inagar == "False"

est replay

* finally we save a copy of this for cross validation analysis which is detailed in the suppelementary methods
save "/Users/timdavies/Dropbox_Temp_Outside/Documents/projects/BLI_inhibitor_paper/modelling_analysis/dphil_redo/final_model/model_output_for_crossval.dta", replace



** Finally pulling out the model coefficients and sending them to a separate file to make figure 6

*This code basically takes statas stored info from the model and turns in into a csv
mat b = e(b)'
mat A = e(V)
mat list A

mat B = A
forval i = 1/`= rowsof(A)' { 
	forval j = 1/`= colsof(A)' { 
		mat B[`i', `j'] = sqrt(A[`i', `j']) 
	}
} 


mat se_mod = vecdiag(B)'
matnames b
local n : word count `r(r)'
lab def coeff_labs -1 "init"
matnames b
forvalues i = 1/`n'{
local x1 : word `i' of `r(r)'
lab def coeff_labs `i' "`x1'", add
}
gen coeff_name = _n if _n <= `n'
svmat b, names(coeffs)
svmat se_mod, names(se_s)
lab values coeff_name coeff_labs

keep coeffs1 coeff_name se_s
rename coeff_name coeff_n
drop if coeffs == .
gen uci = coeffs1 + 1.96 * se_s1 
gen lci = coeffs1 - 1.96 * se_s1 
gen labheight = uci +0.05

export delimited using "/Users/timdavies/Google Drive/aac_code/model_based_resistance_prediction/model_coeffs_for_picture.csv", replace

qui cap log close


