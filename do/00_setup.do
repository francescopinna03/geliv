version 16.0
set more off
* Run from the extracted toolkit directory. No packages installed silently.
adopath ++ "ado"
capture mkdir output
capture mkdir data
* One initialization per Stata session. This explicitly clears Mata objects.
mata: mata clear
do mata/geliv.mata
do mata/weakiv_geometry.mata
display "Core loaded. Run do/01_selftest.do, then do/10_card.do."
* Optional benchmark dependencies (user runs once):
* ssc install ivreg2
* ssc install ranktest
* ssc install avar
* ssc install weakiv
* ssc install bcuse
