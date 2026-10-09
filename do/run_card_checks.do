version 16.0
clear all
set more off
capture log close _all
capture mkdir output
log using output/acceptance.log, text replace
do do/00_setup.do
do do/01_selftest.do
do do/10_card.do
do do/11_weakiv_benchmark.do
display as result "ALL REQUESTED CARD ACCEPTANCE GATES PASSED"
log close
