clear all
log using run_selftest.log, text replace
do do/00_setup.do
do do/01_selftest.do
log close
