version 16.0
clear
capture log close card
log using output/card.log, text replace name(card)
capture confirm file data/card.dta
if _rc {
    capture which bcuse
    if _rc {
        display as error "Run python3 python/download_data.py card, or install bcuse."
        exit 601
    }
    bcuse card, clear
    save data/card.dta, replace
}
else use data/card.dta, clear
local W exper expersq black smsa south smsa66 reg662 reg663 reg664 reg665 reg666 reg667 reg668 reg669
foreach v in lwage educ nearc4 nearc2 `W' {
    confirm numeric variable `v'
}
egen byte nmiss=rowmiss(lwage educ nearc4 nearc2 `W')
keep if nmiss==0
drop nmiss
assert _N==3010
regress lwage educ `W', vce(robust)
ivreg2 lwage `W' (educ=nearc4), robust
assert abs(_b[educ]-.131504)<.00001
scalar b_iv=_b[educ]
foreach method in cue et el {
    geliv lwage, endog(educ) controls(`W') instruments(nearc4) method(`method') weights(p_`method')
    assert abs(_b[educ]-scalar(b_iv))<1e-6
    assert abs(p_`method'-1/_N)<1e-7
}
* Just-identified equality is necessary, not sufficient. Add nearc2.
ivreg2 lwage `W' (educ=nearc4 nearc2), cue robust
scalar b_cue=_b[educ]
matrix b_ref=e(b)
geliv lwage, endog(educ) controls(`W') instruments(nearc4 nearc2) method(cue)
assert abs(_b[educ]-scalar(b_cue))<1e-6
foreach v in educ `W' _cons {
    assert abs(_b[`v']-b_ref[1,colnumb(b_ref,"`v'")])<1e-6
}
foreach method in et el {
    geliv lwage, endog(educ) controls(`W') instruments(nearc4 nearc2) method(`method') weights(p2_`method')
    matrix list e(diagnostics)
}
* These are design globals used by 30_mc.do and 40_geometry.do.
global GEL_Y lwage
global GEL_X educ
global GEL_W `W'
global GEL_SPECS 1 2
global GEL_Z1 nearc4
global GEL_Z2 nearc4 nearc2
global GEL_ZMAX nearc4 nearc2
save output/card_ready.dta, replace
* Independent benchmark, heteroskedastic-robust scalar inference via weakiv.
ivreg2 lwage `W' (educ=nearc4 nearc2), robust
weakiv, null(.1)
log close card
