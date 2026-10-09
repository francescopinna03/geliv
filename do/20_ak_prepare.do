version 16.0
* Minimal official AK extract contains exactly these five fields.
* It is sufficient for a controlled 3/30/180 instrument comparison.
clear
infile double lwklywge educ yob qob pob using data/asciiqob.txt
assert inlist(qob,1,2,3,4)
replace yob=yob+1900 if yob<100
keep if inrange(yob,1930,1939)
drop if missing(lwklywge,educ,yob,qob,pob)
* Valid US-state/DC FIPS only, common sample in every specification.
local states 1 2 4 5 6 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56
generate byte us=0
foreach s of local states {
    replace us=1 if pob==`s'
}
keep if us
drop us
local W
forvalues y=1931/1939 {
    generate byte yr`y'=yob==`y'
    local W `W' yr`y'
}
* Michigan (26) omitted, following Jive180.sas.
foreach s of local states {
    if `s'!=26 {
        generate byte st`s'=pob==`s'
        local W `W' st`s'
    }
}
local Z3
local Z30
local Z150
forvalues q=1/3 {
    generate byte q`q'=qob==`q'
    local Z3 `Z3' q`q'
    forvalues y=1930/1939 {
        generate byte q`q'y`y'=(qob==`q')*(yob==`y')
        local Z30 `Z30' q`q'y`y'
    }
    foreach s of local states {
        if `s'!=26 {
            generate byte q`q's`s'=(qob==`q')*(pob==`s')
            local Z150 `Z150' q`q's`s'
        }
    }
}
local Z180 `Z30' `Z150'
local n3: word count `Z3'
local n30: word count `Z30'
local n180: word count `Z180'
assert `n3'==3 & `n30'==30 & `n180'==180
* Global lists are explicit; no fabricated 180 instruments or random dummies.
global GEL_Y lwklywge
global GEL_X educ
global GEL_W `W'
global GEL_SPECS 3 30 180
global GEL_Z3 `Z3'
global GEL_Z30 `Z30'
global GEL_Z180 `Z180'
global GEL_ZMAX `Z180'
egen int gel_cell=group(yob qob pob)
quietly summarize gel_cell
assert r(max)<=2040
global GEL_CELL gel_cell
display "Instruments and controls are constant within " r(max) " birth year, quarter and state cells."
compress
save output/ak_ready.dta, replace
display "Common sample N=" _N "; 59 controls plus intercept; 3/30/180 excluded instruments."
display "State controls held fixed: this is a controlled extension, not a literal Table V replication."
