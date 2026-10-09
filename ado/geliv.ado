*! GELIV 0.1 07oct2026 - research implementation, Stata 16+
program define geliv, eclass
    version 16.0
    syntax varname(numeric) [if] [in], ENDog(varlist numeric) ///
        INSTruments(varlist numeric) [CONTrols(varlist numeric) ///
        METHOD(string) BLOCK(integer 4096) FROM(name) WEIGHTS(name) NOCONStant]
    if "`method'" == "" local method et
    if !inlist("`method'", "et", "el", "cue") error 198
    local gamma = cond("`method'"=="et",0,cond("`method'"=="el",-1,1))
    marksample touse
    markout `touse' `endog' `instruments' `controls'
    tempvar cons
    local x `endog' `controls'
    local z `instruments' `controls'
    local coef `x'
    if "`noconstant'" == "" {
        quietly generate double `cons' = 1 if `touse'
        local x `x' `cons'
        local z `z' `cons'
        local coef `coef' _cons
    }
    if "`weights'" != "" {
        confirm new variable `weights'
        quietly generate double `weights' = .
    }
    tempname b V D S
    mata: gel_stata("`varlist'","`x'","`z'","`touse'",`gamma',`block', ///
        "`from'","`b'","`V'","`D'","`weights'","`S'")
    matrix colnames `D' = converged q Jgel Jdistance inner_error outer_error ///
        moment_error ess min_nw max_nw negative N m p
    if `D'[1,1] != 1 {
        matrix list `S'
        di as error "No certified numerical solution. Try rescaling or alternative starting values."
        exit 430
    }
    matrix colnames `b' = `coef'
    matrix colnames `V' = `coef'
    matrix rownames `V' = `coef'
    ereturn post `b' `V', esample(`touse') depname(`varlist')
    ereturn matrix diagnostics = `D', copy
    ereturn matrix starts = `S'
    ereturn scalar N = `D'[1,12]
    ereturn scalar moments = `D'[1,13]
    ereturn scalar Jgel = `D'[1,3]
    ereturn scalar Jdistance = `D'[1,4]
    ereturn scalar ess = `D'[1,8]
    ereturn scalar converged = 1
    ereturn local method "`method'"
    ereturn local cmd "geliv"
    di as text "GEL `method': conventional HC0 strong-identification covariance"
    ereturn display
end
