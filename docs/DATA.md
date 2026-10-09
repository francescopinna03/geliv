# Data and comparability

No dataset is shipped with the package. The downloader is explicit. For the
Angrist and Krueger extract it compares the size and the SHA-256 hash of the
archive with the tracked record `data/asciiqob.zip.source.json` and stops on any
difference. For every source it writes the URL, the size, the hash and the date
of the check to an untracked local file ending in `.check.json`, so that a
download never alters a tracked file.

## Card

The Wooldridge dataset is distributed by Boston College at
https://fmwww.bc.edu/ec-p/data/wooldridge/card.dta and can also be loaded in
Stata with `bcuse card, clear`. The benchmark controls are
`exper expersq black smsa south smsa66 reg662-reg669` together with the
intercept. The expected sample has 3010 observations and the 2SLS coefficient
with `nearc4` is about 0.131504. A numerical benchmark validates no economic
identification argument by itself.

## Angrist and Krueger, minimal extract

The original archive of Joshua Angrist at MIT provides the extract at
https://economics.mit.edu/sites/default/files/publications/asciiqob.zip

Its five columns are the log weekly wage, years of education, year of birth,
quarter of birth and state of birth. The layout can be checked against the
original program https://economics.mit.edu/sites/default/files/inline-files/ak91.sas
and the construction of the 180 instruments is documented in
https://economics.mit.edu/sites/default/files/publications/Jive180.sas

The software keeps men born between 1930 and 1939 and the FIPS codes of the 50
states and the District of Columbia. The effective sample must be recorded by
the run itself and never imposed to match a target number of observations. The
instrument matrices must have full rank, and when some groups are missing in a
subsample the software must stop calling a rank-deficient matrix "180
instruments".

The interactions between quarter and year use quarters one to three as the
included categories, whereas Jive180 uses quarters two to four. With year fixed
effects in the model the two sets span the same instrument space, so the choice
of base leaves the model unchanged. Michigan is the omitted state, as in the
original file. In the controlled design the comparisons with 3 and 30
instruments also include state effects, which in the original file enter only
the specification with 180 instruments.

`20_ak_prepare.do` also builds `gel_cell`, the group identifier of year, quarter
and state of birth, and stores its name in the global `GEL_CELL`. Every excluded
instrument and every control is a function of these three variables, so all of
them are constant within the cells, and the estimation scripts pass the
identifier to `geliv` through `cells()`.

## Full extract

`python3 python/download_data.py ak-full` downloads
https://economics.mit.edu/sites/default/files/inline-files/NEW7080_1.rar

The archive must be opened with a RAR utility to obtain the original Stata
dataset. It contains several censuses and cohorts, so it cannot be used as a
whole as the 1980 census sample of men born between 1930 and 1939. For Table V
the official recodes are in
https://economics.mit.edu/sites/default/files/inline-files/QOB%20Table%20V.do

The original program maps v4 to EDUC, v9 to LWKLYWGE, v16 to CENSUS, v18 to QOB
and v27 to YOB, and contains the recodes of the controls. A comparison with
additional covariates requires an explicitly prepared common dataset, an updated
`GEL_W` and an identical sample across methods. The package makes no guess about
the state-of-birth column in the full extract.

## Memory and time

With 330 thousand observations and 240 moments, a full double-precision moment
matrix would occupy about 634 MB. The row-level engine never forms it, and a
block of 4096 rows by 240 columns occupies about 7.9 MB. Stata variables, Mata
views or copies of X and Z and the weight vectors remain in memory. The inner
quadratic products cost work of the order of n m² per Newton iteration, which
blocking bounds in memory and leaves unchanged in arithmetic. On Stata/BE 19.0
a single row-level CUE estimate with 3 excluded instruments, that is 63 moments,
took 362 seconds on the full Angrist and Krueger sample.

The cell-level engine selected by `cells()` stores the instruments and the
controls once per cell, so memory is dominated by the outcome, the endogenous
regressors and a few vectors of length n. Its cost per evaluation is of the
order of n + C m², with C the number of cells, about 2040 in this design. With 3
excluded instruments the same CUE estimate took 52 seconds and ET took 80
seconds. At 63 moments the cell products are already negligible and the
individual passes over the 329509 observations dominate, so the gain relative to
the row-level engine grows with the number of moments.
