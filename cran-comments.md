## New version submission

CRAN archived mmcif on 2025-11-07 because check issues remained after
reminders. Version 0.1.2 declares R 4.1.0 as the minimum for the syntax used
in the examples and adds the missing package anchor for the link to
`alabama::auglag()`.

## Test environment

- Debian 13 (x86_64), R 4.5.3, GCC and GNU Fortran 15.3.0.
  Built the source tarball and ran `R CMD check --as-cran`.

## R CMD check results

No errors or warnings. The only NOTE is the CRAN incoming check for a new
submission of an archived package.
