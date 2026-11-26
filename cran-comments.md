## New version submission

CRAN archived mmcif on 2025-11-07 because check issues remained after
reminders. Version 0.1.2 declares R 4.1.0 as the minimum for the syntax used
in the examples and adds the missing package anchor for the link to
`alabama::auglag()`.

## Test environments

- Debian 13 (x86_64), R 4.5.3, GCC and GNU Fortran 15.3.0.
  Built the source tarball and ran `R CMD check --as-cran`.

Hosted checks used `--no-manual --as-cran` with CRAN incoming checks disabled.

- Ubuntu 26.04.1 LTS (x86_64), GCC and GNU Fortran 15.2.0.
  Checked R 4.5.3, R 4.6.2, and R-devel (2026-11-24 r90980).
- macOS Tahoe 26.7.2 (arm64), R 4.6.2, Apple Clang 21.0.0, and GNU Fortran
  14.2.0.
- Windows Server (x86_64, build 26100), R 4.6.2 (UCRT), GCC and GNU Fortran
  14.3.0.

## R CMD check results

The local check reported no errors or warnings. Its only NOTE identifies
the package as a new submission from the archive.

The hosted checks reported no errors, warnings, or NOTEs. An earlier R-devel
check flagged the deprecated `.Dim` attribute in test fixtures. The fixtures
now use `dim`.
