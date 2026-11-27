## New version submission

CRAN archived mmcif on 2025-11-07 because check issues remained after
reminders. Version 0.1.2 declares R 4.1.0 as the minimum for the syntax used
in the examples and adds the missing package anchor for the link to
`alabama::auglag()`.

## Test environments

- Debian 13 (x86_64), R 4.5.3, GCC and GNU Fortran 15.3.0.
  Built the source tarball and ran `R CMD check --as-cran`.

Standard CI used `--no-manual --as-cran` with CRAN incoming checks disabled.

- Ubuntu 26.04.1 LTS (x86_64), GCC and GNU Fortran 15.2.0.
  Checked R 4.5.3, R 4.6.2, and R-devel (2026-11-24 r90980).
- macOS Tahoe 26.7.2 (arm64), R 4.6.2, Apple Clang 21.0.0, and GNU Fortran
  14.2.0.
- Windows Server (x86_64, build 26100), R 4.6.2 (UCRT), GCC and GNU Fortran
  14.3.0.

The Clang sanitizer checks used Ubuntu 26.04.1 LTS (x86_64), R-devel
(2026-11-25 r90991), and Clang and Flang 23.1.6. Both used
`--no-manual --extra-arch --no-stop-on-test-error`.

- AddressSanitizer and UndefinedBehaviorSanitizer ran together. This run
  disabled leak detection, allocation and deallocation mismatch checks,
  alignment checks, and floating-point division-by-zero checks.
- The separate UndefinedBehaviorSanitizer run disabled function type checks.

The Fedora checks used Fedora 44 (x86_64), R-devel
(2026-11-25 r90991), GCC and GNU Fortran 16.2.1.

- The GCC sanitizer job used AddressSanitizer, UndefinedBehaviorSanitizer,
  and strict bounds checks. Leak detection was disabled. The check used
  `--no-manual --no-build-vignettes --as-cran`.
- The Valgrind 3.27.1 job used one OpenBLAS thread and the R check options
  `--no-manual --use-valgrind --no-stop-on-test-error`.
  Memcheck used `--track-origins=yes --leak-check=full`.

## R CMD check results

The local check reported no errors or warnings. Its only NOTE was for
a new submission from the archive.

Standard CI reported no errors, warnings, or NOTEs. An earlier R-devel
check flagged the deprecated `.Dim` attribute in test fixtures. The fixtures
now use `dim`.

Both Clang checks reported no errors, warnings, NOTEs, or sanitizer
diagnostics. They ran the package examples and tests.

The GCC sanitizer check ran the examples and tests without sanitizer
diagnostics. It reported no errors or warnings. Its only NOTE was that
Ghostscript was unavailable for the PDF size check. The workflow now installs
Ghostscript for that job.

An earlier GCC sanitizer job stopped because a binary dependency required
`libRlapack.so`. The completed check used source builds for the affected
dependencies.

Valgrind completed the examples and tests, and `R CMD check` returned
`Status: OK`. R-hub marked the job as failed because Memcheck classified
one block of OpenMP worker thread-local storage as possibly lost in each
process (368 bytes for the examples and 384 bytes for the tests). Memcheck
reported zero definitely or indirectly lost bytes and no invalid accesses.
A standalone OpenMP program reproduced the allocation trace on Debian. The
R-hub workflow now suppresses that trace for possibly lost blocks only.
The updated GCC sanitizer and Valgrind jobs still need verification on R-hub.
