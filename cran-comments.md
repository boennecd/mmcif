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

The Clang sanitizer checks used Ubuntu 26.04.1 LTS (x86_64), R-devel
(2026-11-25 r90991), and Clang and Flang 23.1.6. Both used
`--no-manual --extra-arch --no-stop-on-test-error`.

- AddressSanitizer and UndefinedBehaviorSanitizer ran together. This run
  disabled leak detection, allocation and deallocation mismatch checks,
  alignment checks, and floating-point division-by-zero checks.
- The separate UndefinedBehaviorSanitizer run disabled function type checks.

The attempted Fedora checks used Fedora 44 (x86_64), R-devel
(2026-11-25 r90991), GCC and GNU Fortran 16.2.1.

- The GCC sanitizer job used AddressSanitizer, UndefinedBehaviorSanitizer,
  and strict bounds checks. Leak detection was disabled.
- The Valgrind 3.27.1 job used
  `--no-manual --use-valgrind --extra-arch --no-stop-on-test-error`.
  Memcheck used `--track-origins=yes --leak-check=full`.

## R CMD check results

The local checks reported no errors or warnings. They reported the expected
NOTE for a new submission from the archive. The most recent check also could
not verify the clock against its time services. All files in the source
archive had timestamps at or before the local clock.

The hosted checks reported no errors, warnings, or NOTEs. An earlier R-devel
check flagged the deprecated `.Dim` attribute in test fixtures. The fixtures
now use `dim`.

Both Clang checks reported no errors, warnings, NOTEs, or sanitizer
diagnostics. They ran the package examples and tests.

The GCC sanitizer job stopped during dependency installation. A prebuilt
dependency required `libRlapack.so`, which the container's R build does
not provide.

The Valgrind job reached the runner's time limit in the parallel likelihood
gradient example after OpenBLAS reported nested threading. The package tests
did not run, and Memcheck did not report a summary.

The GCC sanitizer and Valgrind checks still need a complete run.
