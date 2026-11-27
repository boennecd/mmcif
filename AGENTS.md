# Package documentation

- Edit `vignettes/mmcif.Rmd` to change the README or the main vignette.
- From the package root, run `Rscript vignettes/create-vignette.R` to update
  `README.md`, `vignettes/mmcif.html`, and the README figures.
- To run that script, install mmcif, bench, fastGHQuad, mvtnorm, numDeriv, and
  rmarkdown. Use Pandoc 2.10.1 or later for GitHub math.
- Put each matrix row on a separate source line so GitHub preserves the row
  breaks.
