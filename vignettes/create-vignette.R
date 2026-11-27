(function(){
  if (basename(path = getwd()) != "vignettes")
    setwd(dir = "vignettes")

  rmarkdown::render(input = "mmcif.Rmd", output_format = "html_document")
  rmarkdown::render(
    input = "mmcif.Rmd",
    output_format = rmarkdown::github_document(
      html_preview = FALSE, math_method = "default"),
    output_file = "README.md")
  file.copy(from = "README.md", to = "../README.md", overwrite = TRUE)
  unlink(x = "README.md")

  figures <- list.files(path = "man/figures", full.names = TRUE)
  file.copy(from = figures, to = "../man/figures", overwrite = TRUE)
  unlink(x = "man", recursive = TRUE)
})()
