
library(rmarkdown)

# get rmd files present in coverage_dir
coverage_dir <- "/home/pmonteagudo/workspace/silencing_project/results/summary-coverage"
rmd_files <- list.files(coverage_dir, pattern = "Rmd$", full.names = TRUE)

# render rmd files into pdf_coverage_dir
pdf_coverage_dir <- "/home/pmonteagudo/workspace/silencing_project/results/summary-coverage-pdf"
if (!file.exists(pdf_coverage_dir)){
  dir.create(pdf_coverage_dir, recursive=TRUE, showWarnings = TRUE)
}
for (ff in rmd_files) {
  rmarkdown::render(ff, output_format = "all", output_dir = pdf_coverage_dir, clean = TRUE)
}

# clean tex files
rmd_files <- list.files(coverage_dir, pattern = "tex$", full.names = TRUE)
do.call(file.remove, list(rmd_files))
