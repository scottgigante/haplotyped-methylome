rmd_files <- commandArgs(trailingOnly = TRUE)
purrr::map(rmd_files, 
           ~ rmarkdown::render(., 
                               output_format='html_document', 
                               output_file=sub("Rmd$", "html", .)))