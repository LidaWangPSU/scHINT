# Builds pkgdown/index.md (the site home page): hand-written hero + method figure + model
# cards (pkgdown/hero-*.html, pkgdown/assets/overview.svg) followed by the getting-started
# text knitted from vignettes/getting-started-body.inc. Run from the package root.
library(scHINT)
knitr::knit("pkgdown/index.Rmd", output = "pkgdown/body.md", quiet = TRUE)
rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
fig <- paste0('<div class="sc-figure">\n', rd("pkgdown/assets/overview.svg"), "\n</div>\n")
writeLines(c(rd("pkgdown/hero-top.html"), "", fig, rd("pkgdown/hero-cards.html"), "",
             readLines("pkgdown/body.md")), "pkgdown/index.md")
unlink("pkgdown/body.md")
