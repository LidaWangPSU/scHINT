# Builds pkgdown/index.md (the site home page): hand-written hero (pkgdown/hero-top.html,
# containing the overview figure pkgdown/assets/schint-overview.png) and model cards
# (pkgdown/hero-cards.html), followed by the getting-started text knitted from
# vignettes/getting-started-body.inc. Run from the package root.
library(scHINT)
knitr::knit("pkgdown/index.Rmd", output = "pkgdown/body.md", quiet = TRUE)
rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
writeLines(c(rd("pkgdown/hero-top.html"), "", rd("pkgdown/hero-cards.html"), "",
             readLines("pkgdown/body.md")), "pkgdown/index.md")
unlink("pkgdown/body.md")
