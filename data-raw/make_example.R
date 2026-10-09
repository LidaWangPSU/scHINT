# Builds data/schint_example.rda: fully simulated genotypes (LD via a Gaussian
# copula) and simulated cell-level / sample-level expression. No real individual-level
# data is distributed with the package.
geno <- simulate_geno(400, 250, seed = 2026)
snp <- data.frame(snp = colnames(geno), chr = 22L, pos = 30000000L + 3000L * (seq_len(ncol(geno)) - 1L))
sim <- simulate_schint(geno, seed = 2026)
schint_example <- list(geno = geno, snp = snp, cell = sim$cell, sample = sim$sample, truth = sim$truth)
save(schint_example, file = "../data/schint_example.rda", compress = "xz")
