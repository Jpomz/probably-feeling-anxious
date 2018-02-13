# 3_relative_abundance_matrices

# calculate relative abundance matrices for AMD sites

library(plyr)

# function to calculate relative abundance matrices
get_rel_ab <- function(vec, taxa){
  stopifnot(length(vec) == length(taxa))
  rel.ab <- vec / sum(vec)
  Nij <- matrix(0, length(vec), length(vec))
  for (i in 1:length(vec)){
    for (j in 1:length(vec)){
      Nij[i,j] <- rel.ab[i]*rel.ab[j]
    }
  }
  dimnames(Nij) <- list(taxa, taxa)
  Nij
}

# read in data
dat <- readRDS("data/AMD fish invert dw abundance.RDS")
# sort data by avg.dw so matrices are size ordered below
dat <- dat[with(dat, order(dat$avg.dw)),]
# split into list by site
dat.list <- split(dat, list(dat$site))

rel.ab <- llply(dat.list, function (x){
  get_rel_ab(vec = x$rel.ab, taxa = x$taxa)
})

saveRDS(rel.ab, "data/AMD relative abundance matrices.RDS")
