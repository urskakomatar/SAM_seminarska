######### two class unpaired comparison
# y must take values 1,2

library(samr)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.14")


######### two class unpaired comparison
set.seed(100)
mu <- matrix(100, 1000, 20)
mu[1:100, 11:20] <- 200
mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
x <- matrix(rpois(length(mu), mu), 1000, 20)
y <- c(rep(1, 10), rep(2, 10))
samfit <- SAMseq(x, y, resp.type = "Two class unpaired")
# examine significant gene list
print(samfit)
# plot results
plot(samfit)

#########3
# Generiranje simuliranih podatkov
n_genes <- 1000  # Število genov
n_samples <- 24  # Število vzorcev (12 zdravih + 12 bolnih)

# Naključni izrazi za gene
exprs <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)

# Generiranje oznak skupin (1 = zdravi, 2 = bolni)
group_labels <- rep(1:2, each = n_samples/2)

# Priprava podatkov za SAM analizo
data <- list(x = exprs, y = group_labels, geneid = as.character(1:n_genes), genenames = paste("Gene", 1:n_genes))

samr_obj <- samr(data, resp.type = "Two class unpaired", nperms = 100)
# Prikaz rezultatov
delta <- 0.3  # Prag za FDR