######### two class unpaired comparison
# y must take values 1,2

library(samr)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.14")
set.seed(100)

######### two class unpaired comparison, imam poissonovo porazdeljene in naredim dve grupi. spreminjam prag in gledam kako se fdr spreminja

#funkcija ki za prvih 10 stolpcev da pois(100) ostalih 10 pois(n)
poison.funkcija <- function(n){

  mu <- matrix(100, 1000, 20)
  mu[1:100, 11:20] <- n
  mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
  x <- matrix(rpois(length(mu), mu), 1000, 20)
  y <- c(rep(1, 10), rep(2, 10))
  samfit <- SAMseq(x, y, resp.type = "Two class unpaired")
# examine significant gene list
  print(samfit)
# plot results
  plot(samfit)}

nbin.funkcija <- function(n){
  mu <- matrix(100, 1000, 20)
  mu[1:100, 11:20] <- n
  size <- 1
  mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
  x <- matrix(rnbinom(length(mu), size=size,mu=mu), 1000, 20)
  y <- c(rep(1, 10), rep(2, 10))
  samfit <- SAMseq(x, y, resp.type = "Two class unpaired")
  # examine significant gene list
  print(samfit)
  # plot results
  plot(samfit)
  
}

#rezultati slikce, lepo vidi da poison sprememba veliko večja kot pri nbinomski

poison.funkcija(100)
poison.funkcija(50)
poison.funkcija(150)

nbin.funkcija(100)
nbin.funkcija(50)
nbin.funkcija(150)

#različni treshold ni narjen se

prag.poison <- function(n,delta){
  mu <- matrix(100, 1000, 20)
  mu[1:100, 11:20] <- n
  mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
  x <- matrix(rpois(length(mu), mu), 1000, 20)
  y <- c(rep(1, 10), rep(2, 10))
  samfit <- SAMseq(x, y, resp.type = "Two class unpaired")
  sig_genes <- summary(samfit, delta=delta)
  
}




# Extend the function to include scatter plot
poison.funkcija <- function(n) {
  mu <- matrix(100, 1000, 20)
  mu[1:100, 11:20] <- n
  mu <- scale(mu, center = FALSE, scale = runif(20, 0.5, 1.5))
  x <- matrix(rpois(length(mu), mu), 1000, 20)
  y <- c(rep(1, 10), rep(2, 10))
  
  # Apply SAMseq
  samfit <- SAMseq(x, y, resp.type = "Two class unpaired")
  
  # Print significant gene list
  print(samfit)
  
  # Plot SAMseq results
  plot(samfit)
  
  # Extracting significant genes
  significant_genes <- samfit$siggenes.table$genes.up$GeneID
  
  # Create data frame for scatter plot
  results <- data.frame(
    Gene_ID = significant_genes,
    Specific_Value = apply(x, 1, mean)[significant_genes],
    Relative_Difference = samfit$siggenes.table$genes.up$Score
  )
  
  # Plot gene-specific scatter plot
  ggplot(results, aes(x = Specific_Value, y = Relative_Difference)) +
    geom_point() +
    geom_text(aes(label = Gene_ID), vjust = -1, hjust = 1) +
    labs(title = "Gene-specific Scatter Plot",
         x = "Gene Specific Value",
         y = "Relative Difference") +
    theme_minimal()
}

# Example usage
poison.funkcija(200)  # Here, 200 is the example value for `n`




