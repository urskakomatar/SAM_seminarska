#funkcija simulacije ponovi tolikokrat kot določiš v parametru permutacije in določiš še delto in glede na to dobiš povprečen fdr ven.(poissonova porazdelitev s parametrom 50)

simulacije <- function(ponovitev,delta){
  vsi_fdr <- NULL
  for (i in 1:ponovitev){
    n <- 50
    mu <- matrix(100, 1000, 20)
    mu[1:100, 11:20] <- n
    mu <- scale(mu, center=FALSE, scale=runif(20, 0.5, 1.5))
    x <- matrix(rpois(length(mu), mu), 1000, 20)
    y <- c(rep(1, 10), rep(2, 10))
    
    data=list(x=x,y=y, geneid=as.character(1:nrow(x)),
              genenames=paste("g",as.character(1:nrow(x)),sep=""), logged2=TRUE)
    
    samfit <- samr(data, resp.type = "Two class unpaired")
    delta.table <- samr.compute.delta.table(samfit)
    siggenes.table<-samr.compute.siggenes.table(samfit,delta, data, delta.table)
    
    true_positive_genes <- 1:100  # significant
    false_positive_genes <- 101:1000  # not significant
    
    # total Positives
    total_positives <- siggenes.table$ngenes.lo + siggenes.table$ngenes.up
    
    identified_genes_lo <- as.numeric(siggenes.table$genes.lo[, "Row"])
    identified_genes_up <- as.numeric(siggenes.table$genes.up[, "Row"])
    
    # false positives
    false_positives <- sum(identified_genes_lo > 100) + sum(identified_genes_up > 100)
    
    if (total_positives == 0) {
      fdr <- 0
    } else {
      fdr <- false_positives / total_positives
    }
    vsi_fdr <- append(vsi_fdr,fdr)
    # calculated FDR
  }
  return(vsi_fdr)
  
}


#dela??? :D
#d01 <- simulacije(500, 0.1)
#d025 <- simulacije(500, 0.25)
#d05 <- simulacije(500, 0.5)
#d1 <- simulacije(500, 1)
#d15 <- simulacije(500, 1.5)
#d2 <- simulacije(500, 2)
#d25 <- simulacije(500, 2.5)



#write.csv(d01, file = "d01.csv", row.names = FALSE)
#write.csv(d025, file = "d025.csv", row.names = FALSE)
#write.csv(d05, file = "d05.csv", row.names = FALSE)
#write.csv(d1, file = "d1.csv", row.names = FALSE)
#write.csv(d15, file = "d15.csv", row.names = FALSE)
#write.csv(d2, file = "d2.csv", row.names = FALSE)
#write.csv(d25, file = "d25.csv", row.names = FALSE)
