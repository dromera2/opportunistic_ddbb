GRID.SIZE<- 50 #Grid size for artificial data generation
samples = 100 #Sample size for hmsc calibration
nchains = 2 #number of Markov chain Montecarlo for hmsc calibration
nParallel = 2 
thin= 100 #thinning factor for hmsc calibration
transient = 80*thin  #warm-up value for hmsc
p= c(0.10, 0.25, 0.50, 0.75, 0.90) #percentages for subsampling
sp_names = paste0("sp", 1:10) #species names
set.seed(1)
speciesCor <- cov2cor(solve(rWishart(1, 22, diag(10))[, , 1])) #omega parameters for generate artificial data (2 block)
i <- c(1:10) #number of replicates
pattern_colnames <- c("var1", "var2", "prob", "sp-sp", "sample")