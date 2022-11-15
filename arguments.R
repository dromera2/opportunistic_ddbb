GRID.SIZE<- 50
samples = 100
nchains = 2
nParallel = 2
thin= 100
transient = 80*thin
p= c(0.10, 0.25, 0.50, 0.75, 0.90) 
sp_names = paste0("sp", 1:10)
set.seed(1)
speciesCor <- cov2cor(solve(rWishart(1, 22, diag(10))[, , 1]))
i <- c(1:10)
pattern_colnames <- c("var1", "var2", "prob", "sp-sp", "sample")