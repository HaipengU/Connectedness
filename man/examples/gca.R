# Load cattle data
data(GCcattle)

# Compute genomic relationship matrix
G <- computeG(cattle.W, maf = 0.05, impute = 'mean', method = 'G1')

# Eigendecomposition of genomic relationship matrix
EVD <- eigen(G)

# Estimate variance component
var <- varcomp(y = cattle.pheno$Phenotype, Evector = EVD$vectors, Evalue = EVD$values) 

# Design matrix of fixed effects
## one unit effect
X1 <- model.matrix(~ -1 + factor(cattle.pheno$Unit))
## unit effect and sex effect
X2 <- model.matrix(~ -1 + factor(cattle.pheno$Unit) + factor(cattle.pheno$Sex))

# Calculate CD_IdAve
CD_IdAve <- gcm(Kmatrix = G, Xmatrix = X1, sigma2a = var$Vu, sigma2e = var$Ve, 
                MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CD_IdAve', 
                NumofMU = 'Overall')

# Calculate CDVED1
CDVED1 <- gcm(Kmatrix = G, Xmatrix = X1, sigma2a = var$Vu, sigma2e = var$Ve, 
              MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CDVED1',
              NumofMU = 'Pairwise', diag = TRUE)

# Calculate CDVED2
CDVED2 <- gcm(Kmatrix = G, Xmatrix = X2, sigma2a = var$Vu, sigma2e = var$Ve,
              MUScenario = as.factor(cattle.pheno$Unit), statistic = 'CDVED2', 
              NumofMU = 'Pairwise', Uidx = 5, diag = TRUE)



