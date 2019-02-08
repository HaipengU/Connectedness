##########################################################################
#### Function to calculate genomic relationship matrix (G2; VanRaden) ####
##########################################################################
# Imputation + quality control; Arguments: raw Wmatrix and maf (minor allele freq)
computeG2 <- function(Wmatrix, maf) {
    set.seed(1213)
    p1 <- (colMeans(Wmatrix, na.rm = T)/2)
    # imputation for missing markers with rbinom
    for (j in 1:ncol(Wmatrix)){
      Wmatrix[,j] <- ifelse(is.na(Wmatrix[,j]), rbinom(1,2,p1[j]), Wmatrix[,j])
    }
    # remove snp with less than 0.05 MAF
    p2 <- colMeans(Wmatrix)/2
    maf2 <- pmin(p2,1-p2)
    maf.index <- which(maf2< maf)
    W <- Wmatrix[,-maf.index]#
    W_cs <- scale(W)
    Gmatrix <- tcrossprod(W_cs)/ncol(W_cs)
    return(Gmatrix)
} 

#######################################################
####### Function to estimate variance components ######
#######################################################
## Estimate variance componets with rrBLUP package
var_est <- function(y,K,X){
  library(rrBLUP)
  GBLUP <- mixed.solve(y = y, K = K, X = X)
  varcomp <- list(sigma2a = GBLUP$Vu, sigma2e = GBLUP$Ve)
  return(varcomp)
}

#############################################################################
#### Function to calculate connectedness with statistics of PEVD, CD and r ##
#############################################################################
GCfunc <- function(Kmatrix, Xmatrix, sigma2a, sigma2e, MUScenario, statistic, NumofMU){
  # Check if the management unit is factor
  if (is.factor(MUScenario) != TRUE) stop("Management unit is not factor!")
  # Incidence matrix of breeding value
  Zi <- diag(x = 1, nrow = nrow(Kmatrix),ncol = nrow(Kmatrix))
  Zitr <- t(Zi)
  Kinv <- solve(Kmatrix)
  lamda <- sigma2e/sigma2a
  # Design matrix for fixed effects (e.g., management unit)
  X <- Xmatrix
  Mabs <- diag(nrow(X)) - X %*% solve(crossprod(X)) %*% t(X)
  # C_22 matrix
  CuuK = solve(Zitr %*% Mabs %*% Zi + Kinv*lamda)
  # PEV matrix
  PEVK = CuuK*sigma2e
  # management unit
  Management_Unit <- unique(MUScenario)
  if (statistic == 'PEVD'){
    # set contrast matrix for PEVD matrix and full siblings overlap matrix (Laloe contrast method)
    PEVD.contrast <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
    colnames(PEVD.contrast)<- rownames(PEVD.contrast)<- Management_Unit
    for (i in 1:(length(Management_Unit)-1)) {
      for(j in (i+1):length(Management_Unit)) {
        xcontrast <- matrix(0, ncol = 1, nrow = ncol(Kmatrix))
        indexi <- which(MUScenario == Management_Unit[i])
        indexj <- which(MUScenario == Management_Unit[j])
        xcontrast[indexi] <-  1/length(which(MUScenario == Management_Unit[i]))
        xcontrast[indexj] <- -1/length(which(MUScenario == Management_Unit[j]))
        # PEVD matrix across management units
        PEVD.contrast[i,j]<- PEVD.contrast[j,i] <- (crossprod(xcontrast,PEVK) %*% xcontrast)/sigma2a
      }
    }
    # Overall PEVD across all pairwise PEVD 
    PEVD.across.contrast.overall <- mean(PEVD.contrast[upper.tri(PEVD.contrast)])
    if (NumofMU == 'Pairwise'){
      return(PEVD.contrast)
    } else if (NumofMU == 'Overall'){
      return(PEVD.across.contrast.overall)
    }
  } else if (statistic == 'CD'){
    # set contrast matrix for CD matrix and full siblings overlap matrix (Laloe contrast method)
    CD.contrast <-  matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
    colnames(CD.contrast) <- rownames(CD.contrast) <- Management_Unit
    for (i in 1:(length(Management_Unit)-1)) {
      for(j in (i+1):length(Management_Unit)) {
        xcontrast <- matrix(0, ncol = 1, nrow = ncol(Kmatrix))
        indexi <- which(MUScenario == Management_Unit[i])
        indexj <- which(MUScenario == Management_Unit[j])
        xcontrast[indexi] <- 1/length(which(MUScenario == Management_Unit[i]))
        xcontrast[indexj] <- -1/length(which(MUScenario == Management_Unit[j]))
        CD.contrast[i,j] <- CD.contrast[j,i] <- 1 - (lamda* (crossprod(xcontrast,CuuK) %*% xcontrast)/(crossprod(xcontrast,Kmatrix) %*% xcontrast))
      }
    }
    # Overall CD across all pairwise CD
    CD.across.contrast.overall <- mean(CD.contrast[upper.tri(CD.contrast)])
    if (NumofMU == 'Pairwise'){
      return(CD.contrast)
    } else if (NumofMU == 'Overall'){
      return(CD.across.contrast.overall)
    }
  } else if (statistic == 'r'){
    # convert PEV to correlation matrix. r is calculated with individual-wise contrast (Kennedy and Trus)
    rijK<- cov2cor(PEVK)
    # set contrast matrix for rij matrix and full siblings overlap matrix
    rij.contrast <- matrix(NA,ncol=length(Management_Unit),nrow=length(Management_Unit))
    colnames(rij.contrast) <- rownames(rij.contrast) <- Management_Unit
    for (i in 1:(length(Management_Unit)-1)) {
      for(j in (i+1):length(Management_Unit)) {
        indexi <- which(MUScenario == Management_Unit[i])
        indexj <- which(MUScenario == Management_Unit[j])
        # rij matrix across management unit
        rij.Aacross <- rijK[indexi,indexj]
        rij.contrast[i,j] <- rij.contrast[j,i] <- mean(rij.Aacross, na.rm = TRUE)
      }
    }
    # Overall r across all pairwise r
    rij.across.contrast.overall <- mean(rij.contrast[upper.tri(rij.contrast)]) 
    if (NumofMU == 'Pairwise'){
      return(rij.contrast)
    } else if (NumofMU == 'Overall'){
      return(rij.across.contrast.overall)
    }
  }
}







