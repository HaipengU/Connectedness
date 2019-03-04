#' Measurement of connectedness.
#'
#' The estimates of connectedness across management units.
#' 
#' @param Kmatrix a n by n relationship matrix. 
#' @param Xmatrix a design matrix which associates fixed effects with phenotypes. 
#' @param sigma2a additive genetic variance.
#' @param sigma2e residual variance.
#' @param MUScenario a vector of managment units which will be treatd as a factor. 
#' @param statistic a statistic which will be used to measure connectedness.
#' @param NumofMU number of management unit used to calculate connectedness. 
#' 
#' @return 
#' A value of overall connectedness measurements across management units when NumofMU is set as 'Overall'.
#' A matrix of connectedness measurments with diagnol as NA when when NumofMU is set as 'Pairwise'.
#' 
#' @examples 
#' GC()
#' 
#' @export
#' 
GC <- function(Kmatrix, Xmatrix, sigma2a, sigma2e, MUScenario, statistic, NumofMU){
  if (is.factor(MUScenario) != TRUE) stop("Management unit is not factor!")
  Zi <- diag(x = 1, nrow = nrow(Kmatrix),ncol = nrow(Kmatrix))
  Zitr <- t(Zi)
  Kinv <- solve(Kmatrix)
  lamda <- sigma2e/sigma2a
  X <- Xmatrix
  Mabs <- diag(nrow(X)) - X %*% solve(crossprod(X)) %*% t(X)
  CuuK = solve(Zitr %*% Mabs %*% Zi + Kinv*lamda)
  PEVK = CuuK * sigma2e
  Management_Unit <- unique(MUScenario)
  if (statistic == 'PEVD'){
    PEVD.contrast <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
    colnames(PEVD.contrast) <- rownames(PEVD.contrast) <- Management_Unit
    for (i in 1 : (length(Management_Unit) - 1)) {
      for(j in (i+1) : length(Management_Unit)) {
        xcontrast <- matrix(0, ncol = 1, nrow = ncol(Kmatrix))
        indexi <- which(MUScenario == Management_Unit[i])
        indexj <- which(MUScenario == Management_Unit[j])
        xcontrast[indexi] <-  1 / length(which(MUScenario == Management_Unit[i]))
        xcontrast[indexj] <- -1 / length(which(MUScenario == Management_Unit[j]))
        PEVD.contrast[i,j] <- PEVD.contrast[j,i] <- (crossprod(xcontrast,PEVK) %*% xcontrast) / sigma2a
      }
    }
    if (NumofMU == 'Pairwise') {
      return(PEVD.contrast)
    } else if (NumofMU == 'Overall') {
      PEVD.contrast.overall <- mean(PEVD.contrast[upper.tri(PEVD.contrast)])
      return(PEVD.contrast.overall)
    }
  } else if (statistic == 'CD') {
    CD.contrast <-  matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
    colnames(CD.contrast) <- rownames(CD.contrast) <- Management_Unit
    for (i in 1 : (length(Management_Unit) - 1)) {
      for(j in (i + 1) : length(Management_Unit)) {
        xcontrast <- matrix(0, ncol = 1, nrow = ncol(Kmatrix))
        indexi <- which(MUScenario == Management_Unit[i])
        indexj <- which(MUScenario == Management_Unit[j])
        xcontrast[indexi] <- 1 / length(which(MUScenario == Management_Unit[i]))
        xcontrast[indexj] <- -1 / length(which(MUScenario == Management_Unit[j]))
        CD.contrast[i,j] <- CD.contrast[j,i] <- 1 - (lamda* (crossprod(xcontrast,CuuK) %*% xcontrast) / (crossprod(xcontrast,Kmatrix) %*% xcontrast))
      }
    }
    if (NumofMU == 'Pairwise') {
      return(CD.contrast)
    } else if (NumofMU == 'Overall') {
      CD.contrast.overall <- mean(CD.contrast[upper.tri(CD.contrast)])
      return(CD.contrast.overall)
    }
  } else if (statistic == 'r1') {
    rijK <- cov2cor(PEVK)
    rij.contrast <- matrix(NA, ncol=length(Management_Unit), nrow=length(Management_Unit))
    colnames(rij.contrast) <- rownames(rij.contrast) <- Management_Unit
    for (i in 1 : (length(Management_Unit) - 1)) {
      for(j in (i + 1) : length(Management_Unit)) {
        indexi <- which(MUScenario == Management_Unit[i])
        indexj <- which(MUScenario == Management_Unit[j])
        rij.Aacross <- rijK[indexi, indexj]
        rij.contrast[i, j] <- rij.contrast[j, i] <- mean(rij.Aacross, na.rm = TRUE)
      }
    }
    if (NumofMU == 'Pairwise'){
      return(rij.contrast)
    } else if (NumofMU == 'Overall'){
      rij.contrast.overall <- mean(rij.contrast[upper.tri(rij.contrast)])
      return(rij.contrast.overall)
    }
  } else if (statistic == 'r2'){
    rij.contrast <- matrix(NA, ncol=length(Management_Unit),nrow=length(Management_Unit))
    colnames(rij.contrast) <- rownames(rij.contrast) <- Management_Unit
    for (i in 1 : (length(Management_Unit)-1)) {
      for(j in (i + 1) : length(Management_Unit)) {
        indexi <- which(MUScenario == Management_Unit[i])
        indexj <- which(MUScenario == Management_Unit[j])
        rij.contrast[i, j] <- rij.contrast[j, i]<- 
          sum(PEVK[indexi, indexj]) / sqrt(sum(PEVK[indexi, indexi]) * sum(PEVK[indexj, indexj]))
      }
    }
    if (NumofMU == 'Pairwise'){
      return(rij.contrast)
    } else if (NumofMU == 'Overall'){
      rij.contrast.overall <- mean(rij.contrast[upper.tri(rij.contrast)]) 
      return(rij.contrast.overall)
    }
  } else if (statistic == 'VED') {
    C22_inv <- solve(crossprod(Zi) + Kinv * lamda)
    var_Bhat <- solve(crossprod(X) - (crossprod(X, Zi) %*% C22_inv %*% crossprod(Zi, X))) * sigma2e
    VED <- matrix(NA, ncol = ncol(var_Bhat), nrow = nrow(var_Bhat))
    for (i in 1 : ncol(VED) - 1) {
      for (j in (i + 1) : ncol(VED)) {
        VED[i, j] <- VED[j, i] <- (var_Bhat[i, i] + var_Bhat[j, j] - 2 * var_Bhat[i, j])
      }
    }
    if (NumofMU == 'Pairwise') {
      return(VED)
    } else if (NumofMU == 'Overall') {
      VED.overall <- mean(VED[upper.tri(VED)])
      return(VED.overall)
    } 
  } else if (statistic == 'VED2') {
    C22_inv <- solve(crossprod(Zi) + Kinv * lamda)
    var_Bhat <- (solve(crossprod(X) - (crossprod(X, Zi) %*% C22_inv %*% crossprod(Zi, X))) - solve(crossprod(X))) * sigma2e
    VED <- matrix(NA, ncol = ncol(var_Bhat), nrow = nrow(var_Bhat))
    for (i in 1 : ncol(VED) - 1) {
      for (j in (i + 1) : ncol(VED)) {
        VED[i, j] <- VED[j, i] <- (var_Bhat[i, i] + var_Bhat[j, j] - 2 * var_Bhat[i, j])
      }
    }
    if (NumofMU == 'Pairwise') {
      return(VED)
    } else if (NumofMU == 'Overall') {
      VED.overall <- mean(VED[upper.tri(VED)])
      return(VED.overall)
    } 
  } else if (statistic == 'VED3') {
    C22_inv <- solve(crossprod(Zi) + Kinv * lamda)
    m1 <- dim(Xmatrix)[2]
    m2 <- length(unique(MUScenario))
    X1 <- X[, 1 : m2]
    X2 <- X[, -c(1 : m2)]
    var_Bhat <- (solve(crossprod(X) - (crossprod(X, Zi) %*% C22_inv %*% crossprod(Zi, X))) - solve(crossprod(X))) * sigma2e
    PEVmean <- solve(crossprod(X1)) %*% crossprod(X1, X2) %*% var_Bhat[(m2+1) : m1, (m2+1) : m1] %*% crossprod(X2, X1) %*% 
                solve(crossprod(X1)) + solve(crossprod(X1)) %*% crossprod(X1, X2) %*% var_Bhat[(m2 + 1) : m1, 1 : m2] +
                var_Bhat[1 : m2, (m2 + 1) : m1] %*% crossprod(X2, X1) %*% solve(crossprod(X1)) +
                var_Bhat[1 : m2, 1 : m2] - sigma2e * solve(crossprod(X1))
    VED <- matrix(NA, ncol = ncol(PEVmean), nrow = nrow(PEVmean))
    for (i in 1 : ncol(VED) - 1) {
      for (j in (i + 1) : ncol(VED)) {
        VED[i, j] <- VED[j, i] <- (PEVmean[i, i] + PEVmean[j, j] - 2 * PEVmean[i, j])
      }
    }
    if (NumofMU == 'Pairwise') {
      return(VED)
    } else if (NumofMU == 'Overall') {
      VED.overall <- mean(VED[upper.tri(VED)])
      return(VED.overall)
    } 
  } else if (statistic == 'CR') {
    C22_inv <- solve(crossprod(Zi) + Kinv * lamda)
    var_Bhat <- solve(crossprod(X) - (crossprod(X, Zi) %*% C22_inv %*% crossprod(Zi, X))) * sigma2e
    CR <- matrix(NA, ncol = ncol(var_Bhat), nrow = nrow(var_Bhat))
    for (i in 1 : ncol(CR) - 1){
      for (j in (i + 1) : ncol(CR)){
        CR[i, j] <- CR[j, i] <- (var_Bhat[i, j]/sqrt(var_Bhat[i, i] * var_Bhat[j, j]))
      }
    }
    if (NumofMU == 'Pairwise') {
      return(CR)
    } else if (NumofMU == 'Overall') {
      CR.overall <- mean(CR[upper.tri(CR)])
      return(CR.overall)
    }
  } else if (statistic == 'CR2'){
    C22_inv <- solve(crossprod(Zi) + Kinv * lamda)
    var_Bhat <- (solve(crossprod(X) - (crossprod(X, Zi) %*% C22_inv %*% crossprod(Zi, X))) - solve(crossprod(X))) * sigma2e
    CR <- matrix(NA, ncol = ncol(var_Bhat), nrow = nrow(var_Bhat))
    for (i in 1 : ncol(CR) - 1) {
      for (j in (i + 1) : ncol(CR)) {
        CR[i, j] <- CR[j, i] <- (var_Bhat[i, j] / sqrt(var_Bhat[i, i] * var_Bhat[j, j]))
      }
    } 
    if (NumofMU == 'Pairwise') {
      return(CR)
    } else if (NumofMU == 'Overall') {
      CR.overall <- mean(CR[upper.tri(CR)])
      return(CR.overall)
    } 
  } else if (statistic == 'CR3'){
    C22_inv <- solve(crossprod(Zi) + Kinv * lamda)
    m1 <- dim(Xmatrix)[2]
    m2 <- length(unique(MUScenario))
    X1 <- X[, 1 : m2]
    X2 <- X[, -c(1 : m2)]
    var_Bhat <- (solve(crossprod(X) - (crossprod(X, Zi) %*% C22_inv %*% crossprod(Zi, X))) - solve(crossprod(X))) * sigma2e
    PEVmean <- solve(crossprod(X1)) %*% crossprod(X1, X2) %*% var_Bhat[(m2+1) : m1, (m2+1) : m1] %*% crossprod(X2, X1) %*% 
      solve(crossprod(X1)) + solve(crossprod(X1)) %*% crossprod(X1, X2) %*% var_Bhat[(m2 + 1) : m1, 1 : m2] +
      var_Bhat[1 : m2, (m2+1) : m1] %*% crossprod(X2, X1) %*% solve(crossprod(X1)) +
      var_Bhat[1 : m2, 1 : m2] - sigma2e * solve(crossprod(X1))
    CR <- matrix(NA, ncol = ncol(PEVmean), nrow = nrow(PEVmean))
    for (i in 1 : ncol(CR) - 1) {
      for (j in (i + 1) : ncol(CR)) {
        CR[i, j] <- CR[j, i] <- PEVmean[i, j] / sqrt(PEVmean[i, i] * PEVmean[j, j])
      }
    }
    if (NumofMU == 'Pairwise') {
      return(CR)
    } else if (NumofMU == 'Overall') {
      CR.overall <- mean(CR[upper.tri(CR)])
      return(CR.overall)
    } 
  }
}

