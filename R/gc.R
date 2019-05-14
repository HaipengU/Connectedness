#' Measurement of connectedness.
#'
#' The estimates of connectedness across management units using pedigree or genomic data.
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
#' gc()
#' 
#' @export
#' 
gc <- function(Kmatrix, Xmatrix, sigma2a, sigma2e, MUScenario, statistic, NumofMU) {
  if(is.factor(MUScenario) != TRUE) stop("Management unit is not factor!")
  Zi <- diag(x = 1, nrow = nrow(Kmatrix), ncol = nrow(Kmatrix))
  Zitr <- t(Zi)
  diag(Kmatrix) <- diag(Kmatrix) + 0.00001
  Kinv <- solve(Kmatrix)
  lamda <- sigma2e/sigma2a
  X <- Xmatrix
  Mabs <- diag(nrow(X)) - X %*% solve(crossprod(X)) %*% t(X)
  CuuK = solve(Zitr %*% Mabs %*% Zi + Kinv*lamda)
  PEVK = CuuK * sigma2e
  Management_Unit <- unique(MUScenario)
  switch(statistic,
         PEVD_IdAve = pevd.idAve(Management_Unit = Management_Unit, Kmatrix = Kmatrix, 
                                 MUScenario = MUScenario, PEVK = PEVK, NumofMU = NumofMU),
         PEVD_GrpAve = pevd.grpAve(Management_Unit = Management_Unit, Kmatrix = Kmatrix, 
                                   MUScenario = MUScenario, PEVK = PEVK, NumofMU = NumofMU),
         PEVD_contrast = gc.contrast(statistic = statistic, Management_Unit = Management_Unit, MUScenario = MUScenario, Kmatrix = Kmatrix,
                                     sigma2a = sigma2a, sigma2e = sigma2e, CuuK = CuuK, NumofMU = NumofMU),
         VED = pevd.cor (statistic = statistic, NumofMU = NumofMU , Zi = Zi, Kinv = Kinv, lamda = lamda,
                         X = X, sigma2e = sigma2e, Management_Unit = Management_Unit),
         PEVD.cor.single = pevd.cor (statistic = statistic, NumofMU = NumofMU, Zi = Zi, Kinv = Kinv, lamda = lamda,
                                     X = X, sigma2e = sigma2e, Management_Unit = Management_Unit),
         PEVD.cor.multiple = pevd.cor (statistic = statistic, NumofMU = NumofMU, Zi = Zi, Kinv = Kinv, lamda = lamda,
                                       X = X, sigma2e = sigma2e, Management_Unit = Management_Unit),
         CD_IdAve = cd.idAve(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit, 
                             MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a),
         CD_GrpAve = cd.grpAve(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit, 
                               MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a),
         CD_contrast = gc.contrast(statistic = statistic, Management_Unit = Management_Unit, MUScenario = MUScenario, Kmatrix = Kmatrix,
                                   sigma2a = sigma2a, sigma2e = sigma2e, CuuK = CuuK, NumofMU = NumofMU),
         CD.approx = cd.cor(statistic = statistic, NumofMU, Zi = Zi, Kinv = Kinv, Kmatrix = Kmatrix, lamda = lamda,
                            X = X, sigma2e = sigma2e, sigma2a =sigma2a, Management_Unit = Management_Unit, MUScenario = MUScenario),
         CD.cor.single = cd.cor(statistic = statistic, NumofMU, Zi = Zi, Kinv = Kinv, Kmatrix = Kmatrix, lamda = lamda,
                                X = X, sigma2e = sigma2e, sigma2a =sigma2a, Management_Unit = Management_Unit, MUScenario = MUScenario),
         CD.cor.multiple = cd.cor(statistic = statistic, NumofMU, Zi = Zi, Kinv = Kinv, Kmatrix = Kmatrix, lamda = lamda,
                                  X = X, sigma2e = sigma2e, sigma2a =sigma2a, Management_Unit = Management_Unit, MUScenario = MUScenario),
         r_IdAve = r.idAve(PEVK = PEVK, Management_Unit = Management_Unit, MUScenario = MUScenario, NumofMU = NumofMU),
         r_GrpAve = r.grpAve(PEVK = PEVK, Kmatrix = Kmatrix, Management_Unit = Management_Unit, 
                             MUScenario = MUScenario, NumofMU = NumofMU),
         r_contrast = gc.contrast(statistic = statistic, Management_Unit = Management_Unit, MUScenario = MUScenario, Kmatrix = Kmatrix,
                                  sigma2a = sigma2a, sigma2e = sigma2e, CuuK = CuuK, NumofMU = NumofMU),
         CR = r.cor (statistic = statistic, NumofMU, Zi = Zi, Kinv = Kinv, lamda = lamda,
                     X = X, sigma2e = sigma2e, Management_Unit = Management_Unit),
         r.cor.single = r.cor (statistic = statistic, NumofMU, Zi = Zi, Kinv = Kinv, lamda = lamda,
                               X = X, sigma2e = sigma2e, Management_Unit = Management_Unit),
         r.cor.multiple = r.cor (statistic = statistic, NumofMU, Zi = Zi, Kinv = Kinv, lamda = lamda,
                                 X = X, sigma2e = sigma2e, Management_Unit = Management_Unit)
  )
}


# Contrast for PEVD, CD & r
gc.contrast <- function(statistic = statistic, Management_Unit = Management_Unit, MUScenario = MUScenario, Kmatrix = Kmatrix,
                        sigma2a = sigma2a, sigma2e = sigma2e, CuuK = CuuK, NumofMU = NumofMU){
  Contrast_Matrix <-  matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
  colnames(Contrast_Matrix) <- rownames(Contrast_Matrix) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i + 1) : length(Management_Unit)) {
      xcontrast <- matrix(0, ncol = 1, nrow = ncol(Kmatrix))
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      xcontrast[indexi] <- 1 / length(which(MUScenario == Management_Unit[i]))
      xcontrast[indexj] <- -1 / length(which(MUScenario == Management_Unit[j]))
      if (statistic == 'PEVD_contrast'){
        PEVK <- CuuK * sigma2e
        Contrast_Matrix[i,j] <- Contrast_Matrix[j,i] <- (crossprod(xcontrast, PEVK) %*% xcontrast) 
        
      } else if (statistic == 'CD_contrast'){
        lamda <- sigma2e/sigma2a
        Contrast_Matrix[i,j] <- Contrast_Matrix[j,i] <- 
          1 - (lamda* (crossprod(xcontrast,CuuK) %*% xcontrast) / (crossprod(xcontrast, Kmatrix) %*% xcontrast))
      } else if (statistic == 'r_contrast'){
        rijK <- cov2cor(PEVK)
        Contrast_Matrix[i,j] <- Contrast_Matrix[j,i] <- (crossprod(xcontrast, rijK) %*% xcontrast)
      }
    }
  }
  if (NumofMU == 'Pairwise') {
    return(Contrast_Matrix)
  } else if (NumofMU == 'Overall') {
    Contrast_Matrix.overall <- mean(Contrast_Matrix[upper.tri(Contrast_Matrix)])
    return(Contrast_Matrix.overall)
  }
}


# Individual Ave PEVD
pevd.idAve <- function(PEVK = PEVK, Kmatrix = Kmatrix, 
                       Management_Unit = Management_Unit, 
                       MUScenario = MUScenario, NumofMU = NumofMU){
  PEVD.Pairwise.all <- matrix(NA, ncol = ncol(PEVK), nrow = nrow(PEVK))
  for (i in 1:(nrow(Kmatrix)-1)){
    for (j in (i+1):nrow(Kmatrix)){
      PEVD.Pairwise.all[i,j] <- PEVD.Pairwise.all[j,i] <- (PEVK[i,i] + PEVK[j,j] - 2*PEVK[i,j])
    }
  }
  PEVD.Pairwise.unit <- matrix(NA, ncol = length(Management_Unit), 
                               nrow = length(Management_Unit))
  colnames(PEVD.Pairwise.unit) <- rownames(PEVD.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i+1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      PEVD.Pairwise.unit[i, j] <- PEVD.Pairwise.unit[j, i] <- 
        mean(PEVD.Pairwise.all[indexi, indexj])
    }
  }
  if (NumofMU == 'Pairwise') {
    return(PEVD.Pairwise.unit)
  } else if (NumofMU == 'Overall') {
    PEVD.IdAve.overall <- mean(PEVD.Pairwise.unit[upper.tri(PEVD.Pairwise.unit)])
    return(PEVD.IdAve.overall)
  }
}


# Group Ave PEVD
pevd.grpAve <- function(PEVK = PEVK, Kmatrix = Kmatrix, 
                        Management_Unit = Management_Unit, 
                        MUScenario = MUScenario, NumofMU = NumofMU){
  PEVD.Pairwise.unit <- matrix(NA, ncol = length(Management_Unit), 
                               nrow = length(Management_Unit))
  colnames(PEVD.Pairwise.unit) <- rownames(PEVD.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i+1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      PEV.ii <- PEVK[indexi, indexi]
      PEV.jj <- PEVK[indexj, indexj]
      PEC.ij <- mean(PEVK[indexi, indexj])
      PEVD.Pairwise.unit[i, j] <- PEVD.Pairwise.unit[j, i] <- mean(PEV.ii[upper.tri(PEV.ii)]) +  
        mean(PEV.jj[upper.tri(PEV.jj)]) - 2* PEC.ij
    }
  }
  if (NumofMU == 'Pairwise') {
    return(PEVD.Pairwise.unit)
  } else if (NumofMU == 'Overall') {
    PEVD.GrpAve.overall <- mean(PEVD.Pairwise.unit[upper.tri(PEVD.Pairwise.unit)])
    return(PEVD.GrpAve.overall)
  }
}


# VED; PEVD.cor.single; PEVD.cor.multiple (Xmatrix does not include intercept)
pevd.cor <- function(statistic = statistic, NumofMU, Zi = Zi, Kinv = Kinv, lamda = lamda,
                     X = X, sigma2e = sigma2e, Management_Unit = Management_Unit){
  C22_inv <- solve(crossprod(Zi) + Kinv * lamda)
  var_Bhat <- solve(crossprod(X) - (crossprod(X, Zi) %*% C22_inv %*% crossprod(Zi, X))) * sigma2e
  Summary.Matrix <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
  colnames(Summary.Matrix) <- rownames(Summary.Matrix) <- Management_Unit
  if(statistic == 'VED') {
    PEV_mean <- var_Bhat
  } else if (statistic == 'PEVD.cor.single'){
    cor_factor <- solve(crossprod(X)) * sigma2e
    PEV_mean <- var_Bhat - cor_factor
  } else if (statistic == 'PEVD.cor.multiple'){
    m <- dim(X)[2] # length of total fixed effects
    m1 <- length(Management_Unit) # lengt of MU
    X1 <- X[, 1 : m1] # design matrix for MU
    X2 <- X[, -c(1 : m1)]
    PEV_mean <- solve(crossprod(X1)) %*% crossprod(X1, X2) %*% var_Bhat[(m1 + 1) : m, (m1 + 1) : m] %*% 
      crossprod(X2, X1) %*% solve(crossprod(X1)) + solve(crossprod(X1)) %*% crossprod(X1, X2) %*% 
      var_Bhat[(m1 + 1) : m, 1 : m1] + var_Bhat[1 : m1, (m1 + 1) : m] %*% crossprod(X2, X1) %*% 
      solve(crossprod(X1)) + var_Bhat[1 : m1, 1 : m1] - sigma2e * solve(crossprod(X1))
  }
  for (i in 1 : ncol(Summary.Matrix) - 1) {
    for (j in (i + 1) : ncol(Summary.Matrix)) {
      Summary.Matrix[i, j] <- Summary.Matrix[j, i] <- (PEV_mean[i, i] + PEV_mean[j, j] - 2 * PEV_mean[i, j])
    }
  }
  if (NumofMU == 'Pairwise') {
    return(Summary.Matrix)
  } else if (NumofMU == 'Overall') {
    Summary.Matrix.overall <- mean(Summary.Matrix[upper.tri(Summary.Matrix)])
    return(Summary.Matrix.overall)
  }
}


# Individual Average CD
cd.idAve <- function(PEVK = PEVK, Kmatrix = Kmatrix, 
                     Management_Unit = Management_Unit, 
                     MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a){
  PEVD.Pairwise.all <- K.diff.Pairwise.all <- matrix(NA, ncol = ncol(PEVK), nrow = nrow(PEVK))
  for (i in 1:(nrow(Kmatrix)-1)){
    for (j in (i+1):nrow(Kmatrix)){
      PEVD.Pairwise.all[i,j] <- PEVD.Pairwise.all[j,i] <- (PEVK[i,i] + PEVK[j,j] - 2*PEVK[i,j])
      K.diff.Pairwise.all[i,j] <- K.diff.Pairwise.all[j,i] <- (Kmatrix[i,i] + Kmatrix[j,j] - 2*Kmatrix[i,j])
    }
  }
  CD.Pairwise.unit <- matrix(NA, ncol = length(Management_Unit), 
                             nrow = length(Management_Unit))
  colnames(CD.Pairwise.unit) <- rownames(CD.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i+1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      CD.Pairwise.unit[i, j] <- CD.Pairwise.unit[j, i] <- 
        1- (mean(PEVD.Pairwise.all[indexi, indexj])/(mean(K.diff.Pairwise.all[indexi, indexj]) * sigma2a))
    }
  }
  if (NumofMU == 'Pairwise') {
    return(CD.Pairwise.unit)
  } else if (NumofMU == 'Overall') {
    CD.IdAve.overall <- mean(CD.Pairwise.unit[upper.tri(CD.Pairwise.unit)])
    return(CD.IdAve.overall)
  }
}


# Group Average CD
cd.grpAve <- function(PEVK = PEVK, Kmatrix = Kmatrix, 
                      Management_Unit = Management_Unit, 
                      MUScenario = MUScenario, NumofMU = NumofMU, sigma2a = sigma2a){
  CD.Pairwise.unit <- matrix(NA, ncol = length(Management_Unit), 
                             nrow = length(Management_Unit))
  colnames(CD.Pairwise.unit) <- rownames(CD.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i+1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      PEV.ii <- PEVK[indexi, indexi]
      K.ii <- Kmatrix[indexi, indexi]
      PEV.jj <- PEVK[indexj, indexj]
      K.jj <- Kmatrix[indexj, indexj]
      PEC.ij <- mean(PEVK[indexi, indexj])
      K.ij <- mean(Kmatrix[indexi, indexj])
      PEVD.ij <- mean(PEV.ii[upper.tri(PEV.ii)]) + mean(PEV.jj[upper.tri(PEV.jj)]) - 2* PEC.ij
      K.diff.ij <- mean(K.ii[upper.tri(K.ii)]) + mean(K.jj[upper.tri(K.jj)]) - 2* K.ij
      CD.Pairwise.unit[i, j] <- CD.Pairwise.unit[j, i] <- 1 - (PEVD.ij/(K.diff.ij*sigma2a))
    }
  }
  if (NumofMU == 'Pairwise') {
    return(CD.Pairwise.unit)
  } else if (NumofMU == 'Overall') {
    CD.GrpAve.overall <- mean(CD.Pairwise.unit[upper.tri(CD.Pairwise.unit)])
    return(CD.GrpAve.overall)
  }
}


# CD.approx, CD.cor.single; CD.cor.multiple
cd.cor <- function(statistic = statistic, NumofMU, Zi = Zi, Kinv = Kinv, Kmatrix = Kmatrix, lamda = lamda,
                   X = X, sigma2e = sigma2e, sigma2a =sigma2a, Management_Unit = Management_Unit,
                   MUScenario = MUScenario){
  C22_inv <- solve(crossprod(Zi) + Kinv * lamda)
  var_Bhat <- solve(crossprod(X) - (crossprod(X, Zi) %*% C22_inv %*% crossprod(Zi, X))) * sigma2e
  Summary.Matrix <- K.Pairwise.unit.diff <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
  colnames(Summary.Matrix) <- rownames(Summary.Matrix) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i+1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      K.ii <- Kmatrix[indexi, indexi]
      K.jj <- Kmatrix[indexj, indexj]
      K.ij <- mean(Kmatrix[indexi, indexj])
      K.Pairwise.unit.diff[i,j] <- K.Pairwise.unit.diff[j,i] <- (mean(K.ii[upper.tri(K.ii)]) + mean(K.jj[upper.tri(K.jj)]) - 2* K.ij) * sigma2a
    }
  }
  if(statistic == 'CD.approx') {
    PEV_mean <- var_Bhat
  } else if (statistic == 'CD.cor.single'){
    cor_factor <- solve(crossprod(X)) * sigma2e
    PEV_mean <- var_Bhat - cor_factor
  } else if (statistic == 'CD.cor.multiple'){
    m <- dim(X)[2] # length of total fixed effects
    m1 <- length(Management_Unit) # lengt of MU
    X1 <- X[, 1 : m1] # design matrix for MU
    X2 <- X[, -c(1 : m1)]
    PEV_mean <- solve(crossprod(X1)) %*% crossprod(X1, X2) %*% var_Bhat[(m1 + 1) : m, (m1 + 1) : m] %*% 
      crossprod(X2, X1) %*% solve(crossprod(X1)) + solve(crossprod(X1)) %*% crossprod(X1, X2) %*% 
      var_Bhat[(m1 + 1) : m, 1 : m1] + var_Bhat[1 : m1, (m1 + 1) : m] %*% crossprod(X2, X1) %*% 
      solve(crossprod(X1)) + var_Bhat[1 : m1, 1 : m1] - sigma2e * solve(crossprod(X1))
  }
  for (i in 1 : ncol(Summary.Matrix) - 1) {
    for (j in (i + 1) : ncol(Summary.Matrix)) {
      Summary.Matrix[i, j] <- Summary.Matrix[j, i] <- 1 - (
        (PEV_mean[i, i] + PEV_mean[j, j] - 2 * PEV_mean[i, j])/ 
          K.Pairwise.unit.diff[i, j]
      )
    }
  } 
  if (NumofMU == 'Pairwise') {
    return(Summary.Matrix)
  } else if (NumofMU == 'Overall') {
    Summary.Matrix.overall <- mean(Summary.Matrix[upper.tri(Summary.Matrix)])
    return(Summary.Matrix.overall)
  }
}


# Individual r
r.idAve <- function(PEVK = PEVK, Management_Unit = Management_Unit,
                    MUScenario = MUScenario, NumofMU = NumofMU) {
  rijK <- cov2cor(PEVK)
  rij.Pairwise.unit <- matrix(NA, ncol=length(Management_Unit), nrow=length(Management_Unit))
  colnames(rij.Pairwise.unit) <- rownames(rij.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i + 1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      rij.Aacross <- rijK[indexi, indexj]
      rij.Pairwise.unit[i, j] <- rij.Pairwise.unit[j, i] <- mean(rij.Aacross, na.rm = TRUE)
    }
  }
  if (NumofMU == 'Pairwise'){
    return(rij.Pairwise.unit)
  } else if (NumofMU == 'Overall'){
    rij.Pairwise.unit.overall <- mean(rij.Pairwise.unit[upper.tri(rij.Pairwise.unit)])
    return(rij.Pairwise.unit.overall)
  }
}


# Group Average r
r.grpAve <- function(PEVK = PEVK, Kmatrix = Kmatrix, 
                     Management_Unit = Management_Unit, 
                     MUScenario = MUScenario, NumofMU = NumofMU){
  r.Pairwise.unit <- matrix(NA, ncol = length(Management_Unit), 
                            nrow = length(Management_Unit))
  colnames(r.Pairwise.unit) <- rownames(r.Pairwise.unit) <- Management_Unit
  for (i in 1 : (length(Management_Unit) - 1)) {
    for(j in (i+1) : length(Management_Unit)) {
      indexi <- which(MUScenario == Management_Unit[i])
      indexj <- which(MUScenario == Management_Unit[j])
      PEV.ii.mean <- mean(PEVK[indexi, indexi])
      PEV.jj.mean <- mean(PEVK[indexj, indexj])
      PEC.ij.mean <- mean(PEVK[indexi, indexj])
      r.Pairwise.unit[i, j] <- r.Pairwise.unit[j, i] <- 
        PEC.ij.mean/sqrt(PEV.ii.mean * PEV.jj.mean)
    }
  }
  if (NumofMU == 'Pairwise') {
    return(r.Pairwise.unit)
  } else if (NumofMU == 'Overall') {
    r.GrpAve.overall <- mean(r.Pairwise.unit[upper.tri(r.Pairwise.unit)])
    return(r.GrpAve.overall)
  }
}


# CR; r.cor.single; r.cor.multiple
r.cor <- function(statistic = statistic, NumofMU, Zi = Zi, Kinv = Kinv, lamda = lamda,
                  X = X, sigma2e = sigma2e, Management_Unit = Management_Unit){
  C22_inv <- solve(crossprod(Zi) + Kinv * lamda)
  var_Bhat <- solve(crossprod(X) - (crossprod(X, Zi) %*% C22_inv %*% crossprod(Zi, X))) * sigma2e
  Summary.Matrix <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
  colnames(Summary.Matrix) <- rownames(Summary.Matrix) <- Management_Unit
  if(statistic == 'CR') {
    PEV_mean <- var_Bhat
  } else if (statistic == 'r.cor.single'){
    cor_factor <- solve(crossprod(X)) * sigma2e
    PEV_mean <- var_Bhat - cor_factor
  } else if (statistic == 'r.cor.multiple'){
    m <- dim(X)[2] # length of total fixed effects
    m1 <- length(Management_Unit) # lengt of MU
    X1 <- X[, 1 : m1] # design matrix for MU
    X2 <- X[, -c(1 : m1)]
    PEV_mean <- solve(crossprod(X1)) %*% crossprod(X1, X2) %*% var_Bhat[(m1 + 1) : m, (m1 + 1) : m] %*% 
      crossprod(X2, X1) %*% solve(crossprod(X1)) + solve(crossprod(X1)) %*% crossprod(X1, X2) %*% 
      var_Bhat[(m1 + 1) : m, 1 : m1] + var_Bhat[1 : m1, (m1 + 1) : m] %*% crossprod(X2, X1) %*% 
      solve(crossprod(X1)) + var_Bhat[1 : m1, 1 : m1] - sigma2e * solve(crossprod(X1))
  }
  for (i in 1 : ncol(Summary.Matrix) - 1) {
    for (j in (i + 1) : ncol(Summary.Matrix)) {
      Summary.Matrix[i, j] <- Summary.Matrix[j, i] <- PEV_mean[i, j] / sqrt(PEV_mean[i, i] * PEV_mean[j, j])
    }
  }
  if (NumofMU == 'Pairwise') {
    return(Summary.Matrix)
  } else if (NumofMU == 'Overall') {
    Summary.Matrix.overall <- mean(Summary.Matrix[upper.tri(Summary.Matrix)])
    return(Summary.Matrix.overall)
  }
}
