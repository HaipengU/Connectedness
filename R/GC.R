GCfunc <- function(Kmatrix, Xmatrix, sigma2a, sigma2e, MUScenario, statistic, NumofMU){
  if (is.factor(MUScenario) != TRUE) stop("Management unit is not factor!")
  Zi <- diag(x = 1, nrow = nrow(Kmatrix),ncol = nrow(Kmatrix))
  Zitr <- t(Zi)
  Kinv <- solve(Kmatrix)
  lamda <- sigma2e/sigma2a
  X <- Xmatrix
  Mabs <- diag(nrow(X)) - X %*% solve(crossprod(X)) %*% t(X)
  CuuK = solve(Zitr %*% Mabs %*% Zi + Kinv*lamda)
  PEVK = CuuK*sigma2e
  Management_Unit <- unique(MUScenario)
  if (statistic == 'PEVD'){
    PEVD.contrast <- matrix(NA, ncol = length(Management_Unit), nrow = length(Management_Unit))
    colnames(PEVD.contrast)<- rownames(PEVD.contrast)<- Management_Unit
    for (i in 1:(length(Management_Unit)-1)) {
      for(j in (i+1):length(Management_Unit)) {
        xcontrast <- matrix(0, ncol = 1, nrow = ncol(Kmatrix))
        indexi <- which(MUScenario == Management_Unit[i])
        indexj <- which(MUScenario == Management_Unit[j])
        xcontrast[indexi] <-  1/length(which(MUScenario == Management_Unit[i]))
        xcontrast[indexj] <- -1/length(which(MUScenario == Management_Unit[j]))
        PEVD.contrast[i,j]<- PEVD.contrast[j,i] <- (crossprod(xcontrast,PEVK) %*% xcontrast)/sigma2a
      }
    }
    PEVD.across.contrast.overall <- mean(PEVD.contrast[upper.tri(PEVD.contrast)])
    if (NumofMU == 'Pairwise'){
      return(PEVD.contrast)
    } else if (NumofMU == 'Overall'){
      return(PEVD.across.contrast.overall)
    }
  } else if (statistic == 'CD'){
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
    CD.across.contrast.overall <- mean(CD.contrast[upper.tri(CD.contrast)])
    if (NumofMU == 'Pairwise'){
      return(CD.contrast)
    } else if (NumofMU == 'Overall'){
      return(CD.across.contrast.overall)
    }
  } else if (statistic == 'r'){
    rijK<- cov2cor(PEVK)
    rij.contrast <- matrix(NA,ncol = length(Management_Unit),nrow = length(Management_Unit))
    colnames(rij.contrast) <- rownames(rij.contrast) <- Management_Unit
    for (i in 1:(length(Management_Unit)-1)) {
      for(j in (i+1):length(Management_Unit)) {
        indexi <- which(MUScenario == Management_Unit[i])
        indexj <- which(MUScenario == Management_Unit[j])
        rij.Aacross <- rijK[indexi,indexj]
        rij.contrast[i,j] <- rij.contrast[j,i] <- mean(rij.Aacross, na.rm = TRUE)
      }
    }
    rij.across.contrast.overall <- mean(rij.contrast[upper.tri(rij.contrast)]) 
    if (NumofMU == 'Pairwise'){
      return(rij.contrast)
    } else if (NumofMU == 'Overall'){
      return(rij.across.contrast.overall)
    }
  }
}
