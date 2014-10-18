CCM_boot_par<-function (A, B, E, tau = 1,
                        DesiredL = ((tau * (E - 1) + (E + 1)):length(A) - E + 2),
                        iterations = 100, nclus=NA)
{
  require(doParallel)
  
  if(is.na(nclus)) {
    nclus<-detectCores()
  }
  
  cl <- makeCluster(nclus, timeout=1200)
  registerDoParallel(cl)
  
  out_tmp <- NULL
  plengtht = length(A[!is.na(A)])
  pLibLength = length(A[!is.na(A)])
  Aest = rep(0, length(A))
  rho = Aest
  varrho = Aest
  if (plengtht > pLibLength) {
    plengtht = pLibLength
  }
  gapdist <- tau * (E - 1)
  acceptablelib <- as.numeric(is.finite(A))
  lA <- length(A)
  for (i in 1:gapdist) {
    acceptablelib <- acceptablelib * as.numeric(is.finite(c(rep(NA, 
                                                                i), A[-c((lA - i + 1):lA)])))
  }
  acceptablelib <- which(acceptablelib > 0) - 1
  acceptablelib <- acceptablelib[acceptablelib < ((plengtht - 
                                                     1) - (tau))]
  lengthacceptablelib <- length(acceptablelib)
  DesiredL <- DesiredL + E - 2
  for (i in 1:length(DesiredL)) {
    DesiredL[i] <- acceptablelib[which.min(abs(acceptablelib - 
                                                 DesiredL[i]))]
  }
  DesiredL <- unique(DesiredL)
  A[is.na(A)] <- 0
  B[is.na(B)] <- 0
  if (tau * (E + 1) > lengthacceptablelib) {
    print(paste("Error - too few records to test E =", E, 
                "and tau =", tau))
    return(out = list(A = A, Aest = NA, B = B, rho = NA, 
                      varrho = NA, sdevrho = NA, Lobs = NA, E = out$E, 
                      tau = tau, FULLinfo = NA, rejectedL = NA))
  }
  else {
    parout<-foreach (c(1:iterations)) %dopar% {
      require(multispatialCCM)
      out <- .C("CCM_bootstrap", A = as.double(A), Aest = as.double(Aest), 
                B = as.double(B), rho = as.double(rho), E = as.integer(E), 
                tau = as.integer(tau), plength = as.integer(plengtht), 
                pLibLength = as.integer(pLibLength), DesiredL = as.integer(DesiredL), 
                plengthDesL = as.integer(length(DesiredL)), acceptablelib =
                as.integer(acceptablelib), 
                plengthacceptablelib = as.integer(lengthacceptablelib))
      out$Aest[out$Aest == 0] <- NA
      out_tmpT <- list(Aest = out$Aest, rho = out$rho[out$rho != 
             0], Lobs = (1:length(A))[out$rho != 0] - E + 1)
      lposT = sort(unique(out_tmpT$Lobs))
      datreturn<-list(out_tmpT=out_tmpT, lposT=lposT)
      datreturn
    }
    
    #Unpack
    lpos<-NULL
    out_tmp<-NULL
    for(itcnt in 1:length(parout)) {
      lpos<-sort(unique(c(lpos, parout[[itcnt]]$lposT)))
      out_tmp[[itcnt]]<-parout[[itcnt]]$out_tmpT
    }
    
    
    Aest_mat <- matrix(nrow = length(A), ncol = iterations)
    rho_mat <- matrix(nrow = length(lpos), ncol = iterations)
    rhosq_mat <- matrix(nrow = length(lpos), ncol = iterations)
    for (itlst in 1:iterations) {
      Aest_mat[, itlst] <- out_tmp[[itlst]]$Aest
      rho_mat[, itlst] <- out_tmp[[itlst]]$rho[match(lpos, 
                                                     out_tmp[[itlst]]$Lobs)]
    }
    return(list(A = A, Aest = rowMeans(Aest_mat, na.rm = T), 
                B = B, rho = rowMeans(rho_mat, na.rm = T), sdevrho = apply(rho_mat, 
           1, function(x) sd(x, na.rm = T)), Lobs = lpos, 
                E = E, tau = tau, FULLinfo = rho_mat))
  }
  
  stopCluster(cl)
  rm(cl)
}
  
  
  