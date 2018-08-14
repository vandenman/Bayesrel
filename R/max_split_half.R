#' this function uses either an exhaustive search to determine the best possible splithalf reliability
#' or uses a 12x12 hadamard matrix for possible starting points for 12 splits,
#' the latter is much faster when the number of items increases e.g. 10 and up
#' ref: Benton 2013
#'

MaxSplitExhaustive <- function(M){
  #data – matrix of items scores (row=candidates,column=items)
  cov1 <- M
  nite <- ncol(M)
  mat1 <- (bin.combs2(nite)+1)/2
  res1 <- 0
  for (jjz in 1:length(mat1[,1])){
    xal <- mat1[jjz,]
    gutt1 <- 4*(t(xal)%*%cov1%*%(1-xal))/sum(cov1)
    resrand <- gutt1
    if (resrand > res1){
      res1 <- resrand
    }
  }
  return(res1)
}

bin.combs2 <- function (p) {
  retval <- matrix(0, nrow = 2^p, ncol = p)
  for (n in 1:p) {
    retval[, n] <- rep(c(rep(-1, (2^p/2^n)), rep(1, (2^p/2^n))),
                       length = 2^p)
  }
  len <- (nrow(retval)/2)
  comb <- retval[1:len, ]
  comb
}

# this is the solution that is supposed to be faster.
# yet when the covariance matrix has lots of negative entries, it doesnt produce a solution
MaxSplitHalfHad12 <- function(M){
  #data – matrix of items scores (row=candidates,column=items)
  #start with odd vs even
  nite <- ncol(M)
  sequence <- 1:nite
  xal <- (sequence%%2)
  res1 <- MaxSplitHalf(M, xal)
  #now try 12 further splits based on 12*12 Hadamard matrix
  had <- hadamard2(11)
  for (iz in 1:12){
    nextra <- max(nite-12,0)
    resrand <- MaxSplitHalf(M, c(had[,iz], rep(0,nextra))[1:nite])
    if (resrand>res1){res1<-resrand}
  }
  return(res1)
}

#Function to find best split half from a given starting split
MaxSplitHalf = function(M, xal){
  #data – matrix of items scores (row=candidates,column=items)
  #xal – vector of 0s and 1s specifying initial split
  nite = ncol(M)
  cov1 = M
  v = diag(cov1)
  yal = 1-xal
  ones = rep(1,nite)
  covxy = t(xal)%*%cov1%*%yal
  #Code to examine all possible swaps
  maxchg1=9
  while(maxchg1>0){
    #Calculate change for swapping items in X and Y;
    #This is equal to 2covxiyj+covxix+covyyj-vx-vy-covxiy-covxyj;
    covxiyj = cov1
    covxix = (cov1%*%xal)%*%t(ones)
    covyyj = ones%*%(yal%*%cov1)
    vx = v%*%t(ones)
    vy = t(vx)
    covxiy = (cov1%*%yal)%*%t(ones)
    covxyj = ones%*%(xal%*%cov1)
    result = 2*covxiyj+covxix+covyyj-vx-vy-covxiy-covxyj
    for (i in 1:nite){
      for (j in 1:nite){
        if (xal[i]==xal[j]){
          result[i,j]=0}}}
    #Add bits for swapping with no other item
    result = cbind(result,as.vector(cov1%*%xal-cov1%*%yal-v)*xal)
    result = rbind(result,c(as.vector(cov1%*%yal-cov1%*%xal-v)*yal,0))
    #find indices of maximum change;
    maxchg=0
    maxchgx=0
    maxchgy=0
    which1=which(result==max(result),arr.ind=TRUE)[1,]
    if (result[which1[1],which1[2]]>0){maxchgx=which1[1]
    maxchgy=which1[2]
    maxchg=result[which1[1],which1[2]]}
    maxchg1 = maxchg
    if (maxchgx>0 & maxchgx<(nite+1)) {xal[maxchgx]=0}
    if (maxchgy>0 & maxchgy<(nite+1)) {xal[maxchgy]=1}
    if (maxchgx>0 & maxchgx<(nite+1)) {yal[maxchgx]=1}
    if (maxchgy>0 & maxchgy<(nite+1)) {yal[maxchgy]=0}
    covxy = t(xal)%*%cov1%*%yal}
  guttman = 4*covxy/sum(cov1)
  # pites = sum(xal)/nite
  # raju = covxy/(sum(cov1)*pites*(1-pites))
  # v1 = t(xal)%*%cov1%*%xal
  # v2 = t(yal)%*%cov1%*%yal
  # feldt = 4*covxy/(sum(cov1)-((v1-v2)/sqrt(sum(cov1)))**2);
  res = as.vector(guttman)
             # raju=as.vector(raju),
             # feldt=as.vector(feldt),
             # xal=xal)
  return(res)
}

