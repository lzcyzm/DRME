
# Generate simualted data
simulateData <- function(n_Sites=20000,replicate=3,min_expression=1,max_expression=3,
                         dif_express=1/4,per_me=1/3,dif_me=1/4,lib_s=1/10,var=2){
  #expression
  q <- 10^runif(n_Sites,min_expression,max_expression)
  p <- runif(n_Sites,0,1)
  
  #difference
  #   x <- rnorm(n_Sites/2,0,dif_me)
  #   dif <- exp(x)
  #   p[round(n_Sites*(1-per_me))+1:n_Sites]=p[round(n_Sites*(1-per_me))+1:n_Sites]*dif
  
  #e
  y <- rnorm(n_Sites,0,dif_express)
  e1 <- exp(y)
  y <- rnorm(n_Sites,0,dif_express)
  e2 <- exp(y)
  T <- rnorm(4*replicate,0,lib_s)
  S <- 2^T
  
  #s
  s_t1 <- S[1:replicate]
  s_t2 <- S[replicate+1:2*replicate]
  s_c1 <- S[2*replicate+1:3*replicate]
  s_c2 <- S[3*replicate+1:4*replicate]
  
  #mu
  mu_t1 <- (e1*q*p)%*%t(as.numeric(s_t1))
  mu_t2 <- (e2*q*p)%*%t(as.numeric(s_t2))
  mu_c1 <- (e1*q*(1-p))%*%t(as.numeric(s_c1))
  mu_c2 <- (e2*q*(1-p))%*%t(as.numeric(s_c2))
  
  #different mu
  x <- rnorm(n_Sites/2,0,dif_me)
  dif <- exp(x)
  mu_t1[round(n_Sites*(1-per_me))+1:n_Sites] <- mu_t1[round(n_Sites*(1-per_me))+1:n_Sites]/dif
  mu_t2[round(n_Sites*(1-per_me))+1:n_Sites] <- mu_t2[round(n_Sites*(1-per_me))+1:n_Sites]*dif
  
  #variance
  # get all variance
  var_t1 <- mu_t1+var*p*q
  var_t2 <- mu_t2+var*p*q
  var_c1 <- mu_c1+var*p*q
  var_c2 <- mu_c2+var*p*q
  
  size_t1=mu_t1^2/(var_t1-mu_t1)
  size_t2=mu_t2^2/(var_t2-mu_t2)
  size_c1=mu_c1^2/(var_c1-mu_c1)
  size_c2=mu_c2^2/(var_c2-mu_c2)
  
  
  #differential expression
  #   y <- rnorm(n_Sites,0,1/4)
  #   t <- 2^y
  #   mu_c1 <- input_mu/t
  #   mu_c2 <- input_mu*t
  #   
  #   #different methylation
  #   mu_t1 <- vector(mode="numeric",length=n_Sites)
  #   mu_t2 <- vector(mode="numeric",length=n_Sites)
  #   
  #   #no difference
  #   mu_t1[1:round(n_Sites*(1-per_me))] <- mu_c1
  #   mu_t2[1:round(n_Sites*(1-per_me))] <- mu_c2
  #   
  #   #difference
  #   z <- rnorm(n_Sites-round(n_Sites*(1-per_me)),0,dif_me)
  #   d <- 2^z
  #   mu_t1[round(n_Sites*(1-per_me)):n_Sites] <- mu_c1/d
  #   mu_t2[round(n_Sites*(1-per_me)):n_Sites] <- mu_c2*d
  #   
  #   #different library size
  #   
  #   mu_t1 <- mu_t1*S
  #   mu_t2 <- mu_t2*S
  #   mu_c1 <- mu_c1*S
  #   mu_c2 <- mu_c2*S
  #   #variance
  #   var_t1 <- mu_t1*var
  #   var_t2 <- mu_t2*var
  #   var_c1 <- mu_c1*var
  #   var_c2 <- mu_c2*var
  
  meth1=matrix(0,n_Sites,3)
  meth2=matrix(0,n_Sites,3)
  unmeth1=matrix(0,n_Sites,3)
  unmeth2=matrix(0,n_Sites,3)
  for(i in 1:n_Sites){
    meth1[i,]=rnbinom(3,mu=mu_t1[i],size=size_t1)
    meth2[i,]=rnbinom(3,mu=mu_t2[i],size=size_t2)
    unmeth1[i,]=rnbinom(3,mu=mu_c1[i],size=size_c1)
    unmeth2[i,]=rnbinom(3,mu=mu_c2[i],size=size_c2)
  }
  res=list(4)
  res[[1]]=meth1
  res[[2]]=meth2
  res[[3]]=unmeth1
  res[[4]]=unmeth2
  
  return (res)
}

# Fisher's test
fisher <- function(meth1,meth2,unmeth1,unmeth2){
  t1 <- rowSums(meth1)
  t2 <- rowSums(meth2)
  c1 <- rowSums(unmeth1)
  c2 <- rowSums(unmeth2)
  s <- .sizeFactor(cbind(t1,t2,c1,c2), useTotal=FALSE)
  t1 <- round(t1/s[1])
  t2 <- round(t2/s[2])
  c1 <- round(c1/s[3])
  c2 <- round(c2/s[4])
  
  # pvalue
  pvalue <- rep(1,length(t1))
  m <- cbind(t1,t2,c1,c2)
  
  for (i in 1:length(t1)) {
    print(i)
    Convictions <-
      matrix(m[i,],
             nrow = 2,
             dimnames =
               list(c("Dizygotic", "Monozygotic"),
                    c("Convicted", "Not convicted")))
    temp <- fisher.test(Convictions)
    pvalue[i] <- temp[[1]]
  }
  return(pvalue)
}

# DRME
DMEseq <- function(meth1,meth2,unmeth1,unmeth2) {
  # don't improve robustness, we don't add 1 to reads count
  s <- .sizeFactor2(cbind(meth1,meth2,unmeth1,unmeth2))
  s_t1 <- s[1:length(meth1[1,])]
  s_t2 <- s[(length(meth1[1,])+1):(length(cbind(meth1,meth2)[1,]))]
  s_c1 <- s[(length(cbind(meth1,meth2)[1,])+1):(length(cbind(meth1,meth2,unmeth1)[1,]))]
  s_c2 <- s[(length(cbind(meth1,meth2,unmeth1)[1,])+1):(length(cbind(meth1,meth2,unmeth1,unmeth2)[1,]))]
  
  # estimate probability of methylation under a condition
  p1 <- .estimateP(meth1, unmeth1, s_t1, s_c1)
  p2 <- .estimateP(meth2, unmeth2, s_t2, s_c2)
  p0 <- .estimateP(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2))
  
  # estimate the abundance of feature
  q0 <- .estimateQ(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2),p0)
  q1 <- .estimateQ(meth1,unmeth1,s_t1,s_c1,p0)
  q2 <- .estimateQ(meth2,unmeth2,s_t2,s_c2,p0)
  
  # estimate size e
  e1 <- q1/q0
  e2 <- q2/q0
  
  
  # calculate methylation reads count variance on common scale, condition 1
  w_t1 <-.calculateW(meth1,s_t1,e1)
  w_t2 <-.calculateW(meth2,s_t2,e2)
  w_c1 <-.calculateW(unmeth1,s_c1,e1)
  w_c2 <-.calculateW(unmeth2,s_c2,e2)
  
  # locfit 
  fit_t1 <- .locfitW(p1,q0,w_t1)
  fit_t2 <- .locfitW(p2,q0,w_t2)
  fit_c1 <- .locfitW(p1,q0,w_c1)
  fit_c2 <- .locfitW(p2,q0,w_c2)
  .plotDispersion(fit_t1,fit_t2)
  
  # calculate z
  z_t1 <- .calculateZ(q0,p1,s_t1,e1)
  z_t2 <- .calculateZ(q0,p2,s_t2,e2)
  z_c1 <- .calculateZ(q0,(1-p1),s_c1,e1)
  z_c2 <- .calculateZ(q0,(1-p2),s_c2,e2)
  
  # get estimate w
  w_fit_t1 <- .fittedW(p0,q0,fit_t1)
  w_fit_t2 <- .fittedW(p0,q0,fit_t2)
  w_fit_c1 <- .fittedW(p0,q0,fit_c1)
  w_fit_c2 <- .fittedW(p0,q0,fit_c1)
  
  # get estimate of upi
  ups_t1 <- pmax(w_fit_t1 - z_t1, 1e-8)
  ups_t2 <- pmax(w_fit_t2 - z_t2, 1e-8)
  ups_c1 <- pmax(w_fit_c1 - z_c1, 1e-8)
  ups_c2 <- pmax(w_fit_c2 - z_c2, 1e-8)
  
  # get all means
  mu_t1 <- (e1*q0*p0)%*%t(as.numeric(s_t1))
  mu_t2 <- (e2*q0*p0)%*%t(as.numeric(s_t2))
  mu_c1 <- (e1*q0*(1-p0))%*%t(as.numeric(s_c1))
  mu_c2 <- (e2*q0*(1-p0))%*%t(as.numeric(s_c2))
  
  # get all variance
  raw_t1 <- (e1%*%t(s_t1))^2*ups_t1
  raw_t2 <- (e2%*%t(s_t2))^2*ups_t2
  raw_c1 <- (e1%*%t(s_c1))^2*ups_c1
  raw_c2 <- (e2%*%t(s_c2))^2*ups_c2
  
  # put mu together
  mu1_t <- rowSums(mu_t1)
  mu2_t <- rowSums(mu_t2)
  mu1_c <- rowSums(mu_c1)
  mu2_c <- rowSums(mu_c2)
  
  # put size together
  size1_t <- (mu1_t^2)/rowSums(raw_t1)
  size2_t <- (mu2_t^2)/rowSums(raw_t2)
  size1_c <- (mu1_c^2)/rowSums(raw_c1)
  size2_c <- (mu2_c^2)/rowSums(raw_c2)
  
  # observation together
  t1 <- rowSums(meth1)
  t2 <- rowSums(meth2)
  c1 <- rowSums(unmeth1)
  c2 <- rowSums(unmeth2)
  t <- t1 + t2
  n1 <- t1 + c1
  n2 <- t2 + c2
  
  # go to test
  res <- .quadNBtest(t1,t,n1,n2,mu1_t,mu2_t,mu1_c,mu2_c,size1_t,size2_t,size1_c,size2_c)
  
  # add fc
  fc <- log(p1/p2)
  m1 <- rowSums(t(t(meth1)/s_t1))
  m2 <- rowSums(t(t(meth2)/s_t2))
  u1 <- rowSums(t(t(unmeth1)/s_c1))
  u2 <- rowSums(t(t(unmeth2)/s_c2))
  mfc <- log(m1)-log(m2)
  ufc <- log(u1)-log(u2)
  
  res2 <- data.frame(res,fc,mfc,ufc,q0)
  res <- res2[,c(1,3:7)]
  colnames(res) <- c("pvalue","mu","fc","meth.fc","unmeth.fc","q0")
  return(res)}

# Other Sub-functions
.sizeFactor <- function(n, useTotal=FALSE) {
  if (useTotal) {
    temp <- log(colSums(n))
    temp <- temp - mean(temp)
    s_size <- exp(temp)
  } else {
    n <- pmax(n,1e-5)
    log_n <- log(n)
    pseudo <- rowMeans(log_n)
    ratio <- log_n-pseudo
    s_size <- exp(apply(ratio,2,median)) }
  return(s_size)
}
.sizeFactor2 <- function(n) {
  temp <- log(colSums(n))
  temp <- temp - mean(temp)
  s_size <- exp(temp)
  return(s_size)
}
.estimateP <- function(meth, unmeth, size_t, size_c) {
  temp_t <- t(t(meth)/size_t)
  temp_c <- t(t(unmeth)/size_c)
  temp_n <- temp_t+temp_c
  p <- rowSums(temp_t)/rowSums(temp_n)
  p[is.na(p)] <- 0.5
  p[is.infinite(p)] <- 0.5
  return(p)
  # which(is.na(p))
  # which(is.infinite(p))  
}
.estimateQ <- function(meth,unmeth,size_t,size_c,p,useAll=FALSE){
  if (useAll) {
    temp_t <- t(t(meth)/size_t)
    temp_c <- t(t(unmeth)/size_c)
    temp_n <- temp_t+temp_c
    q <- rowSums(temp_n)/length(size_t)
  } else {
    temp_c <- t(t(unmeth)/size_c)
    qc <- rowMeans(temp_c)
    q <- qc/(1-p)
  }
  
  q[is.na(q)] <- 0
  return(q)
  # which(is.na(q))
  # which(is.infinite(q))  
}
.calculateW <- function(meth,size_t,e_t){
  temp_t <- t(t(meth)/size_t)
  q <- temp_t/e_t
  #  q[is.na(q)] <- 0
  resi <- q-rowMeans(q)
  w <- rowSums(resi^2)/(length(size_t)-1)
  w <- pmax(w,1e-8)
  return(w)
}
.locfitW <- function(p,q,w) {
  l <- log(q+1)
  data=data.frame(cbind(p,l,w))
  ID <- which(rowSums(is.na(data))>0)
  #  data <- data[-ID,]
  fit=locfit(w~lp(p,l),data=data,family="gamma")
  return(fit)
}
.calculateZ <- function(q,p,size,e){
  temp <- p*q/length(size)
  
  norow <- length(q)
  nocol <- length(size)
  temp2 <- matrix(1,nrow=norow, ncol=nocol )
  temp3 <- t(t(temp2)/size)/e
  
  z <- rowSums(temp3)*temp
  
  #  z[is.infinite(z)] <- 0
  #  z[is.na(z)] <- 0
  return(z)
}
.fittedW <- function(p,q,fit){ 
  library(locfit)
  l <- log(q+1)
  data=data.frame(cbind(p,l))
  w_fit <- predict(fit,data)
  return(w_fit)
}
.quadNBtest <- function(t1,t,n1,n2,mu1_t,mu2_t,mu1_c,mu2_c,size1_t,size2_t,size1_c,size2_c){
  nrows <- length(t)
  pval4 <- rep(1,nrows)
  pval2 <- rep(1,nrows)
  
  t2 <- t-t1
  for (irow in 1:nrows) {
    print(irow)
    trip <- t[irow]
    
    if (trip<1) {p <- NA} else {
      
      trip_t1 <- 0:t[irow]
      trip_t2 <- t[irow] - trip_t1
      trip_c1 <- n1[irow] - trip_t1
      trip_c2 <- n2[irow] - trip_t2
      
      
      p1 <- dnbinom(x=trip_t1, size=size1_t[irow], mu=mu1_t[irow], log = TRUE)
      p2 <- dnbinom(x=trip_t2, size=size2_t[irow], mu=mu2_t[irow], log = TRUE)
      p3 <- dnbinom(x=trip_c1, size=size1_c[irow], mu=mu1_c[irow], log = TRUE)
      p4 <- dnbinom(x=trip_c2, size=size2_c[irow], mu=mu2_c[irow], log = TRUE)
      
      #p4
      p <- p1+p2+p3+p4
      p <- p-max(p)
      p <- exp(p)/sum(exp(p))
      p <- min(1,2*min(sum(p[1:(t1[irow]+1)]),sum(p[(t1[irow]+1):(t[irow]+1)])))- p[(t1[irow]+1)]/2  
      pval4[irow] <- p
      
      #p2
      p <- p1+p2
      p <- p-max(p)
      p <- exp(p)/sum(exp(p))
      p <- min(1,2*min(sum(p[1:(t1[irow]+1)]),sum(p[(t1[irow]+1):(t[irow]+1)])))- p[(t1[irow]+1)]/2  
      pval2[irow] <- p
    }
    
  }
  
  mu <- (mu1_t+mu2_t+mu1_c+mu2_c)/4
  
  res <- data.frame(pval2,pval4,mu)
  return(res)
}
.plotDispersion <-function(fit_t1,fit_t2) {
  p <- rep(seq(from=0,to=1,by=0.01),10)
  q <- rep(seq(from=0,to=10,by=0.1),10)
  p <- matrix(p,nrow = 101, ncol = 101,byrow=TRUE)
  q <- matrix(q,nrow = 101, ncol = 101,byrow=FALSE)
  
  p <- p[1:10201]
  q <- q[1:10201]
  
  w1 <- .fittedW(p,q,fit_t1)
  w1 <- matrix(log(w1),nrow = 101, ncol = 101,byrow=FALSE)
  w2 <- .fittedW(p,q,fit_t2)
  w2 <- matrix(log(w2),nrow = 101, ncol = 101,byrow=FALSE)
  
  pdf("dispersion.pdf",height=4,width=7)
  par(mfrow=c(1,2))
  contour(seq(from=0,to=1,by=0.01),seq(from=0,to=10,by=0.1),w1)
  contour(seq(from=0,to=1,by=0.01),seq(from=0,to=10,by=0.1),w2)
  dev.off()
}

