
# define fisher
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

test2 <- function(){
  ##################################################
  # var = 0.1
  data <- simulateData(n_Sites=20000,replicate=3,min_expression=1,max_expression=4,
                       dif_express=1/2,per_me=1/2,dif_me=1,lib_s=1/5,var=0.1)
  res <- DMEseq(data[[1]],data[[2]],data[[3]],data[[4]])
  pvalue <- fisher(data[[1]],data[[2]],data[[3]],data[[4]])
  
  
  library(pROC)
  pdf("rocJM0.1.pdf",height=3.5,width=5)
  label=c(rep(0,10000),rep(1,10000))
  roc1_1 <- plot.roc(label,res[,1], percent=TRUE, col="red",print.auc=TRUE)
  roc2_1 <- lines.roc(label,pvalue, percent=TRUE,col="black",print.auc=TRUE)
  auc1_1 <-auc(roc(label,res[,1]))
  auc2_1 <-auc(roc(label,pvalue))
  dev.off()
  
  ##################################################
  # var = 0.2
  data <- simulateData(n_Sites=20000,replicate=3,min_expression=1,max_expression=4,
                       dif_express=1/2,per_me=1/2,dif_me=1,lib_s=1/5,var=0.2)
  res <- DMEseq(data[[1]],data[[2]],data[[3]],data[[4]])
  pvalue <- fisher(data[[1]],data[[2]],data[[3]],data[[4]])
  
  
  library(pROC)
  pdf("rocJM0.2.pdf",height=3.5,width=5)
  label=c(rep(0,10000),rep(1,10000))
  roc1_2 <- plot.roc(label,res[,1], percent=TRUE, col="red",print.auc=TRUE)
  roc2_2 <- lines.roc(label,pvalue, percent=TRUE,col="black",print.auc=TRUE)
  auc1_2 <-auc(roc(label,res[,1]))
  auc2_2 <-auc(roc(label,pvalue))
  dev.off()
  
  ##################################################
  # var = 0.4
  data <- simulateData(n_Sites=20000,replicate=3,min_expression=1,max_expression=4,
                       dif_express=1/2,per_me=1/2,dif_me=1,lib_s=1/5,var=0.3)
  res <- DMEseq(data[[1]],data[[2]],data[[3]],data[[4]])
  pvalue <- fisher(data[[1]],data[[2]],data[[3]],data[[4]])
  
  
  library(pROC)
  pdf("rocJM0.3.pdf",height=3.5,width=5)
  label=c(rep(0,10000),rep(1,10000))
  roc1_3 <- plot.roc(label,res[,1], percent=TRUE, col="red",print.auc=TRUE)
  roc2_3 <- lines.roc(label,pvalue, percent=TRUE,col="black",print.auc=TRUE)
  auc1_3 <-auc(roc(label,res[,1]))
  auc2_3 <-auc(roc(label,pvalue))
  dev.off()
  
  ##################################################
  # var = 0.8
  data <- simulateData(n_Sites=20000,replicate=3,min_expression=1,max_expression=4,
                       dif_express=1/2,per_me=1/2,dif_me=1,lib_s=1/5,var=0.4)
  res <- DMEseq(data[[1]],data[[2]],data[[3]],data[[4]])
  pvalue <- fisher(data[[1]],data[[2]],data[[3]],data[[4]])
  
  
  library(pROC)
  pdf("rocJM0.4.pdf",height=3.5,width=5)
  label=c(rep(0,10000),rep(1,10000))
  roc1_4 <- plot.roc(label,res[,1], percent=TRUE, col="red",print.auc=TRUE)
  roc2_4 <- lines.roc(label,pvalue, percent=TRUE,col="black",print.auc=TRUE)
  auc1_4 <-auc(roc(label,res[,1]))
  auc2_4 <-auc(roc(label,pvalue))
  dev.off()
  
  
  
  ##################################################
  # var = 1.6
  data <- simulateData(n_Sites=20000,replicate=3,min_expression=1,max_expression=4,
                       dif_express=1/2,per_me=1/2,dif_me=1,lib_s=1/5,var=0.5)
  res <- DMEseq(data[[1]],data[[2]],data[[3]],data[[4]])
  pvalue <- fisher(data[[1]],data[[2]],data[[3]],data[[4]])
  
  
  library(pROC)
  pdf("rocJM0.5.pdf",height=3.5,width=5)
  label=c(rep(0,10000),rep(1,10000))
  roc1_5 <- plot.roc(label,res[,1], percent=TRUE, col="red",print.auc=TRUE)
  roc2_5 <- lines.roc(label,pvalue, percent=TRUE,col="black",print.auc=TRUE)
  auc1_5 <-auc(roc(label,res[,1]))
  auc2_5 <-auc(roc(label,pvalue))
  dev.off()
  
  ##################################################
  # var = 0.6
  data <- simulateData(n_Sites=20000,replicate=3,min_expression=1,max_expression=4,
                       dif_express=1/2,per_me=1/2,dif_me=1,lib_s=1/5,var=0.6)
  res <- DMEseq(data[[1]],data[[2]],data[[3]],data[[4]])
  pvalue <- fisher(data[[1]],data[[2]],data[[3]],data[[4]])
  
  
  library(pROC)
  pdf("rocJM0.6.pdf",height=3.5,width=5)
  label=c(rep(0,10000),rep(1,10000))
  roc1_6 <- plot.roc(label,res[,1], percent=TRUE, col="red",print.auc=TRUE)
  roc2_6 <- lines.roc(label,pvalue, percent=TRUE,col="black",print.auc=TRUE)
  auc1_6 <-auc(roc(label,res[,1]))
  auc2_6 <-auc(roc(label,pvalue))
  dev.off()
  
  ##################################################
  
  dr <- c(auc1_1,auc1_2,auc1_3,auc1_4,auc1_5,auc1_6)
  fi <- c(auc2_1,auc2_2,auc2_3,auc2_4,auc2_5,auc2_6)
  dis <- c(0.1,0.2,0.3,0.4,0.5,0.6)
  res1 <- data.frame(dispersion=dis,auc=dr,method="DRME")
  res2 <- data.frame(dispersion=dis,auc=fi,method="Fisher")
  res <- rbind(res1,res2)
  return(res)
}

# duplicate the results
re1 <- test2()
re2 <- test2()
re3 <- test2()
re4 <- test2()
re5 <- test2()
re6 <- test2()
re7 <- test2()
re8 <- test2()
re9 <- test2()
re10 <- test2()
res <- rbind(re1,re2,re3,re4,re5,re6,re7,re8,re9,re10)


# put result together
p <- ggplot(res, aes(factor(dispersion), auc))
p + + geom_jitter() + geom_boxplot(aes(fill = factor(method)))







