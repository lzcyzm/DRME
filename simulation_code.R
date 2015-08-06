library(DESeq)
library(ggplot2)

test2 <- function(){
  ##################################################
  # var = 0.1
  data <- simulateData(n_Sites=2000,replicate=3,min_expression=1,max_expression=4,
                       dif_express=1/2,per_me=1/2,dif_me=1,lib_s=1/5,var=0.1)
  res <- DMEseq(data[[1]],data[[2]],data[[3]],data[[4]])
  pvalue <- fisher(data[[1]],data[[2]],data[[3]],data[[4]])
  
  
  library(pROC)
  pdf("rocJM0.1.pdf",height=3.5,width=5)
  label=c(rep(0,1000),rep(1,1000))
  roc1_1 <- plot.roc(label,res[,1], percent=TRUE, col="red",print.auc=TRUE)
  roc2_1 <- lines.roc(label,pvalue, percent=TRUE,col="black",print.auc=TRUE)
  auc1_1 <-auc(roc(label,res[,1]))
  auc2_1 <-auc(roc(label,pvalue))
  dev.off()
  
  ##################################################
  # var = 0.2
  data <- simulateData(n_Sites=2000,replicate=3,min_expression=1,max_expression=4,
                       dif_express=1/2,per_me=1/2,dif_me=1,lib_s=1/5,var=0.2)
  res <- DMEseq(data[[1]],data[[2]],data[[3]],data[[4]])
  pvalue <- fisher(data[[1]],data[[2]],data[[3]],data[[4]])
  
  
  library(pROC)
  pdf("rocJM0.2.pdf",height=3.5,width=5)
  label=c(rep(0,1000),rep(1,1000))
  roc1_2 <- plot.roc(label,res[,1], percent=TRUE, col="red",print.auc=TRUE)
  roc2_2 <- lines.roc(label,pvalue, percent=TRUE,col="black",print.auc=TRUE)
  auc1_2 <-auc(roc(label,res[,1]))
  auc2_2 <-auc(roc(label,pvalue))
  dev.off()
  
  ##################################################
  # var = 0.4
  data <- simulateData(n_Sites=2000,replicate=3,min_expression=1,max_expression=4,
                       dif_express=1/2,per_me=1/2,dif_me=1,lib_s=1/5,var=0.4)
  res <- DMEseq(data[[1]],data[[2]],data[[3]],data[[4]])
  pvalue <- fisher(data[[1]],data[[2]],data[[3]],data[[4]])
  
  
  library(pROC)
  pdf("rocJM0.3.pdf",height=3.5,width=5)
  label=c(rep(0,1000),rep(1,1000))
  roc1_3 <- plot.roc(label,res[,1], percent=TRUE, col="red",print.auc=TRUE)
  roc2_3 <- lines.roc(label,pvalue, percent=TRUE,col="black",print.auc=TRUE)
  auc1_3 <-auc(roc(label,res[,1]))
  auc2_3 <-auc(roc(label,pvalue))
  dev.off()
  
  ##################################################
  
  dr <- c(auc1_1,auc1_2,auc1_3)
  fi <- c(auc2_1,auc2_2,auc2_3)
  dis <- c(0.1,0.2,0.4)
  res1 <- data.frame(dispersion=dis,auc=dr,method="DRME")
  res2 <- data.frame(dispersion=dis,auc=fi,method="Fisher")
  res <- rbind(res1,res2)
  return(res)
}

# duplicate the results
re1 <- test2()
for (i in 1:99) {
  temp <- test2()
  re1 <- rbind(re1,temp)
}

# put result together
pdf(file="AUC.pdf",height = 3,width=5)
p <- ggplot(re1, aes(factor(dispersion), auc))
p + geom_boxplot(aes(fill = factor(method)))
dev.off()