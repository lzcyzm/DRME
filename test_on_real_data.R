library(DESeq)
library(ggplot2)
load("MeRIPseqOriginal_DAA.RData")

# re=removeSites(meth1,meth2,unmeth1,unmeth2,threshold=30)
# meth1=re[[1]]
# meth2=re[[2]]
# unmeth1=re[[3]]
# unmeth2=re[[4]]
res <- DMEseq(meth1,meth2,unmeth1,unmeth2)
#res <- DMEseq(t_meth1,t_meth2,c_unmeth1,c_unmeth2)

# plot ########################################################################################################

hist(res[,1])
hist(p.adjust(res[,1],method="BH"))
temp1 <- (p.adjust(res[,1],"BH") < 0.05); table(sign(res[temp1,3]))
write.table(res,file="output.xls",sep="\t",col.names=TRUE,row.names =FALSE)

# 
pdf("mu_fc_p.pdf",height=3.5,width=5)
p <- ggplot(res, aes(mu, fc))
p +  scale_x_log10() +  ylim(c(-1,1) ) + geom_point(aes(colour = (pvalue>0.05),alpha=1/100 ),size=1)
# p + geom_point(aes(colour = (pval4>0.10),alpha=0.001 ),size=0.001)+scale_x_log10() +  ylim(c(-1,1) ) 
dev.off()

# additional dispersion plot
pdf("mfc_ufc_p.pdf",height=3.5,width=5)
p <- ggplot(res, aes(meth.fc, unmeth.fc))
p +   ylim(c(-4,4) ) + geom_point(aes(colour = (pvalue>0.05),alpha=1/100 ),size=1)+   
  xlim(c(-4,4) )+ geom_abline(intercept = 0)
# p + geom_point(aes(colour = (pval4>0.10),alpha=0.001 ),size=0.001)+scale_x_log10() +  ylim(c(-1,1) ) 
dev.off()

pdf("pvalue.pdf",height=7,width=5)
#p=pvalue=res[,1]
p=res[,1]
par( mfrow = c( 2, 1 ) )
hist(p, breaks=100, col="skyblue", border="slateblue", main="")
hist(p.adjust(p,method="BH"), breaks=100, col="skyblue", border="slateblue", main="")
dev.off()

pdf("condi_nb.pdf",height=3.5,width=10)
p1=dnbinom(0:50, mu=20, size=2)*dnbinom(50-(0:50), mu=200, size=4); p1=p1/sum(p1);
p2=dnbinom(0:50, mu=50, size=4)*dnbinom(50-(0:50), mu=20, size=400);p2=p2/sum(p2);
p3=dnbinom(0:50, mu=20, size=20)*dnbinom(50-(0:50), mu=50, size=1600);p3=p3/sum(p3);
p4=dnbinom(0:50, mu=50, size=5)*dnbinom(50-(0:50), mu=50, size=5);p4=p4/sum(p4);
p5=dnbinom(0:50, mu=50, size=50)*dnbinom(50-(0:50), mu=50, size=50);p5=p5/sum(p5);
p6=dnbinom(0:50, mu=50, size=500)*dnbinom(50-(0:50), mu=50, size=500);p6=p6/sum(p6);
a1<-data.frame(t=0:50,p=p1,label="nb(20,2)*nb(200,4)")
a2<-data.frame(t=0:50,p=p2, label="nb(50,4)*nb(20,400)")
a3<-data.frame(t=0:50,p=p3, label="nb(20,20)*nb(50,1600)")
a4<-data.frame(t=0:50,p=p4, label="nb(50,5)*nb(50,5)")
a5<-data.frame(t=0:50,p=p5, label="nb(50,50)*nb(50,50)")
a6<-data.frame(t=0:50,p=p6, label="nb(50,500)*nb(50,500)")
data <- rbind(a1,a2,a3,a4,a5,a6)
p <- ggplot(data, aes(x=t, y=p, group=label))
p + geom_line(aes(colour = label))
dev.off()




#
pdf("distribution.pdf",height=3.5,width=5)
fdr_DRME=p.adjust(res[,1],method="BH")
fdr_fisher=p.adjust(pvalue,method="BH")
den1=data.frame(q0=res[,6],type="all genes")
den2=data.frame(q0=res[,6][fdr_fisher<0.1],type="Fisher")
den3=data.frame(q0=res[,6][fdr_DRME<0.1],type="DRME")
den=rbind(den1,den2,den3)

#p <- ggplot(den, aes(x=q0,group=type))+scale_x_log10()
p <- ggplot(den, aes(x=log10(q0),group=type))
p + geom_density(aes(colour = type))

dev.off()



library(pROC)
pdf("roc.pdf",height=3.5,width=5)
# reads1=rowSums(unmeth1)
# reads2=rowSums(unmeth2)
label=c(rep(0,10000),rep(1,10000))
# row1=which(reads1>20&reads2>20)
# label=label[row1]
roc1 <- plot.roc(label,res[,1], percent=TRUE, col="1",print.auc=TRUE)
#roc1 <- plot.roc(label,res[,1], percent=TRUE, col="1")
roc2 <- lines.roc(label,pvalue, percent=TRUE,col="2",print.auc=TRUE)
auc1<-auc(roc(label,res[,1]))
auc2<-auc(roc(label,pvalue))
dev.off()



