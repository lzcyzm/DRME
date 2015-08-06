# Fisher's exact test
# reads1=rowSums(unmeth1)
# reads2=rowSums(unmeth2)
# re_meth1=meth1[reads1>30&reads2>30,]
# re_meth2=meth2[reads1>30&reads2>30,]
# re_unmeth1=unmeth1[reads1>30&reads2>30,]
# re_unmeth2=unmeth2[reads1>30&reads2>30,]
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



pf <- pvalue
pd <- res[,1]
p <- data.frame(pf,pd)
lg.p <- -log10(p)

library(ggplot2)
pdf("p_comp.pdf",height=3.5,width=5)
p <- ggplot(lg.p, aes(pf, pd))
p +  geom_point(aes(alpha=1/100),size=1)+ xlim(c(0,8)) +ylim(c(0,8))+geom_abline(intercept = 0, colour = "red", size=0.05)
dev.off()

# another test
lg.p[is.na(lg.p)] <- 0
pthresh = seq(from=100,to=1,by=-1)
count <- matrix(rep(0,300),nrow=100)
colnames(count) <- c("fisher","drme","both")
for (i in 1:100) {
  th <- pthresh[i]
  count[i,1:2] <- colSums(lg.p > th)
  count[i,3] <- sum(rowSums(lg.p > th)==2)
}
dd <- data.frame(count,pthresh)

d1 <- data.frame(count=dd[,1],threshold=dd[,4],type="Fisher")
d2 <- data.frame(count=dd[,2],threshold=dd[,4],type="DRME")
d3 <- data.frame(count=dd[,3],threshold=dd[,4],type="Overlap")
d <- rbind(d1,d2,d3)


pdf("overlap2.pdf",height=3.5,width=5)
p <- ggplot(d, aes(x=threshold, y=count, group=type))+scale_x_log10()+ xlim(c(0,10))
p + geom_line(aes(colour = type))
dev.off()




