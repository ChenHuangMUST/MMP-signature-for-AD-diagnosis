load("D:/RProject/AD/AD_train.RData")
load("D:/RProject/AD/NMF/afterbatch.marker.select.3.11")
dat_nondense=dat_rmbatch
dat_nondense[dat_nondense<0] <- -1
dat_nondense[dat_nondense>=0] <- 1
index=order(pd[,2])
pd=pd[index,]

index=intersect(pd[,1],colnames(dat_nondense))
pd=pd[index,]
table(pd[,2])
i=90
out=matrix(nrow = 3486,ncol = 5)
colnames(out)=c("MPP","X-squared","df","p-value","p.adj")
for (i in 1:3486) {
  a=table(dat_nondense[i,1:487])
  b=table(dat_nondense[i,488:975])
  
  a=if (length(a) == 2) {
    a
  } else {
    c(a, 487 - a)
  }
  b=if (length(b) == 2) {
    b
  } else {
    c(b, 488 - b)
  }
  a
  b
  m=chisq.test(rbind(a,b))
  out[i,]=c(rownames(dat_nondense)[i],m$statistic,m$parameter,m$p.value,NA)
  
}
out[,5]=p.adjust(out[,4])
imp_MPP=out[which(as.numeric(out[,4])<0.05),]
save(imp_MPP,file="imp_MPP.RData")
save(dat_rmbatch,file="MPP.RData")
save.image()
load(".RData")