get.thr <-
function(C, n, p, max.num=1, alpha=NULL, step=NULL){
corr<-abs(C[upper.tri(C)])
sc<-sort(corr, decreasing=FALSE)
thr.seq<-NULL
len<-p*(p-1)/2
if(is.null(alpha)) alpha<-1/sqrt(p)
if(is.null(step)) step<-max(1, floor(len/p/n))

for(k in 1:max.num){
D<-rmvnorm(n, sigma=diag(1, p))
D<-t(t(D)/col.norm(D))
D<-t(D)%*%D
ref<-abs(D[upper.tri(D)])
i<-1
while(i<=len){
c<-sc[i]
prob<-sum(ref>c)/len
if(prob<=(len-i+1)/len*alpha) break
i<-i+step
}
thr.seq<-c(thr.seq, c)
}
thr<-median(thr.seq)
list(thr.seq=thr.seq, thr=thr)
}

