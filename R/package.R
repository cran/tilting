#library(mvtnorm)
#package.skeleton(name="tilting")
tilting <-
function(X, y, thr.step=NULL, thr.rep=1, max.size=NULL, max.count=NULL, op=2, bic.gamma=1, eps=1e-10){

n<-nrow(X); p<-ncol(X)
if(is.null(max.size)) max.size<-floor(n/2)
if(is.null(max.count)) max.count<-floor(n/2)
if(!is.null(thr.step)) step<-thr.step

active<-NULL; inactive<-setdiff(1:p, active)
Z<-X; res<-y
term<-0; count<-0
thr.seq<-Score<-score.seq<-bic.seq<-NULL

while(term==0){
count<-count+1
if(count>max.count-1) term<-1

temp<-rep(0, p)
colnorm<-col.norm(Z)
Zstar<-t(t(Z)/colnorm)
C<-t(Zstar)%*%Zstar
if(is.null(thr.step)) step<-max(1, floor((p-length(active))*(p-1-length(active))/(n-length(active))/10))
thr<-get.thr(C, n, p-length(active), alpha=min(.05, 1/sqrt(p)), step=step, max.num=thr.rep)$thr
thr.seq<-c(thr.seq, thr)
C1<-thresh(C, thr, eps)

corr<-t(Zstar)%*%res
k<-which.max(abs(corr))
rsv<-setdiff(1:ncol(Z), k)[abs(C1[-k,k])>0]
if(length(rsv)==0){
temp[inactive]<-corr
active<-c(active, inactive[k])
inactive<-setdiff(1:p, active)
score.seq<-c(score.seq, max(abs(corr)))
} else{
J<-c(k, rsv)
score<-NULL
m<-max(2, min(max.size, round(max(apply(C1[J,], 1, function(x){sum(x!=0)-1})))))
for(i in J){
z<-Zstar[, i, drop=FALSE]
rsv<-setdiff(1:ncol(Z), i)[abs(C1[-i,i])>0]
num<-length(rsv)
if(num>m){
if(i!=k){
rsv<-c(k, (rsv[rsv!=k])[sort(abs(C[rsv[rsv!=k],i]), decreasing=TRUE, index=TRUE)$ix[1:(m-1)]])
}else{
rsv<-rsv[sort(abs(C[rsv,i]), decreasing=TRUE, index=TRUE)$ix[1:m]]
}
}
rsv.len<-length(rsv)
while((rsv.len>0) && (min(eigen(C[c(rsv, i), c(rsv, i)])$values)<eps)){ 
rsv<-rsv[-rsv.len]
rsv.len<-rsv.len-1
}
if(rsv.len==0){
score<-c(score, 0)
next
}

PI<-projection(Zstar, rsv)
a<-drop(t(PI%*%z)%*%z)
if(op==1){
score<-c(score, t(z-PI%*%z)%*%res/(1-a))
} else{
ay<-drop(t(PI%*%res)%*%res/(t(res)%*%res))
score<-c(score, t(z-PI%*%z)%*%res/sqrt((1-a)*(1-ay)))
}
}
temp[inactive[J]]<-score
Score<-cbind(Score, temp)
if(sum(score!=0)==0) break
sc<-sort(abs(score), decreasing=TRUE, index=TRUE)
k<-inactive[J[sc$ix[1]]]
active<-c(active, k)
inactive<-setdiff(1:p, active)
score.seq<-c(score.seq, sc$x[1])
}
Score<-cbind(Score, temp)
values<-eigen(t(X[,active, drop=FALSE])%*%X[,active, drop=FALSE])$values
if(abs(max(values)/min(values))>1/eps) break

Z<-X[,-active, drop=FALSE]-projection(X, active)%*%X[,-active, drop=FALSE]
res<-y-projection(X, active)%*%y
bic<-log(sum(res^2)/n)+length(active)/n*(log(n)+2*bic.gamma*log(p))
bic.seq<-c(bic.seq, rep(bic, length(k)))
}

active.hat<-select.model(bic.seq=bic.seq, active=active)

list(Score=Score, active=active, thr.seq=thr.seq, score.seq=score.seq, bic.seq=bic.seq, active.hat=active.hat)

}

select.model <-
function(bic.seq, active){

temp<-which.min(bic.seq)
return(active[1:min(temp)])

}

projection <-
function(X, active=NULL){
if(!is.null(active)){
Xa<-X[,active, drop=FALSE]
PI<-Xa%*%solve(t(Xa)%*%Xa)%*%t(Xa)
} else{
PI<-matrix(0, nrow(X), nrow(X)) 
}
return(PI)
}

lse.beta <-
function(X, y, active=NULL){
bhat<-rep(0, ncol(X))
if(!is.null(active)){
Xa<-X[,active, drop=FALSE]
bhat[active]<-drop(solve(t(Xa)%*%Xa)%*%t(Xa)%*%y)
}
return(bhat)
}

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

thresh <-
function(C, alph, eps=1e-10){
C[abs(C)<alph+eps]<-0
return(C)
}

col.norm <-
function(X){
return(apply(X, 2, function(x){sqrt(sum(x^2))}))
}
