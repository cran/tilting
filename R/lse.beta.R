lse.beta <-
function(X, y, active=NULL){
bhat<-rep(0, ncol(X))
if(!is.null(active)){
Xa<-X[,active, drop=FALSE]
bhat[active]<-drop(solve(t(Xa)%*%Xa)%*%t(Xa)%*%y)
}
return(bhat)
}

