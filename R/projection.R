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

