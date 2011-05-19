thresh <-
function(C, alph, eps=1e-10){
C[abs(C)<alph+eps]<-0
return(C)
}

