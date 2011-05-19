col.norm <-
function(X){
return(apply(X, 2, function(x){sqrt(sum(x^2))}))
}

