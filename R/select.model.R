select.model <-
function(bic.seq, active){

temp<-which.min(bic.seq)
return(active[1:min(temp)])

}

