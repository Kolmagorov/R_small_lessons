# Transforms an amino acid Sequence string into a pretty table

splitSeq<-function(x, nseg=10, ncol=5, space="", filler="-"){
  
  n<-nchar(x)
  ind <- 1:n
  
  if (any(c(nseg, ncol) <= 0)){
    stop("nseg and ncol must be greater than zero\n")
  } 
  
  # Compute padding number
  k<-ceiling(n/nseg/ncol) * nseg * ncol - n
  cat("\nSequence Length:", n,"AA\n")
  cat("\nPadding:", k,"AA\n")
  
  # Padding
  if(k != 0){
    pad<-paste0(rep(filler, k), collapse="")
    x<-paste0(x, pad, collapse="")
  }
  
  segmentSEQ <- substring(text = x
                          , first = seq(1, nchar(x)-1, nseg)
                          , last = seq(nseg, nchar(x),nseg)
                          )
  # Controls spacer character between letters in a segment
  for(j in 1:length(segmentSEQ)){

    temp<-unlist(strsplit(segmentSEQ[j], split=""))
    segmentSEQ[j] <- paste0(temp, collapse = space)
    
  }
  # Building an output matrix with AA sequence distributed
  out<-matrix(segmentSEQ, ncol=ncol, byrow=T)
  rnames<-NULL
  
  
  # Constructing Row names , i.e. ranges
  for (i in 1:nrow(out)){
    rn<-range((1:(n+k))[1:(nseg*ncol) + (i-1)*(nseg*ncol)])
    print(rn)
    # Adjust the last row index subtracting k-padding
    if(i == nrow(out)) {rn[2] <- rn[2] - k}
    rnames <- cbind(rnames,paste("(",rn[1],"-",rn[2],")",sep =""))
  }
  row.names(out)<-rnames
  out
}
