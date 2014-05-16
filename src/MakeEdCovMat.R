

EditCovMat <- function(filename){
    npops<-length(scan(filename, skip=20, nlines=1))
    np.text <- paste("-n", npops*3, sep="")
    outname <- paste(filename, "ED", sep="")
    system(paste("tail", np.text, filename, ">", outname))
}

