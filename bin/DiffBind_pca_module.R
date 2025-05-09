pv.pcmask <- function(pv,numSites, mask, sites,removeComps,cor=F,bLog=T){
   
   if(missing(numSites)) numSites <- nrow(pv$binding)
   if(is.null(numSites)) numSites <- nrow(pv$binding)  
   numSites <- min(numSites,nrow(pv$binding))
   
   if(missing(sites)) sites <- 1:numSites
   if(is.null(sites)) sites <- 1:numSites
   
   if(missing(mask)) mask <- rep(T,ncol(pv$class))
   for(i in which(mask)) {
      if(nrow(pv$peaks[[i]])==0) {
         mask[i]=F
      }
   }
   if(sum(mask)<2) {
      stop('Need at least two samples for PCA.')
   }
   
   res <- NULL   
   res$class <- pv$class
   pv$values <- pv$binding[sites,c(F,F,F,mask)]
   active   <- apply(pv$values,1,pv.activefun)
   numSites <- min(numSites,sum(active))
   
   pv$values <- pv$values[active,][1:numSites,]
   
   if(!missing(removeComps)) {
      pv$values <- pv.removeComp(pv$values,numRemove=removeComps)
   }
   
   if(bLog) {
      if(max(pv$values)>1) {
         pv$values[pv$values<=1] <- 1
         pv$values <- log2(pv$values)
      }
   }
   
   if(nrow(pv$values) >= sum(mask)) {
      res$pc <- prcomp(pv$values) #,cor=cor)
   }
   res$mask <- mask
   
   return(res)
}

pv.activefun <- function(x){
   if(sum(x>0)>0){
      return(TRUE)
   } else {
      return(FALSE)
   }	
}

PCAbasic <- function(pv, outfile) {
  numSites <- nrow(pv$binding)
  mask <- rep(T,ncol(pv$class))
  pv <- pv.pcmask(pv, numSites, mask, sites=NULL, cor=F, bLog=T)

  pc <- pv$pc

  vr <- rep(0,length(pc$sdev))
  for (i in 1:length(vr)) {
      vr[i] <- pc$sdev[i] ^ 2
  }

  comps=1:3
  c1 <- comps[1]
  c2 <- comps[2]
  c3 <- comps[3]

  c1p <- vr[c1] / sum(vr) * 100
  c2p <- vr[c2] / sum(vr) * 100
  c3p <- vr[c3] / sum(vr) * 100

  plotData <- as.data.frame(pc$rotation[,c(c1,c2,c3)])
  colnames(plotData) <- c("PC1","PC2","PC3")
  write.table(plotData, outfile,quote=F,append=T,sep="\t")
  write.table(paste0("PC1 %: ", c1p), outfile,quote=F,row.names=F,col.names=F,append=T)
  write.table(paste0("PC2 %: ", c2p), outfile,quote=F,row.names=F,col.names=F,append=T)
  write.table(paste0("PC3 %: ", c3p), outfile,quote=F,row.names=F,col.names=F,append=T)  
}
