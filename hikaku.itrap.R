###############################################################################
#
# April, 26 2016 
#
# hikaku.itrap.R
#
# R function for processing label-free comparative proteomics data.
#
# Help: hikaku.itrap.Rd
# 
# Required packages: gplots, RColorBrewer, numbers
#
# If you have any question or suggestions please contact the author:
#
# Author: Tiago A. de Souza & Fernando G. de Almeida
# tiagoantonio@gmail.com
# github.com/tiagoantonio
#
# Licensed under MIT License (license.txt)
#
#################################################################################

hikaku.itrap= function(filename="filename",samplen, nrep, expl=FALSE, pep=0, psm=0, heat.qt=0.1, pnames=TRUE, housek){
 
#####libraries#####
  
  if (!require("gplots")) { #package required for heatmap graph
    install.packages("gplots", dependencies = TRUE) #install gplots if package is not present
    library(gplots) #load library
  }

  if (!require("RColorBrewer")) { #package required for heatmap
    install.packages("RColorBrewer", dependencies = TRUE) #install RColorBrewer if package is not present
    library(RColorBrewer) #load library
  }

  if (!require("numbers")) { #package required for prime calculations
    install.packages("numbers", dependencies = TRUE) #install numbers if package is not present
    library(RColorBrewer) #load library
  }
  
####acessory functions######
  
  PrimeSumRec=function(x){ # grid for graphs
    if (x==2 | x==1){ # if n combinations = 1 or 2
      x=c(2,1)
      return(x)
    }
    if(isPrime(x)==TRUE){ # if prime sum 1 then factorize
      x=x+1
      x=primeFactors(x)
    }
    else {
      x=primeFactors(x) # if not prime factorize
    }
    if (length(x)==2){ # if prime factors number =2 return vector
      return(x)
    }
    if (length(x)==3){ # if prime factors number equal =3
      x1=length(x)-1
      x.sum=x[length(x)]*x[x1]
      x=x[-length(x)]
      x=x[-(x1)]
      x=append(x, x.sum, after = length(x))
    }
    else { # multiply borders of prime factors
      i=2
      x1=i+1
      x.sum=x[i]*x[i+1]
      x=x[-i]
      x=x[-(x1)]
      x=append(x, x.sum)
      x1=length(x)-1
      x.sum=x[length(x)]*x[x1]
      x=x[-length(x)]
      x=x[-(x1)]
      x=append(x, x.sum, after = length(x))
    }
    return(x) # return vector
  }
  
  count.zero=function(x){ #function to count the sum of replicates without 0s of conditions/samples
    sum(x!=0)
  }

#### parameter check ######
  
  if (samplen==1| samplen>=16){
    stop("Please use samplen (number of conditions) > 1 and < 16")
  }
#######csv#########
  
  report=read.csv(filename, sep="\t", dec = ".", header=T, as.is=T) # reading main input table tab-delimited csv with dec=,
  
####parsing areas, nsample and nrep #######
  
  area.samples=report[8:(7+(samplen*nrep))] #separating sample area data
  rep.matrix=matrix(NA, ncol = samplen, nrow=length(area.samples[,1]), byrow = F) #matrix NA for sample area data
  
  i.sample=seq(from=1,to=(samplen*nrep), by=nrep) #counter for number of samples
  

###explo analysis##  
  
  if (expl==TRUE){ #expl=TRUE for exploratory analysis
    
    x11() #first window for pep, psm, score and hist of pep abundance
    par(mfrow=c(2,2)) #4x4 graph window
    
    pep.column=seq(from=(10+(samplen*nrep)), to=(6+(samplen*nrep))+((samplen*nrep)*4), by=4) #indexes of peptide columns for each sample
    boxplot(report[pep.column], outline=F, col=brewer.pal(name="Set3", n=12), main="Peptides") #boxplot for peptides
    
    psm.column=seq(from=(11+(samplen*nrep)), to=(7+(samplen*nrep))+((samplen*nrep)*4), by=4) #indexes of PSM columns for each sample
    boxplot(report[psm.column], outline=F, col=brewer.pal(name="Set3", n=12), main="PSM") #boxplot for psm
    
    score.column=seq(from=(8+(samplen*nrep)), to=(5+(samplen*nrep))+((samplen*nrep)*4), by=4) #indexes of score columns for each sample
    boxplot(report[score.column], outline=F, col=brewer.pal(name="Set3", n=12), main="Score") #boxplot for scores
    
    hist(report[5][report[5]!=0],xlab="", main="Unique Peptides") #histogram for unique peptides
    
    par(mfrow=c(1,1)) #reseting par
    
    # number of proteins which appears in 0, 1, 2 or 3 replicas
    
    ii=0 # ii for number of conditions
    for (i in i.sample){ #loop for absence of protein in replicates
      ii=ii+1
      rep.sample=apply(area.samples[i:(i+2)],1,count.zero)#1.2omiting NAs, NA=0 if NA mean with others
      rep.matrix[,ii]=rep.sample #matrix containing reps from samples
    }
    
    x11() # second graph for proteins which appears in replicas
    # checking number of samples to build panel
    
    if (samplen<=4){
      par(mfrow=c(2,2))
    }
    if (samplen<=6 & samplen>4){
      par(mfrow=c(2,3))
    }
    if (samplen<=9 & samplen>6){
      par(mfrow=c(3,3))
    }
    if(samplen>9) {
      par(mfrow=c(5,5))
    }
    
    for (i in 1:samplen){ #hist graph for each condition 
      histname=paste("Sample",i)
      sum0=sum(rep.matrix[,i]==0) # sum of presence of each protein in conditions
      sum1=sum(rep.matrix[,i]==1)
      sum2=sum(rep.matrix[,i]==2)
      sum3=sum(rep.matrix[,i]==3)
      
      hist(rep.matrix[,i], breaks = 6,xlab="", main=histname, xaxt='n')
      
      axis(1,at=c(0.25,0.75,1.75,2.75), labels=c(0,1,2,3))
      text(x = 0.25,y=30, sum0) # text in histograms
      text(x = 0.75,y=30, sum1)
      text(x = 1.75,y=30, sum2)
      text(x = 2.75,y=30, sum3)
    }
    par(mfrow=c(1,1))
    return() # explo analysis ends here with no output.
  }
  
  
  ###Comparative Proteomics Analysis#######
  else{
    
    #receiving parameters for peptide and psm filtering
    
    #peptide
    pep.column=seq(from=(10+(samplen*nrep)), to=(6+(samplen*nrep))+((samplen*nrep)*4), by=4)
    ii=7 #start column
    for(i in pep.column){
      ii=ii+1
      pep1=report[,i]
      pep1[pep1<pep]=NA
      report[,i]=pep1
      report[,ii][is.na(report[,i])]=0
    }
    
    #psm
    psm.column=seq(from=(11+(samplen*nrep)), to=(7+(samplen*nrep))+((samplen*nrep)*4), by=4)
    ii=7 # start column
    for(i in psm.column){
      ii=ii+1
      psm1=report[,i]
      psm1[psm1<psm]=NA
      report[,i]=psm1
      report[,ii][is.na(report[,i])]=0
    }
    
    ### Checking NAs in area.samples ####
    
    #areas=0 -> NA
    area.samples=report[8:(7+(samplen*nrep))]
    area.samples[area.samples==0]=NA
    
    ## matrix for mean of conditions
    mean.matrix=matrix(NA, ncol = samplen, nrow=length(area.samples[,1]), byrow = F)
    
    if (pnames==TRUE){
      rownames(mean.matrix)=report$Description
    }
    else{rownames(mean.matrix)=report$Accession}
    
    
    #matrix for sd of conditions #NOT USED
    sd.matrix=matrix(NA, ncol = samplen, nrow=length(area.samples[,1]), byrow = F)
    if (pnames==TRUE){
      rownames(sd.matrix)=report$Description
    }
    else{rownames(sd.matrix)=report$Accession}
    
    #counter conditions
    i.sample=seq(from=1,to=(samplen*nrep), by=nrep)
    
    
    #mean apply for each sample in nrep
    ii=0
    for (i in i.sample){
      ii=ii+1
      mean.sample=apply(area.samples[i:(i+2)],1,mean, na.rm=TRUE)##1.2omiting NAs, NA=0 if NA mean with others
      mean.matrix[,ii]=mean.sample
    }
    
    
    #sd aplly for each sample in nrep #NOT USED
    ii=0
    for (i in i.sample){
      ii=ii+1
      sd.sample=apply(area.samples[i:(i+2)],1,sd, na.rm=TRUE)##1.2omiting NAs, NA=0 if NA mean with others
      sd.matrix[,ii]=sd.sample
    }
    
    #comb 2x2 of conditions
    combinat=combn(samplen,2)
    ncomb=(factorial(samplen)/factorial(samplen-2))/2
    # grid for combinations
    gridn=PrimeSumRec(ncomb)

    #graph1
    x11() # graph of unsorted log2mean
    par(mfrow=c(gridn[2],gridn[1]))
    for (i in 1:ncomb){
      ratio=log2(mean.matrix[,combinat[2,i]]/mean.matrix[,combinat[1,i]])
      plotname=paste(combinat[2,i],":", combinat[1,i] )
      ratio=matrix(ratio)
      median.ratio=apply(ratio,2,median,na.rm=T)
      ratio=ratio-median.ratio # median norm
      ratio=sort(ratio, decreasing=FALSE)
      plot(ratio, main=plotname, col=ifelse(ratio>0, "#66C2A5", ifelse(ratio<0, "#FC8D62", "black")), ylim=c(-10,6))
      abline(h=0)
    }
    par(mfrow=c(1,1))
    
   # Protein names or Accession numbers?
     ratio.matrix=matrix(NA, ncol = ncomb, nrow=length(area.samples[,1]), byrow = F)
    if (pnames==TRUE){
      rownames(ratio.matrix)=report$Description
    }
    else{rownames(ratio.matrix)=report$Accession}
    
    #graph2
    x11() # graph of sorted log2mean
    par(mfrow=c(gridn[2],gridn[1]))
    for (i in 1:ncomb){
      ratio=log2(mean.matrix[,combinat[2,i]]/mean.matrix[,combinat[1,i]])
      plotname=paste(combinat[2,i],":", combinat[1,i] )
      ratio=matrix(ratio)
      median.ratio=apply(ratio,2,median,na.rm=T)
      ratio=ratio-median.ratio # median norm
      ratio.matrix[,i]=ratio
      plot(ratio, main=plotname, col=ifelse(ratio>0, "#66C2A5", ifelse(ratio<0, "#FC8D62", "black")),ylim=c(-8,8))
    }
    par(mfrow=c(1,1))
    
    #graph3
    x11()#heatmap
    #housekeeping normalization & pnames check
    if (pnames==FALSE & hasArg(housek)){
      housek=paste(housek)
      hknorm=mean.matrix[housek,]
      mean.matrix=sweep(mean.matrix,2,hknorm, "/")
    }
    if (pnames==TRUE & hasArg(housek)){
      stop("Please use pnames=FALSE to consider housekeeping norm")
    }
    else{}
    
    # normalization and quantile shown in heatmap
    mean.matrix.heat=mean.matrix
    mean.matrix.heat[is.nan(mean.matrix.heat)]=0
    percqt=heat.qt*quantile(rowSums(mean.matrix.heat))[5] #perc of normalized mean of last quantile
    
    mean.matrix.heat=mean.matrix.heat[rowSums(mean.matrix.heat)>percqt,] #perc of normalized mean of last quantile to not appear in heatmap
    
    heatmap.2(mean.matrix.heat,
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(12,15),     # widens margins around plot
              col=brewer.pal(9,"YlOrRd"),       # use on color palette defined earlier
              Colv="NA")
  }
  if (pnames==FALSE & hasArg(housek)){
	return(mean.matrix) #return a matrix with mean of each protein of each sample = NaNs means Na in all reps.
  }
  else{
    return(ratio.matrix) # return ratio if housekeeping normalization was used
  } 
}