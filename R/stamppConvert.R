###########################################################
#
# Import genotype data and convert it into StAMPP format
# 
# Luke Pembleton
# luke.pembleton@dpi.vic.gov.au
#
###########################################################

stamppConvert <-
function (genotype.file, nclusters=1, type="csv"){
    
    if(type=="csv" | type=="r"){ #import genotype data from csv file or R workspace
    
      if(type=="csv"){
        geno <- read.csv(genotype.file)  #import from csv file
      }else{
        geno <- genotype.file #import from R workspace
      }  
      
      library(parallel)
      library(doParallel)
      library(foreach)
      
      cl <- makeCluster(nclusters)
      registerDoParallel(cl) #establish clusters for multithreaded computing 
      
      totalind <- nrow(geno) #number of individuals
      nloc <- ncol(geno)-4 #number of loci/markers
      
      pops <- unique(geno[,2]) #population names
      npops <- length(pops) #number of populations
      
      pop.num <- vector(length=totalind) #create vector of population ID numbers
      
      for (i in 1:totalind){
        pop.num[i]=which(geno[i,2]==pops) #assign population ID numbers to individuals
      }
      
      format <- geno[,4] #genotype format
      
      ploidy <- geno[,3] #ploidy levels
      
      geno <- cbind(geno[,1:2], pop.num, ploidy, format, geno[,5:(4+nloc)]) #combine genotype data with labels to form stampp geno file
      
      comb.geno <- matrix(NA, ncol=nloc, nrow=length(geno[,1])) 
      comb.geno <- cbind(geno[,(1:5)], comb.geno) #matrix to store allele frequency formated genotypes
      
      colnames(comb.geno) <- colnames(geno)
      
      ab.geno <- subset(geno, geno[,5]=="BiA") #subset individuals with AB coded genotypes
      nind.ab.geno <- length(ab.geno[,2])
      
      freq.geno <- subset(geno, geno[,5]=="freq") #subset individuals with genotypes stored as allele frequencies
      nind.freq.geno <- length(freq.geno[,2])
      
      res <- foreach(i = 1:nloc, .combine=cbind, .inorder=TRUE) %dopar% {  
        
        #convert AB format genotypes to allele frequency format genotypes
        
        notmval <- ab.geno[,(5+i)]!=-9
        
        a <- gsub("B", "", ab.geno[,(5+i)])
        a <- nchar(a)*notmval
        b <- gsub("A", "", ab.geno[,(5+i)])
        b <- nchar(b)*notmval
        (a/(a+b))
        
      } 
      
      colnames(res) <- colnames(freq.geno[,6:(5+nloc)])
      
      comb.geno[,(6:(5+nloc))]=apply(rbind(res, freq.geno[,6:(5+nloc)]), 2, as.numeric) #store allele frequency genotypes in matrix
      
      comb.geno[comb.geno==-9] <- NA #replace missing marker genotypes with NaN to be ignored in future calculations
      
      stopCluster(cl)
      
      return(comb.geno)
    
    }
    
    if(type=="genlight"){
      
      library(adegenet)
            
      geno2 <- genotype.file
      
      geno <- as.matrix(geno2) #extract genotype data from genlight object
            
      sample <- row.names(geno) #individual names
      pop.names <- pop(geno2) #population names
      ploidy <- ploidy(geno2) #ploidy level
      geno=geno*(1/ploidy) #convert genotype data (number of allele B) to precentage allele frequency
      geno[is.na(geno)]=NaN 
      format <- vector(length=length(geno[,1])) 
      format[1:length(geno[,1])]="genlight"
          
      
      pops <- unique(pop.names) #population names
            
      pop.num <- vector(length=length(geno[,1])) #create vector of population ID numbers
      
      for (i in 1:length(geno[,1])){
        pop.num[i]=which(pop.names[i]==pops) #assign population ID numbers to individuals
      }
      
      genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, ploidy, format))
      
      geno <- cbind(genoLHS, geno) #combine genotype data with labels to form stampp geno file
      
      geno[,2]=as.character(pop.names)
      geno[,4]=geno2@ploidy
        
      row.names(geno)=NULL
      
      
      return(geno)
      
      
    }
    
    
  }
