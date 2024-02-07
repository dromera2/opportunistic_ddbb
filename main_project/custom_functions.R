# The communitySimmul function has been obtained from https://github.com/guiblanchet/HMSC/blob/master/R/communitySimul.R

#' @title Simulate a species community
#'
#' @description Simulates a species community assuming a linear model with the proposed descriptors, the species traits and one or more random effects.
#'
#' @param X Matrix of descriptors defining the species. Rows are sites and the columns are descriptors.
#' @param Tr Matrix defining species traits. Rows are traits and columns are species.
#' @param Phylo Square correlation matrix defining the phylogenetic relationship among pairs of species.
#' @param Random A factor or a data.frame that includes many factors defining a random effect on the sites.
#' @param Auto A data.frame that includes a factor as its first column and spatial/temporal coordinates for the other columns. If multiple autocorrelated random effects are considered a list is given where each element of the list is a data.frame as defined previously.
#' @param nsp Numeric. Number of species to be simulated in the community.
#' @param family A character string defining the error distribution and link function to be used when simulating the data. Only "gaussian","probit", "poisson" and "nbinomial" have been implemented so far.
#' @param paramDist Numeric. A value defining the standard deviation (if \code{family="gaussian"}) or the dispersal (if  \code{family="nbinomial"}) matrix of model parameters defining how each species (rows) is characterized by each descriptors (columns).
#' @param paramX A matrix of model parameters defining how each species (rows) is characterized by each descriptors (columns).
#' @param paramTr A matrix of model parameters defining how each descriptors (rows) characterizes species traits (columns).
#' @param paramPhylo A numeric value defining the importance of phylogeny in structuring species relationship with the environment.
#' @param meansParamX Vector of means used to generate the model parameters. (See details)
#' @param varX Covariance matrix used to generate the model parameters. (See details)
#' @param latent A list of matrices where each of set defines a random effect on the sites (rows) characterized by different latent variables (columns).
#' @param paramLatent A list of matrices where each of set defines the parameters of a random effect on the sites (rows) characterized by different latent variables (columns).
#' @param spCor A species by species correlation matrix defining the correlations among species in the community.
#' @param shrinkLocal A list of matrices defining local shrinkage parameters for each parameter of the latent variables associated to each random effect.
#' @param paramShrinkGlobal A list of vectors defining the independent global shrinkage parameters for each latent variable.
#' @param latentAuto A list of matrices where each of set defines an autocorrelated random effect on the sites (rows) characterized by different autocorrelated latent variables (columns).
#' @param paramLatentAuto A list of matrices where each of set defines the parameters of an autocorrelated random effect on the sites (rows) characterized by different autocorrelated latent variables (columns).
#' @param paramAuto A list of numerical values defining the importance of the autocorrelation for each autocorrelated latent variables \code{latentAuto}. These values can range from 0 to the largest distance between samples in the autocorrelated level. (See details)
#' @param shrinkLocalAuto A list of matrices defining local shrinkage parameters for each parameter of the autocorrelated latent variables associated to each autocorrelated random effect.
#' @param paramShrinkGlobalAuto A list of vectors defining the independent global shrinkage parameters for each autocorrelated latent variable.
#'
#' @details
#'
#' This function can be used to simulate a community from randomly proposed or fixed parameters.
#'
#' The descriptors in \code{X} are used without any modifications (or additions) to simulate the species community. As such, a column of 1 should be included in \code{X} for the model used to construct the community to include an intercept.
#'
#' The values in \code{meansParamX} and \code{varX} are used as parameter of a multivariate normal distribution to generate the model's parameters (\code{\link[MASS]{mvrnorm}} is used in the function). When \code{paramX} is set to \code{NULL}, the \code{meansParamX} and the \code{varX} will be randomly generated if they are also set to \code{NULL}. When values are given to \code{paramX} the values of the \code{meansParamX} and the \code{varX} are not used and if either is set to \code{NULL}, no data will be generated for either set of parameter. When generated, the values for \code{meansParamX} are randomly sampled from \code{\link{rnorm}} and the values for \code{varX} are randomly sampled from an inverse (\code{solve})  Wishart distribution (\code{\link{rWishart}}).
#'
#' Note that \code{meansParamX} can be calculated directly by \code{Tr} with \code{paramTr}. As such, \code{meansParamX} is made available as an extension. If it is given to the function but \code{paramTr} is not, than \code{paramTr} is calculated from \code{meansParamX} and vice versa. If both \code{meansParamX} and \code{paramTr} are given to the function and there is a mismatch in the parameters calculated, \code{meansParamX} will take precedent.
#'
#' All the parameters associated to the autocorrelated random effect (\code{Auto},\code{latentAuto},\code{paramLantentAuto} and \code{paramAuto}) use the distance between groups (levels of the factors in \code{Auto}) to define the parameters and latent variables. When multiple samples are within a group, the coordinates related to the different samples are averaged both calculating the distances between this group and other groups.
#'
#' Currently, this function simulates four types of community:
#'
#' \itemize{
#' 	\item{\code{gaussian}}{ For normally distributed data.}
#' 	\item{\code{probit}}{ For presence/absence of species using a probit link.}
#' 	\item{\code{poisson}}{ For count data of species using a log link.}
#' 	\item{\code{nbinomial}}{ For count data with many zeros.}
#' }

#' @return
#'
#' The functions \code{communitySimulH} and \code{communitySimulHT} return an object of class \code{communitySimul} with the following components:
#' \item{data}{an object of class HMSCdata}
#' \item{param}{an object of class HMSCparam}
#' \item{sd}{The standard deviation used to simulate normally distributed data. This value is the same as paramDist and appears only when the data is simulated using \code{gaussian}}
#' \item{size}{The dispersal parameter used when simulating data distributed following a negative binomial distribution. This value is the same as paramDist and appears only when the data is simulated using \code{nbinomial}}
#' \item{probMat}{A matrix that defines the occurrence probability of each species for each site.}
#'
#' @author F. Guillaume Blanchet
#'
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rWishart
#' @importFrom stats rnbinom
#' @importFrom stats rpois
#' @importFrom stats pnorm
#' @examples
#'
#' ### Construct some descriptors
#' #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#' ### Simulate community from random parameters
#' #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#' desc <- cbind(1, scale(1:50), scale(1:50)^2)
#' traits <- rbind(1, 1:20/5)
#' random1 <- as.factor(1:50)
#' random2 <- as.factor(rep(letters[1:2], each = 25))
#' randEff <- data.frame(rand1 = random1, rand2 = random2)
#'
#' ### Simulate presence-absence community data
#' commDescTraitRandEffProbit <- communitySimul(X = desc, Tr = traits, Random = randEff, nsp = 20)
#'
#' #--------------------
#' ### Fixing parameters
#' #--------------------
#' ## ProbMat does not change
#' TrueParamX <- commDescTraitRandEffProbit$param$paramX
#' commParamX <- communitySimul(X = desc, paramX = TrueParamX)
#'
#' ### ProbMat changes
#' TrueParamTr <- commDescTraitRandEffProbit$param$paramTr
#' commParamTr <- communitySimul(X = desc, Tr = traits, paramTr = TrueParamTr)
#'
#' ### ProbMat changes
#' commParamLatent <- communitySimul(X = desc, Random = randEff, paramLatent = commDescTraitRandEffProbit$param$paramLatent, latent = commDescTraitRandEffProbit$param$latent)
#'
#' ### ProbMat changes
#' speciesCor <- cov2cor(solve(rWishart(1, 22, diag(20))[, , 1]))
#' commParamSpCor <- communitySimul(X = desc, spCor = speciesCor)
#'
#' @keywords datagen
#' @keywords htest
#' @importFrom MASS mvrnorm
#' @importFrom stats rgamma
#' @importFrom stats dist
#' @export
communitySimul<-function(X=NULL,Tr=NULL,Phylo=NULL,Random=NULL,Auto=NULL,nsp=NULL,family="probit",paramDist=NULL,paramX=NULL,paramTr=NULL,paramPhylo=NULL,meansParamX=NULL,varX=NULL,latent=NULL,paramLatent=NULL,spCor=NULL,shrinkLocal=NULL,paramShrinkGlobal=NULL,latentAuto=NULL,paramLatentAuto=NULL,paramAuto=NULL,shrinkLocalAuto=NULL,paramShrinkGlobalAuto=NULL){
  #### F. Guillaume Blanchet - July 2014, October 2014, January 2016, May 2016
  ##########################################################################################
  ### This makes the outlierSp inactive
  outlierSp<-NULL
  outlierWeight<-0.001
  
  #================================================
  ### Standard format of the X matrix if X is given
  #================================================
  if(!is.null(X)){
    X<-as.matrix(X)
    nsite<-nrow(X)
    nparamX<-ncol(X)
    
    if(is.null(colnames(X))){
      colnames(X)<-paste("x",1:nparamX,sep="")
    }
    
    if(is.null(rownames(X))){
      rownames(X)<-paste("site",1:nsite,sep="")
    }
    
    #=====================
    ### If paramX is given
    #=====================
    if(!is.null(paramX)){
      if(is.vector(paramX)){
        nsp<-length(paramX)
        paramX<-matrix(paramX,ncol=1)
      }else{
        nsp<-nrow(paramX)
      }
    }
  }
  
  #======================
  ### If paramTr is given
  #======================
  if(!is.null(paramTr)){
    if(is.null(Tr)){
      stop("'Tr' needs to be given if 'paramTr' is given")
    }
    nsp<-ncol(Tr)
  }
  
  #=========================
  ### If paramPhylo is given
  #=========================
  if(!is.null(paramPhylo)){
    if(is.null(Phylo)){
      stop("'Phylo' needs to be given if 'paramPhylo' is given")
    }
    nsp<-ncol(Phylo)
  }
  
  #==========================
  ### If paramLatent is given
  #==========================
  if(!is.null(paramLatent)){
    if(is.null(latent) && is.null(Random)){
      stop("If 'paramlatent' is given, 'latent' and 'Random' also needs to be given")
    }
    nsp<-nrow(paramLatent[[1]])
  }
  
  #=====================
  ### If latent is given
  #=====================
  if(!is.null(latent)){
    if(is.null(paramLatent) && is.null(Random)){
      stop("If 'latent' is given, 'paramLatent' and 'Random' also needs to be given")
    }
  }
  
  #====================
  ### If spCor is given
  #====================
  if(!is.null(spCor)){
    if(!is.null(latent) && !is.null(paramLatent)){
      stop("If 'spCor' is given, 'latent' and 'paramLatent' do not need to be given")
    }
    if(!is.null(Random)){
      stop("If 'spCor' is given, 'Random' does not need to be given")
    }
    nsp<-nrow(spCor)
    
    ### Basic check
    if(min(spCor) < -1 | max(spCor) > 1){
      stop("'spCor' is a correlation matrix, so that values in it should range from -1 to 1")
    }
    
    if(nrow(spCor)!=ncol(spCor)){
      stop("'spCor' needs to be a correlation matrix with as many rows as columns")
    }
    
    if(!is.matrix(spCor)){
      spCor<-as.matrix(spCor)
    }
    
    ### Add row and column names to spCor
    if(is.null(colnames(spCor))){
      colnames(spCor)<-paste("sp",1:nsp,sep="")
    }
    
    if(is.null(rownames(spCor))){
      rownames(spCor)<-paste("sp",1:nsp,sep="")
    }
  }
  
  #==========================
  ### If shrinkLocal is given
  #==========================
  if(!is.null(shrinkLocal)){
    if((is.null(paramLatent) & is.null(latent)) && is.null(Random)){
      stop("If 'shrinkLocal' is given, 'paramLatent', 'latent' and 'Random' also needs to be given")
    }
  }
  
  #================================
  ### If paramShrinkGlobal is given
  #================================
  if(!is.null(paramShrinkGlobal)){
    if((is.null(paramLatent) & is.null(latent)) && is.null(Random)){
      stop("If 'paramShrinkGlobal' is given, 'paramLatent', 'latent' and 'Random' also needs to be given")
    }
  }
  
  #========================
  ### If paramAuto is given
  #========================
  if(!is.null(paramAuto)){
    if(is.null(latentAuto) && is.null(paramLatentAuto) && is.null(Auto)){
      stop("If 'paramAuto' is given 'latentAuto' and 'paramLatentAuto' also needs to be given")
    }
    
    nparamAuto<-lapply(paramAuto,length)
    nparamLatentAuto<-lapply(paramLatentAuto,ncol)
    nLatentAuto<-lapply(latentAuto,ncol)
    
    if(!((length(nparamAuto)==length(nparamLatentAuto))==(length(nLatentAuto)==length(nparamAuto)))){
      stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same length")
    }else{
      if(!(nparamAuto%in%nparamLatentAuto)==(nparamAuto%in%nLatentAuto)){
        stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same values")
      }
    }
  }
  
  #==============================
  ### If paramLatentAuto is given
  #==============================
  if(!is.null(paramLatentAuto)){
    if(is.null(latentAuto) && is.null(paramAuto) && is.null(Auto)){
      stop("If 'paramLatentAuto' is given, 'paramAuto', 'latentAuto' and 'Auto' also needs to be given")
    }
    nsp<-nrow(paramLatentAuto[[1]])
    
    nparamAuto<-lapply(paramAuto,length)
    nparamLatentAuto<-lapply(paramLatentAuto,ncol)
    nLatentAuto<-lapply(latentAuto,ncol)
    
    if(!((length(nparamAuto)==length(nparamLatentAuto))==(length(nLatentAuto)==length(nparamAuto)))){
      stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same length")
    }else{
      if(!(nparamAuto%in%nparamLatentAuto)==(nparamAuto%in%nLatentAuto)){
        stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same values")
      }
    }
  }
  
  #=========================
  ### If latentAuto is given
  #=========================
  if(!is.null(latentAuto)){
    if(is.null(paramLatentAuto) && is.null(paramAuto) && is.null(Auto)){
      stop("If 'latentAuto' is given, 'paramLatentAuto', 'paramAuto' and 'Auto' also needs to be given")
    }
    
    nparamAuto<-lapply(paramAuto,length)
    nparamLatentAuto<-lapply(paramLatentAuto,ncol)
    nLatentAuto<-lapply(latentAuto,ncol)
    
    if(!((length(nparamAuto)==length(nparamLatentAuto))==(length(nLatentAuto)==length(nparamAuto)))){
      stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same length")
    }else{
      if(!(nparamAuto%in%nparamLatentAuto)==(nparamAuto%in%nLatentAuto)){
        stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same values")
      }
    }
  }
  
  #==============================
  ### If shrinkLocalAuto is given
  #==============================
  if(!is.null(shrinkLocalAuto)){
    if((is.null(paramLatentAuto) & is.null(latentAuto)) && is.null(Auto)){
      stop("If 'shrinkLocalAuto' is given, 'paramLatentAuto', 'latentAUto' and 'Auto' also needs to be given")
    }
  }
  
  #====================================
  ### If paramShrinkGlobalAuto is given
  #====================================
  if(!is.null(paramShrinkGlobalAuto)){
    if((is.null(paramLatentAuto) & is.null(latentAuto)) && is.null(Auto)){
      stop("If 'paramShrinkGlobalAuto' is given, 'paramLatentAuto', 'latentAUto' and 'Auto' also needs to be given")
    }
  }
  
  #===========================================
  ### If X, latent and latentAuto are all NULL
  #===========================================
  if(is.null(X) & is.null(latent) & is.null(Random) & is.null(Auto) & is.null(latentAuto)){
    stop("At least one of 'X', 'Random', 'latent', 'Auto' or 'latentAuto' needs to be given")
  }
  
  #===================================
  ### Standard format of the Tr matrix
  #===================================
  if(!is.null(Tr)){
    if(!is.vector(Tr)){
      Tr<-as.matrix(Tr)
    }else{
      Tr<-matrix(Tr,nrow=1)
    }
    ### Define nsp if it is not yet defined
    if(!is.null(nsp)){
      if(ncol(Tr)!=nsp){
        stop("'nsp' should be the same as the number of column of Tr")
      }
    }else{
      nsp<-ncol(Tr)
    }
    
    ### Organize Tr
    nTr<-nrow(Tr)
    
    if(is.null(colnames(Tr))){
      colnames(Tr)<-paste("sp",1:nsp,sep="")
    }
    if(is.null(rownames(Tr))){
      rownames(Tr)<-paste("t",1:nTr,sep="")
    }
  }
  
  #======================================
  ### Standard format of the Phylo matrix
  #======================================
  if(!is.null(Phylo)){
    if(ncol(Phylo)!=nrow(Phylo)){
      stop("'Phylo' should be a square correlation matrix")
    }
    
    ### Check if symmetry
    if(!isSymmetric(Phylo)){
      stop("'Phylo' need to be a symmetric matrix")
    }
    
    ### Check positive definiteness
    if(any(eigen(Phylo)$value<0)){
      stop("'Phylo' need to be a positive definite matrix")
    }
    
    ### Check if Phylo is a covariance matrix and convert to correlation if it is
    if(any(diag(Phylo)!=1)){
      stop("'Phylo' need to be a correlation matrix")
    }
    
    ### Define nsp if it is not yet defined
    if(!is.null(nsp)){
      if(ncol(Phylo)!=nsp){
        stop("'nsp' should be the same as the number of column (and rows) of 'Phylo'")
      }
    }else{
      nsp<-ncol(Phylo)
    }
    
    if(is.null(colnames(Phylo))){
      colnames(Phylo)<-paste("sp",1:nsp,sep="")
    }
    if(is.null(rownames(Phylo))){
      rownames(Phylo)<-paste("sp",1:nsp,sep="")
    }
  }
  
  #=======================================
  ### Standard format of the Random effect
  #=======================================
  if(!is.null(Random)){
    if(is.factor(Random)){
      ### Number of sites
      nsite<-length(Random)
      
      Random<-data.frame(random1=Random)
      rownames(Random)<-paste("site",1:nsite,sep="")
      ### Number or random effects
      nRandom<-1
    }else{
      if(is.data.frame(Random)){
        if(!all(mapply(is.factor,Random))){
          stop("If 'Random' is a data.frame, it should only include factors")
        }
        ### Number or random effects
        nRandom<-ncol(Random)
        ### Number of sites
        nsite<-nrow(Random)
        
        colnames(Random)<-paste("random",1:nRandom,sep="")
        rownames(Random)<-paste("site",1:nsite,sep="")
      }else{
        stop("'Random' should be a factor or a data.frame")
      }
    }
    
    ### Number of levels in each random effect considered
    nLevelRandom<-mapply(nlevels,Random)
  }
  
  #======================================================
  ### Standard format of the autocorrelated Random effect
  #======================================================
  if(!is.null(Auto)){
    if(is.data.frame(Auto)){
      nsite<-nrow(Auto)
      
      if(!is.factor(Auto[,1])){
        stop("When Auto is a data.frame, its first column should be a factor and the other columns should be coordinates")
      }
      if(!all(apply(as.matrix(Auto[,-1]),2,is.numeric))){
        stop("When Auto is a data.frame, its first column should be a factor and the other columns should be coordinates")
      }
      
      ### Name stuff
      rownames(Auto)<-paste("site",1:nsite,sep="")
      colnames(Auto)[1]<-paste("autoRandom",1)
      colnames(Auto)[-1]<-paste("coord",1:(ncol(Auto)-1))
      
      Auto<-list(auto1=Auto)
      ### Number or autocorrelated random effects
      nAuto<-1
    }else{
      if(is.list(Auto)){
        nsite<-nrow(Auto[[1]])
        
        if(!all(mapply(is.data.frame,Auto))){
          stop("If 'Auto' is a list, it should only include data.frames")
        }
        
        if(!all(sapply(Auto,function(x) is.factor(x[,1])))){
          stop("When Auto is a list, the first column of the data.frame in each list should be a factor and the other columns should be coordinates")
        }
        
        ### Find the number of columns with coordinates
        AutoCoord<-lapply(lapply(Auto,function(x) x[,-1]),as.matrix)
        
        if(!all(unlist(lapply(AutoCoord,function(x) apply(x,2,is.numeric))))){
          stop("When Auto is a list, the first column of the data.frame in each list should be a factor and the other columns should be coordinates")
        }
        
        ### Number or autocorrelated random effects
        nAuto<-length(Auto)
        
        ### Name stuff
        for(i in 1:nAuto){
          rownames(Auto[[i]])<-paste("site",1:nsite,sep="")
          colnames(Auto[[i]])[1]<-"autoRandom1"
          colnames(Auto[[i]])[-1]<-paste("coord",1:(ncol(Auto[[i]])-1),sep="")
        }
        
        names(Auto)<-paste("auto",1:nAuto,sep="")
        
      }else{
        stop("'Auto' should be a data.frame or a list")
      }
    }
    
    ### Number of levels in each random effect considered
    nLevelAuto<-sapply(Auto,function(x) nlevels(x[,1]))
  }
  
  #=========================================================
  ### If X, Random, paramX, paramLatent and latent are given
  #=========================================================
  if(!is.null(paramX)){
    if(!is.null(paramLatent)){
      if(!is.null(latent)){
        if(!is.null(X)){
          if(!is.null(Random)){
            for(i in 1:nRandom){
              nLatent <- ncol(latent[[1,i]])
              
              ### Check if the levels match between Random and latent
              matchLev <- which(levels(Random[,i]) %in% rownames(latent[[1,i]]))
              locLev <- which(rownames(latent[[1,i]]) %in% Random[matchLev,i])
              
              ### Build new latent variables
              latentNew<-matrix(NA,nrow=nLevelRandom[i],ncol=nLatent)
              ### If none of the levels match
              if(length(matchLev) == 0){
                latentNew <- matrix(rnorm(nLevelRandom[i]*nLatent),nrow=nLevelRandom[i],ncol=nLatent)
              }
              
              ### If the levels match
              if(length(matchLev) == nLevelRandom[i]){
                latentNew <- latent[[1,i]]
              }
              
              ### If only some of the levels match
              if(length(matchLev) > 0 && length(matchLev) < nLevelRandom[i]){
                latentNew[matchLev,] <- latent[[1,i]][locLev,]
                latentNew[-matchLev,] <- matrix(norm((nLevelRandom[i]-length(matchLev))*nLatent),nrow=nLevelRandom[i],ncol=nLatent)
              }
              ### Assign latentNew to latent
              latent[[1,i]] <- latentNew
            }
          }
        }
      }
    }
  }
  
  #===============================================================
  ### If X, Auto, paramX, paramLatentAuto and latentAuto are given
  #===============================================================
  if(!is.null(paramX)){
    if(!is.null(paramLatentAuto)){
      if(!is.null(latentAuto)){
        if(!is.null(X)){
          if(!is.null(Auto)){
            for(i in 1:nAuto){
              nLatentAuto <- ncol(latentAuto[[1,i]])
              
              ### Check if the levels match between Auto and latentAuto
              matchLev <- which(levels(Auto[[i]][,1]) %in% rownames(latentAuto[[1,i]]))
              locLev <- which(rownames(latentAuto[[1,i]]) %in% Auto[[i]][matchLev,1])
              
              ### Build new latentAuto variables
              latentAutoNew<-matrix(NA,nrow=nLevelAuto[i],ncol=nLatentAuto)
              ### If none of the levels match
              if(length(matchLev) == 0){
                latentAutoNew <- matrix(rnorm(nLevelAuto[i]*nLatentAuto),nrow=nLevelAuto[i],ncol=nLatentAuto)
              }
              
              ### If the levels match
              if(length(matchLev) == nLevelAuto[i]){
                latentAutoNew <- latentAuto[[1,i]][locLev,]
              }
              
              ### If only some of the levels match
              if(length(matchLev) > 0 && length(matchLev) < nLevelAuto[i]){
                latentAutoNew[matchLev,] <- latentAuto[[1,i]][locLev,]
                latentAutoNew[-matchLev,] <- matrix(norm((nLevelAuto[i]-length(matchLev))*nLatentAuto),nrow=nLevelAuto[i],ncol=nLatentAuto)
              }
              ### Assign latentAutoNew to latentAuto
              latentAuto[[1, i]] <- latentAutoNew
            }
          }
        }
      }
    }
  }
  
  #========================================
  ### Set parameters for Tr and meansParamX
  #========================================
  if(!is.null(X)){
    if(!is.null(Tr)){
      if(is.null(paramTr)){
        if(is.null(meansParamX)){
          paramTr<-matrix(rnorm(nTr*nparamX),nrow=nparamX,ncol=nTr)
          rownames(paramTr)<-paste("p",1:nparamX,sep="")
          colnames(paramTr)<-paste("t",1:nTr,sep="")
          
          ### Community Mean (projector of traits and parameters)
          meansParamX<-matrix(NA,ncol=nparamX,nrow=nsp)
          for(i in 1:nsp){
            meansParamX[i,]<-tcrossprod(Tr[,i],paramTr)
          }
          rownames(meansParamX)<-paste("sp",1:nsp,sep="")
          colnames(meansParamX)<-paste("p",1:nparamX,sep="")
        }else{
          paramTrAll<-array(dim=c(nparamX,nTr,nsp))
          for(i in 1:nsp){
            paramTrAll[,,i]<-meansParamX[i,]%*%matrix(Tr[,i],nrow=1) # matrix() is to account for the case where there is only 1 trait
          }
          ### Approximation of paramTr
          paramTr<-apply(paramTrAll,1:2,mean)
        }
      }else{
        if(is.null(meansParamX)){
          ### Community Mean (projector of traits and parameters)
          meansParamX<-matrix(NA,ncol=nparamX,nrow=nsp)
          for(i in 1:nsp){
            meansParamX[i,]<-tcrossprod(Tr[,i],paramTr)
          }
          rownames(meansParamX)<-paste("sp",1:nsp,sep="")
          colnames(meansParamX)<-paste("p",1:nparamX,sep="")
        }else{
          ### Check if meansParamX and meansParamTr (means calculated using Tr and paramTr) match
          meansParamTr<-matrix(NA,ncol=nparamX,nrow=nsp)
          for(i in 1:nsp){
            meansParamTr[i,]<-tcrossprod(Tr[,i],paramTr)
          }
          
          if(abs(sum(meansParamX-meansParamTr)) > 10^(-4)){
            warnings("There might be a mismatch between 'meansParamX' and the 'means' calculated using 'Tr' and 'paramTr'")
            warning("The values in 'meansParamX' will be used for the remaining data simulation")
          }
          
          rownames(meansParamX)<-paste("sp",1:nsp,sep="")
          colnames(meansParamX)<-paste("p",1:nparamX,sep="")
        }
      }
    }else{
      meansParamX<-matrix(rnorm(nparamX),ncol=1)
    }
  }
  
  #=========================
  ### Define outlier species
  #=========================
  outlierWeightVec<-rep(1,nsp)
  if(!is.null(outlierSp)){
    noutlierSp<-round(nsp*outlierSp)
    outlierWeightVec[1:noutlierSp]<-outlierWeight
  }
  
  #==========================
  ### Set parameter for Phylo
  #==========================
  if(!is.null(Phylo)){
    if(is.null(paramPhylo)){
      paramPhylo<-runif(1,0,1)
    }
  }
  
  #=======================
  ### Set parameters for X
  #=======================
  if(!is.null(X)){
    if(is.null(Phylo)){
      ### Mean of the community paramX
      if(is.null(paramX)){
        ### Sigma matrix of the community
        if(is.null(varX)){
          varX<-chol2inv(chol(rWishart(1,nparamX+2,diag(nparamX))[,,1]))
        }else{
          if(!isSymmetric(varX)){
            stop("'varX' is not a symmetric matrix")
          }
        }
        if(is.null(colnames(varX))){
          colnames(varX)<-paste("p",1:nparamX,sep="")
        }
        if(is.null(rownames(varX))){
          rownames(varX)<-paste("p",1:nparamX,sep="")
        }
        
        ### paramX for each species
        paramX<-matrix(NA,nrow=nsp,ncol=nparamX)
        if(!is.null(Tr)){
          for(i in 1:nsp){
            paramX[i,]<-MASS::mvrnorm(1,mu=meansParamX[i,],Sigma=(1/outlierWeightVec[i])*varX)
          }
        }else{
          for(i in 1:nsp){
            paramX[i,]<-MASS::mvrnorm(1,mu=meansParamX,Sigma=(1/outlierWeightVec[i])*varX)
          }
        }
        colnames(paramX)<-paste("p",1:nparamX,sep="")
        rownames(paramX)<-paste("sp",1:nsp,sep="")
      }else{
        if(is.null(Tr)){
          if(is.null(meansParamX)){
            meansParamX<-matrix(colMeans(paramX),ncol=1)
          }
        }
        if(is.null(varX)){
          varX<-crossprod(sweep(paramX,2,colMeans(paramX),FUN="-"))
        }
      }
    }else{
      ### Sigma matrix of the community
      if(is.null(varX)){
        varX<-chol2inv(chol(rWishart(1,nparamX+2,diag(nparamX))[,,1]))
      }else{
        if(!isSymmetric(varX)){
          stop("'varX' is not a symmetric matrix")
        }
      }
      if(is.null(colnames(varX))){
        colnames(varX)<-paste("p",1:nparamX,sep="")
      }
      if(is.null(rownames(varX))){
        rownames(varX)<-paste("p",1:nparamX,sep="")
      }
      
      ### Weighted phylogeny
      if(paramPhylo>=0){
        wPhylo<-paramPhylo * Phylo + (1-paramPhylo)*diag(nsp)
      }else{
        iPhylo<-cov2cor(chol2inv(chol(Phylo)))
        wPhylo<-(-paramPhylo) * iPhylo + (1+paramPhylo)*diag(nsp)
      }
      
      varXPhylo <- kronecker(varX,wPhylo);
      
      if(!is.null(Tr)){
        paramX<-MASS::mvrnorm(1,as.vector(meansParamX),varXPhylo)
        paramX<-matrix(paramX,nrow=nsp)
      }else{
        paramX<-MASS::mvrnorm(1,rep(meansParamX,nsp),varXPhylo)
        paramX<-matrix(paramX,nrow=nsp,byrow=TRUE)
      }
    }
    #===================
    ### Model estimation
    #===================
    EstModel<-tcrossprod(X,paramX)
  }
  
  #==================================================
  ### Set parameters for the latent part of the model
  #==================================================
  if(!is.null(Random)){
    if(is.null(latent)){
      ### Object storing the latent variables and the parameters
      latent<-vector("list",length=nRandom)
      dim(latent)<-c(1,ncol(Random))
      colnames(latent)<-paste("random",1:nRandom,sep="")
    }
    
    if(is.null(paramLatent)){
      paramLatent<-vector("list",length=nRandom)
      dim(paramLatent)<-c(1,ncol(Random))
      names(paramLatent)<-paste("random",1:nRandom,sep="")
    }
    
    if(is.null(shrinkLocal)){
      ### Object storing shrinkLocal
      shrinkLocal<-vector("list",length=nRandom)
      dim(shrinkLocal)<-c(1,ncol(Random))
      colnames(shrinkLocal)<-paste("random",1:nRandom,sep="")
    }else{
      if(length(shrinkLocal)!=nRandom){
        stop("'shrinkLocal' needs to have the same length Random")
      }
    }
    
    if(is.null(paramShrinkGlobal)){
      ### Object storing shrinkLocal
      paramShrinkGlobal<-vector("list",length=nRandom)
      dim(paramShrinkGlobal)<-c(1,ncol(Random))
      colnames(paramShrinkGlobal)<-paste("random",1:nRandom,sep="")
    }else{
      if(length(paramShrinkGlobal)!=nRandom){
        stop("'paramShrinkGlobal' needs to have the same length Random")
      }
    }
    
    latentCov<-array(dim=c(nsp,nsp,nRandom))
    dimnames(latentCov)[[1]]<-paste("sp",1:nsp,sep="")
    dimnames(latentCov)[[2]]<-paste("sp",1:nsp,sep="")
    dimnames(latentCov)[[3]]<-paste("random",1:nRandom,sep="")
    
    ### Construct the model
    for(i in 1:nRandom){
      ### Define the latent variables and their associated parameters
      if(is.null(latent[[1,i]])){
        if(is.null(paramLatent[[1,i]])){
          if(is.null(paramShrinkGlobal[[1,i]])){
            if(is.null(shrinkLocal[[1,i]])){
              ### Define the number of latent variables in the model (it will vary with the number of random effects)
              nLatent<-sample(2:5,nRandom)
            }else{
              nLatent<-sapply(shrinkLocal,ncol)
            }
          }else{
            nLatent<-sapply(paramShrinkGlobal,length)
          }
        }else{
          if(is.null(paramShrinkGlobal[[1,i]])){
            if(is.null(shrinkLocal[[1,i]])){
              ### Define the number of latent variables in the model (it will vary with the number of random effects)
              nLatent<-sapply(paramLatent[[1,i]],ncol)
            }else{
              nLatent<-sapply(shrinkLocal,ncol)
            }
          }else{
            nLatent<-sapply(paramShrinkGlobal,length)
          }
        }
        ### Construct latent variables
        latent[[1,i]]<-matrix(rnorm(nLevelRandom[i]*nLatent[i]),nrow=nLevelRandom[i],ncol=nLatent[i])
        rownames(latent[[1,i]])<-levels(Random[,i])
        colnames(latent[[1,i]])<-paste("latent",1:nLatent[i],sep="")
      }
      
      ### Define paramLatent and their associated parameters
      if(is.null(paramLatent[[1,i]])){
        if(is.null(shrinkLocal[[1,i]])){
          shrinkLocal[[1,i]]<-matrix(rgamma(nsp*nLatent[i],1.5,1.5),nsp,nLatent[i])
        }else{
          nLatentTest<-sapply(shrinkLocal,ncol)
          if(!all(nLatentTest%in%nLatent)){
            stop("'shrinkLocal' needs to have the same number levels as Random")
          }
        }
        
        if(is.null(paramShrinkGlobal[[1,i]])){
          paramShrinkGlobal[[1,i]]<-c(rgamma(1,50,1),rgamma(nLatent[i]-1,50,1))
        }else{
          nLatentTest<-sapply(paramShrinkGlobal,length)
          if(!all(nLatentTest%in%nLatent)){
            stop("'paramShrinkGlobal' needs to have the same number levels as Random")
          }
        }
        
        ### Construct latent variables
        shrinkGlobal<-cumprod(paramShrinkGlobal[[1,i]])
        normSD<-as.vector(1/matrix(rep(shrinkGlobal,nsp),nrow=nsp,byrow=TRUE)*shrinkLocal[[1,i]])
        
        paramLatent[[1,i]]<-matrix(rnorm(length(normSD),0,normSD),nsp,nLatent[i])
        colnames(paramLatent[[1,i]]) <- colnames(latent[[1,i]])
        rownames(paramLatent[[1,i]]) <- paste("sp",1:nsp,sep="")
      }
      if(!is.null(X)){
        ### Add the random effect to the estimated model
        EstModel<-EstModel+tcrossprod(latent[[1,i]][Random[,i],],paramLatent[[1,i]])
      }else{
        ### Estimated model
        EstModel<-tcrossprod(latent[[1,i]][Random[,i],],paramLatent[[1,i]])
      }
      
      ### Calculate the variance-covariance matrix from the latent variables
      latentCov[,,i]<-tcrossprod(paramLatent[[1,i]])
    }
  }
  
  if(!is.null(X)){
    if(!is.null(spCor)){
      EstModel<-EstModel+MASS::mvrnorm(n=nsite,mu=rep(0,nsp),Sigma=spCor)
    }
  }else{
    if(!is.null(spCor)){
      EstModel<-MASS::mvrnorm(n=nsite,mu=rep(0,nsp),Sigma=spCor)
    }
  }
  
  #===========================================================
  ### Set parameters for the autocorrelation part of the model
  #===========================================================
  if(!is.null(Auto)){
    #__________________________________________________________
    ### Calculate the distance between sampled within one level
    #__________________________________________________________
    ### Number of levels in each auto effect considered
    nLevelAuto<-sapply(Auto,function(x) nlevels(x[,1]))
    AutoDist<-vector("list",length=nAuto)
    
    for(i in 1:nAuto){
      nAutoCoord<-ncol(Auto[[i]])-1
      AutoCoordMean<-matrix(NA,nrow=nLevelAuto[i],ncol=nAutoCoord)
      
      for(j in 1:nAutoCoord){
        AutoCoordMean[,j]<-tapply(Auto[[i]][,j+1],Auto[[i]][,1],mean)
      }
      AutoDist[[i]]<-dist(AutoCoordMean)
    }
    
    if(is.null(latentAuto)){
      ### Construct latentAuto (no values yet)
      latentAuto<-vector("list",length=nAuto)
      dim(latentAuto)<-c(1,nAuto)
      colnames(latentAuto)<-paste("auto",1:nAuto,sep="")
    }
    
    if(is.null(paramLatentAuto)){
      paramLatentAuto<-vector("list",length=nAuto)
      dim(paramLatentAuto)<-c(1,nAuto)
      names(paramLatentAuto)<-paste("auto",1:nAuto,sep="")
    }
    
    if(is.null(shrinkLocalAuto)){
      ### Object storing shrinkLocalAuto
      shrinkLocalAuto<-vector("list",length=nAuto)
      dim(shrinkLocalAuto)<-c(1,nAuto)
      colnames(shrinkLocalAuto)<-paste("auto",1:nAuto,sep="")
    }else{
      if(length(shrinkLocalAuto)!=nAuto){
        stop("'shrinkLocalAuto' needs to have the same length Auto")
      }
    }
    
    if(is.null(paramShrinkGlobalAuto)){
      ### Object storing paramShrinkGlobalAuto
      paramShrinkGlobalAuto<-vector("list",length=nAuto)
      dim(paramShrinkGlobalAuto)<-c(1,nAuto)
      colnames(paramShrinkGlobalAuto)<-paste("auto",1:nAuto,sep="")
    }else{
      if(length(paramShrinkGlobalAuto)!=nAuto){
        stop("'paramShrinkGlobalAuto' needs to have the same length Auto")
      }
    }
    
    ### Number of levels in each autocorrelated random effect considered
    nLevelAuto<-sapply(Auto,function(x) nlevels(x[,1]))
    
    ### Construct the model
    for(i in 1:nAuto){
      ### Define the latent variables and their associated parameters
      if(is.null(latentAuto[[1,i]])){
        if(is.null(paramLatentAuto[[1,i]])){
          if(is.null(paramShrinkGlobalAuto[[1,i]])){
            if(is.null(shrinkLocalAuto[[1,i]])){
              ### Define the number of autocorrelated latent variables in the model (it will vary with the number of random effects)
              nLatentAuto<-sample(2:5,nAuto)
            }else{
              nLatentAuto<-sapply(shrinkLocalAuto,ncol)
            }
          }else{
            nLatentAuto<-sapply(paramShrinkGlobalAuto,length)
          }
        }else{
          if(is.null(paramShrinkGlobalAuto[[1,i]])){
            if(is.null(shrinkLocalAuto[[1,i]])){
              ### Define the number of autocorrelated latent variables in the model (it will vary with the number of random effects)
              nLatentAuto<-sapply(paramLatentAuto[[1,i]],ncol)
            }else{
              nLatentAuto<-sapply(shrinkLocalAuto,ncol)
            }
          }else{
            nLatentAuto<-sapply(paramShrinkGlobalAuto,length)
          }
        }
        
        #--------------------------------------
        ### Define paramAuto if it is not given
        #--------------------------------------
        if(is.null(paramAuto)){
          paramAuto<-vector("list",length=nAuto)
          for(i in 1:nAuto){
            paramAuto[[i]]<-runif(nLatentAuto[i])*max(AutoDist[[1]])
          }
        }else{
          if(length(paramAuto)!=nAuto){
            stop("The length of 'paramAuto' should be the same as the length of Auto")
          }
          if(all(sapply(paramAuto,length)!=nLatentAuto)){
            stop("Each part of 'paramAuto' should have the same length as the number of latent variables in each part of Auto")
          }
        }
        #-----------------------------------------------------------
        ### Construct object to weighted the distance with paramAuto
        #-----------------------------------------------------------
        ### Use Exponential function
        wAutoDist<-vector("list",length=nAuto)
        for(i in 1:nAuto){
          wAutoDist[[i]]<-vector("list",length=nLatentAuto[i])
        }
        
        for(i in 1:nAuto){
          AutoDistMat<-as.matrix(AutoDist[[i]])
          for(j in 1:nLatentAuto[i]){
            wAutoDist[[i]][[j]]<-exp(-AutoDistMat/paramAuto[[i]][j])
          }
        }
        
        #__________________________________
        ### Autocorrelated latent variables
        #__________________________________
        for(i in 1:nAuto){
          latentAuto[[1,i]]<-matrix(NA,nrow=nLevelAuto[i],ncol=nLatentAuto[i])
          for(j in 1:nLatentAuto[i]){
            latentAuto[[1,i]][,j]<-rmvnorm(1,rep(0,nrow(wAutoDist[[i]][[j]])),wAutoDist[[i]][[j]])
          }
          rownames(latentAuto[[1,i]])<-levels(Auto[[i]][,1])
          colnames(latentAuto[[1,i]])<-paste("latentAuto",1:ncol(latentAuto[[1,i]]),sep="")
        }
      }
      
      ### Define paramLatentAuto and their associated parameters
      if(is.null(paramLatentAuto[[1,i]])){
        if(is.null(shrinkLocalAuto[[1,i]])){
          shrinkLocalAuto[[1,i]]<-matrix(rgamma(nsp*nLatentAuto[i],1.5,1.5),nsp,nLatentAuto[i])
        }else{
          nLatentAutoTest<-sapply(shrinkLocalAuto,ncol)
          if(!all(nLatentAutoTest%in%nLatentAuto)){
            stop("'shrinkLocalAuto' needs to have the same number levels as Auto")
          }
        }
        
        if(is.null(paramShrinkGlobalAuto[[1,i]])){
          paramShrinkGlobalAuto[[1,i]]<-c(rgamma(1,50,1),rgamma(nLatentAuto[i]-1,50,1))
        }else{
          nLatentAutoTest<-sapply(paramShrinkGlobalAuto,length)
          if(!all(nLatentAutoTest%in%nLatentAuto)){
            stop("'paramShrinkGlobalAuto' needs to have the same number levels as Random")
          }
        }
        
        ### Construct latent variables
        shrinkGlobalAuto<-cumprod(paramShrinkGlobalAuto[[1,i]])
        normSD<-as.vector(1/matrix(rep(shrinkGlobalAuto,nsp),nrow=nsp,byrow=TRUE)*shrinkLocalAuto[[1,i]])
        
        paramLatentAuto[[1,i]]<-matrix(rnorm(length(normSD),0,normSD),nsp,nLatentAuto[i])
        colnames(paramLatentAuto[[1,i]]) <- colnames(latentAuto[[1,i]])
        rownames(paramLatentAuto[[1,i]]) <- paste("sp",1:nsp,sep="")
      }
      
      if(!is.null(Random)){
        if(!is.null(X)){
          ### Add the auto effect to the estimated model
          EstModel<-EstModel+tcrossprod(latentAuto[[1,i]][Auto[[i]][,1],],paramLatentAuto[[1,i]])
        }
      }else{
        if(!is.null(X)){
          ### Add the auto effect to the estimated model
          EstModel<-EstModel+tcrossprod(latentAuto[[1,i]][Auto[[i]][,1],],paramLatentAuto[[1,i]])
        }else{
          ### Estimated model
          EstModel<-tcrossprod(latentAuto[[1,i]][Auto[[i]][,1],],paramLatentAuto[[1,i]])
        }
      }
    }
  }
  
  #======================================
  ### Construct species occurrence matrix
  #======================================
  
  if(family=="gaussian"){
    Y<-matrix(rnorm(nsite*nsp,mean=as.vector(EstModel),sd=paramDist),nrow=nsite,ncol=nsp)
  }
  if(family=="nbinomial"){
    Y<-matrix(rnbinom(nsite*nsp,mu=as.vector(EstModel),size=paramDist),nrow=nsite,ncol=nsp)
  }
  if(family=="probit"){
    Ylatent<-matrix(rnorm(nsite*nsp,mean=as.vector(EstModel),sd=1),nrow=nsite,ncol=nsp)
    Y<-Ylatent
    Y[Y>0]<-1 # * much faster than ifelse()
    Y[Y<0]<-0
  }
  if(family=="poisson"){
    Y<-matrix(rpois(nsite*nsp,lambda=as.vector(exp(EstModel))),nrow=nsite,ncol=nsp)
  }
  
  rownames(Y)<-paste("site",1:nsite,sep="")
  colnames(Y)<-paste("sp",1:nsp,sep="")
  
  #=============================
  ### Return the results objects
  #=============================
  ### Data
  if(is.null(Auto)){
    if(!is.null(X)){
      if(is.null(Phylo)){
        if(!is.null(Random)){
          if(!is.null(Tr)){
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Random=Random)
            attributes(data)<-list(names=c("Y","X","Tr","Random"),Ypattern="sp")
          }else{
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Random=Random)
            attributes(data)<-list(names=c("Y","X","Random"),Ypattern="sp")
          }
        }else{
          if(!is.null(Tr)){
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr))
            attributes(data)<-list(names=c("Y","X","Tr"),Ypattern="sp")
          }else{
            data<-list(Y=as.matrix(Y),X=as.matrix(X))
            attributes(data)<-list(names=c("Y","X"),Ypattern="sp")
          }
        }
      }else{
        if(!is.null(Random)){
          if(!is.null(Tr)){
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Phylo=as.matrix(Phylo),Random=Random)
            attributes(data)<-list(names=c("Y","X","Tr","Phylo","Random"),Ypattern="sp")
          }else{
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Phylo=as.matrix(Phylo),Random=Random)
            attributes(data)<-list(names=c("Y","X","Phylo","Random"),Ypattern="sp")
          }
        }else{
          if(!is.null(Tr)){
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Phylo=as.matrix(Phylo))
            attributes(data)<-list(names=c("Y","X","Tr","Phylo"),Ypattern="sp")
          }else{
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Phylo=as.matrix(Phylo))
            attributes(data)<-list(names=c("Y","X","Phylo"),Ypattern="sp")
          }
        }
      }
    }else{
      if(!is.null(Random)){
        data<-list(Y=as.matrix(Y),Random=Random)
        attributes(data)<-list(names=c("Y","Random"),Ypattern="sp")
      }
    }
  }else{
    if(!is.null(X)){
      if(is.null(Phylo)){
        if(!is.null(Random)){
          if(!is.null(Tr)){
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Random=Random,Auto=Auto)
            attributes(data)<-list(names=c("Y","X","Tr","Random","Auto"),Ypattern="sp")
          }else{
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Random=Random,Auto=Auto)
            attributes(data)<-list(names=c("Y","X","Random","Auto"),Ypattern="sp")
          }
        }else{
          if(!is.null(Tr)){
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Auto=Auto)
            attributes(data)<-list(names=c("Y","X","Tr","Auto"),Ypattern="sp")
          }else{
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Auto=Auto)
            attributes(data)<-list(names=c("Y","X","Auto"),Ypattern="sp")
          }
        }
      }else{
        if(!is.null(Random)){
          if(!is.null(Tr)){
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Phylo=as.matrix(Phylo),Random=Random,Auto=Auto)
            attributes(data)<-list(names=c("Y","X","Tr","Phylo","Random","Auto"),Ypattern="sp")
          }else{
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Phylo=as.matrix(Phylo),Random=Random,Auto=Auto)
            attributes(data)<-list(names=c("Y","X","Phylo","Random","Auto"),Ypattern="sp")
          }
        }else{
          if(!is.null(Tr)){
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Phylo=as.matrix(Phylo),Auto=Auto)
            attributes(data)<-list(names=c("Y","X","Tr","Phylo","Auto"),Ypattern="sp")
          }else{
            data<-list(Y=as.matrix(Y),X=as.matrix(X),Phylo=as.matrix(Phylo),Auto=Auto)
            attributes(data)<-list(names=c("Y","X","Phylo","Auto"),Ypattern="sp")
          }
        }
      }
    }else{
      if(!is.null(Random)){
        data<-list(Y=as.matrix(Y),Random=Random,Auto=Auto)
        attributes(data)<-list(names=c("Y","Random","Auto"),Ypattern="sp")
      }else{
        data<-list(Y=as.matrix(Y),Auto=Auto)
        attributes(data)<-list(names=c("Y","Auto"),Ypattern="sp")
      }
    }
    
  }
  class(data)<-"HMSCdata"
  
  ### Parameters
  if(is.null(Auto)){
    if(!is.null(X)){
      if(is.null(Phylo)){
        if(!is.null(Random)){
          if(!is.null(Tr)){
            allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal)
          }else{
            allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal)
          }
        }else{
          if(!is.null(spCor)){
            if(!is.null(Tr)){
              allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,spCor=spCor)
            }else{
              allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,spCor=spCor)
            }
          }else{
            if(!is.null(Tr)){
              allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec)
            }else{
              allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec)
            }
          }
        }
      }else{
        if(!is.null(Random)){
          if(!is.null(Tr)){
            allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal)
          }else{
            allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal)
          }
        }else{
          if(!is.null(spCor)){
            if(!is.null(Tr)){
              allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,spCor=spCor)
            }else{
              allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,spCor=spCor)
            }
          }else{
            if(!is.null(Tr)){
              allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec)
            }else{
              allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec)
            }
          }
        }
      }
    }else{
      if(!is.null(Random)){
        allparam<-list(outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal)
      }else{
        if(!is.null(spCor)){
          allparam<-list(outlier=outlierWeightVec,spCor=spCor)
        }
      }
    }
  }else{
    if(!is.null(X)){
      if(is.null(Phylo)){
        if(!is.null(Random)){
          if(!is.null(Tr)){
            allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
          }else{
            allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
          }
        }else{
          if(!is.null(spCor)){
            if(!is.null(Tr)){
              allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,spCor=spCor,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
            }else{
              allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,spCor=spCor,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
            }
          }else{
            if(!is.null(Tr)){
              allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
            }else{
              allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
            }
          }
        }
      }else{
        if(!is.null(Random)){
          if(!is.null(Tr)){
            allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
          }else{
            allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
          }
        }else{
          if(!is.null(spCor)){
            if(!is.null(Tr)){
              allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,spCor=spCor,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
            }else{
              allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,spCor=spCor,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
            }
          }else{
            if(!is.null(Tr)){
              allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
            }else{
              allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
            }
          }
        }
      }
    }else{
      if(!is.null(Random)){
        allparam<-list(outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
      }else{
        if(!is.null(spCor)){
          allparam<-list(outlier=outlierWeightVec,spCor=spCor,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
        }else{
          allparam<-list(outlier=outlierWeightVec,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
        }
      }
    }
  }
  
  class(allparam)<-"HMSCparam"
  
  ### Make sure allparam are all matrices
  allparam<-lapply(allparam,as.matrix)
  
  ### Final result object
  if(family=="gaussian"){
    res<-list(data=data,param=allparam,sd=paramDist,probMat=EstModel)
  }
  if(family=="nbinomial"){
    res<-list(data=data,param=allparam,size=paramDist,probMat=EstModel)
  }
  if(family=="probit"){
    res<-list(data=data,param=allparam,probMat=pnorm(EstModel))
  }
  if(family=="poisson"){
    res<-list(data=data,param=allparam,probMat=exp(EstModel))
  }
  class(res)<-"communitySimul"
  
  return(res)
}

###########################################.....OWN FUNCTIONS........##################################################################

#FOLDS
#Transform dataframe into raster (requirement to use blockCV)

df2raster <- function(df) {
  df<- cbind(df[,c(1,2,4,5)])
  #create raster
  library(raster)
  dfr <- rasterFromXYZ(df)  #Convert first two columns as lon-lat and third as value                
  
  extent(dfr) <- c(0, 5000, 0, 5000)
  terra::crs(dfr) <- "EPSG:4978"
  dfr
}

#Group our data into different blocks (they will be used in the cross validation)

custom_blockCV <- function(df, n = 5, i = 100) {
  library(blockCV)
  library(sf)
  awt <- df2raster(df)
  pb_data <- sp::SpatialPoints(rasterToPoints(awt))
  crs(pb_data) <- "EPSG:4978"
  pb_data <- st_zm(st_as_sf(pb_data))
  
#Spatial blocking by specified range with random assignment
  sb <- spatialBlock(speciesData = pb_data,
                     rasterLayer = awt,
                     theRange = 1500,
                     k = n,
                     selection = "random",
                     iteration = i, # find evenly dispersed folds
  )
  data.frame(partition = sb$foldID)
}

#SUbsampling the original community including false absences, i.e., 0)

get.sample <- function(coo_matrix, p, sp_names = paste0("sp", 1:10), coord_names = c("x_coord", "y_coord"), var_names = c("x1", "x2"), size = FALSE, i) {
  #Create twin matrix with NA
  na_matrix<- coo_matrix
  na_matrix[,]<- NA
  
  #Identificate presence points of different species by colums (vectors) and charge them into empty vector
  sampling <- function(coo_v, p) {
    na_v <- coo_v
    na_v[] <- NA
    sp <- which(coo_v == 1)
    sp_p <- round(length(sp) * p)
    sp_muestra <- sample(sp, sp_p)
    na_v[sp_muestra] <- 1
    na_v
  }
  na_matrix <- apply(coo_matrix[, sp_names], 2, FUN = sampling, p = p)
  na_matrix[is.na(na_matrix)] <- 0
  no_nulo<- which(apply(na_matrix, 1, sum) != 0)
  na_matrix <- cbind(coo_matrix[no_nulo, c(coord_names, "interc", var_names)], na_matrix[no_nulo, ], coo_matrix[no_nulo, "partition"])
  colnames(na_matrix)[16] <- "partition"
  if(size != FALSE){
    set.seed(i)
    sampled_matrix <- na_matrix[sample(nrow(na_matrix), size), ]
  }else{
    sampled_matrix <- na_matrix
  }
  sampled_matrix
}

#Generate artificial data, considering or not the species correlations effect

generate.artificial.data <- function(GRID.SIZE,  spCor = NULL) {
  ENVIRONMENT <- melt(array(NA, dim=c(GRID.SIZE, GRID.SIZE)))[,-3]
  COORDS <- ENVIRONMENT
  colnames(COORDS) <- c("x_coord", "y_coord")
  ENVIRONMENT <- as.data.frame(scale(ENVIRONMENT))
  ENVIRONMENT.NAMES <- c("x1", "x2")
  colnames(ENVIRONMENT) <- ENVIRONMENT.NAMES
  
  desc <- cbind(1,ENVIRONMENT)
  colnames(desc) <- c("interc",  "x1", "x2")
  
  ### Simulate communities with 10 species
  set.seed(1)
  if(is.null(spCor)){
    comm <- communitySimul(desc, nsp = 10)
  }else{
    comm <- communitySimul(desc, nsp = 10, spCor = spCor)
  }
  comm$data$COORDS <- COORDS
  
  ALL <- cbind(comm$data$COORDS, comm$data$X, comm$data$Y)
  comm$data$partition <- custom_blockCV(ALL)
  comm
}


#Building and calibrating hmsc model

calibrate.hmsc <- function(data, thin, samples, transient, nChains, nParallel, sp_names = paste0("sp", 1:10), var_names = c("x1", "x2"), coord_names = c("x_coord", "y_coord")) {
  
  set.model <- function(data, sp_names, var_names, coord_names) {
    Y <- data[, sp_names] 
    XData<- data[, var_names]
    xy<- data[, coord_names]
    Y <- as.matrix(Y)
    rownames(Y) <- rownames(XData)
    XFormula = ~x1 + x2
    studyDesign <- data.frame(units=rownames(Y))
    studyDesign$units <- as.factor(studyDesign$units)
    rl <- HmscRandomLevel(sData=xy, sMethod = 'NNGP', nNeighbours = 10)
    rl = setPriors(rl, nfMin=1, nfMax=1)
    Hmsc(Y=Y,
         XData = XData,
         XFormula=XFormula,
         studyDesign = studyDesign,
         ranLevels = list("units"=rl),
         distr="probit")
  }
  
  model.pa<- set.model(data, sp_names, var_names, coord_names)
  
  # fitting the model
  sampleMcmc(model.pa, thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = nParallel)
}

#Checking the effective sample size (ess) and the Gelman diagnostics (gd)
plot.diagnostics <- function(model.pa, p) {
  mpost = convertToCodaObject(model.pa)
  
  # Note that Hmsc uses coda only to examine convergence
  
  
  # We will explore the MCMC convergence for the beta-parameters (species niches)and
  # V-parameters (variation in species niches),
  # We first examine the effective size of the posterior sample.
  ess.beta = effectiveSize(mpost$Beta)
  
  ess.V = effectiveSize(mpost$V)
  
  # We then examine the Gelman diagnostics, i.e. the Potential scale reduction factors
  # This diagnostic compares if different chains (here we have 2 chains) give consistent results
  # Ideally the value of this diagnostic would be close to one.
  # As you increase thinning, you should observe the values getting closer to one.
  psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
  
  psrf.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
 
  png(paste0("./diagnostics/diagnostics_", p, ".png"))
  par(mfrow = c(2, 2))
  hist(ess.beta, xlab = expression("Effective sample size" ~ beta ~ ""), main=NULL, col= "light green", breaks= 6)
  hist(ess.V, xlab = expression("Effective sample size" ~ v ~ ""),  main=NULL, breaks= 6)
  hist(psrf.beta, xlab = expression("Potential scale reduction factor" ~ beta ~ ""), main=NULL, col= "light green", breaks= 6)
  hist(psrf.V, xlab = expression("Potential scale reduction factor" ~ v ~ ""),  main=NULL, breaks= 6)
  dev.off()
}


median2 = function(x) {
  return(median(x, na.rm = TRUE))
}

mean2 = function(x) {
  return(mean(x, na.rm = TRUE))
}

#Required for sensitivity analysis in performance metrics calculation

mean_prediction <- function(predY, hM) {

  mPredY = matrix(NA, nrow = nrow(predY), ncol = ncol(predY))
  sel = hM$distr[, 1] == 3
  if (sum(sel) > 0) {
    mPredY[, sel] = as.matrix(apply(abind(predY[, sel, , 
                                                drop = FALSE], along = 3), c(1, 2), median2))
  }
  sel = !hM$distr[, 1] == 3
  if (sum(sel) > 0) {
    mPredY[, sel] = as.matrix(apply(abind(predY[, sel, , 
                                                drop = FALSE], along = 3), c(1, 2), mean2))
  }
  colnames(mPredY) <- colnames(predY)
  rownames(mPredY) <- rownames(predY)
  mPredY
}

#modification of Hmsc function, making it more flexible for validation with external data

custom_evaluateModelFit <- function (obsY, predY, hM) {
  computeRMSE = function(obs, pred) {
    ns = dim(obs)[2]
    RMSE = rep(NA, ns)
    for (i in 1:ns) {
      RMSE[i] = sqrt(mean((obs[, i] - pred[, i])^2, na.rm = TRUE))
    }
    return(RMSE)
  }
  computeR2 = function(obs, pred, method = "pearson") {
    ns = dim(obs)[2]
    R2 = rep(NA, ns)
    for (i in 1:ns) {
      co = cor(obs[, i], pred[, i], method = method, use = "pairwise")
      R2[i] = sign(co) * co^2
    }
    return(R2)
  }
  computeAUC = function(obs, pred) {
    ns = dim(obs)[2]
    AUC = rep(NA, ns)
    obsY <- ifelse(obs > 0, 1, 0)
    for (i in 1:ns) {
      sel = !is.na(obs[, i])
      if (length(unique(obs[sel, i])) == 2) 
        AUC[i] = pROC::auc(obs[sel, i], pred[sel, i], levels = c(0, 
                                                          1), direction = "<")
    }
    return(AUC)
  }
  computeTjurR2 = function(obs, pred) {
    ns = dim(obs)[2]
    R2 = rep(NA, ns)
    for (i in 1:ns) {
      R2[i] = mean(pred[which(obs[, i] == 1), i]) - mean(pred[which(obs[, 
                                                                      i] == 0), i])
    }
    return(R2)
  }

  mPredY = matrix(NA, nrow = nrow(predY), ncol = ncol(predY))
  sel = hM$distr[, 1] == 3
  if (sum(sel) > 0) {
    mPredY[, sel] = as.matrix(apply(abind(predY[, sel, , 
                                                drop = FALSE], along = 3), c(1, 2), median2))
  }
  sel = !hM$distr[, 1] == 3
  if (sum(sel) > 0) {
    mPredY[, sel] = as.matrix(apply(abind(predY[, sel, , 
                                                drop = FALSE], along = 3), c(1, 2), mean2))
  }
  RMSE = computeRMSE(obsY, mPredY)
  R2 = NULL
  AUC = NULL
  TjurR2 = NULL
  SR2 = NULL
  O.AUC = NULL
  O.TjurR2 = NULL
  O.RMSE = NULL
  C.SR2 = NULL
  C.RMSE = NULL
  sel = hM$distr[, 1] == 1
  if (sum(sel) > 0) {
    R2 = rep(NA, hM$ns)
    R2[sel] = computeR2(obsY[, sel, drop = FALSE], mPredY[, 
                                                          sel, drop = FALSE])
  }
  sel = hM$distr[, 1] == 2
  if (sum(sel) > 0) {
    AUC = rep(NA, hM$ns)
    TjurR2 = rep(NA, hM$ns)
    AUC[sel] = computeAUC(obsY[, sel, drop = FALSE], mPredY[, 
                                                            sel, drop = FALSE])
    TjurR2[sel] = computeTjurR2(obsY[, sel, drop = FALSE], 
                                mPredY[, sel, drop = FALSE])
  }
  sel = hM$distr[, 1] == 3
  if (sum(sel) > 0) {
    SR2 = rep(NA, hM$ns)
    O.AUC = rep(NA, hM$ns)
    O.TjurR2 = rep(NA, hM$ns)
    O.RMSE = rep(NA, hM$ns)
    C.SR2 = rep(NA, hM$ns)
    C.RMSE = rep(NA, hM$ns)
    SR2[sel] = computeR2(obsY[, sel, drop = FALSE], mPredY[, 
                                                           sel, drop = FALSE], method = "spearman")
    predO = 1 * (predY[, sel, , drop = FALSE] > 0)
    mPredO = as.matrix(apply(abind(predO, along = 3), c(1, 
                                                        2), mean2))
    O.AUC[sel] = computeAUC(1 * (obsY[, sel, drop = FALSE] > 
                                   0), mPredO)
    O.TjurR2[sel] = computeTjurR2(1 * (obsY[, sel, drop = FALSE] > 
                                         0), mPredO)
    O.RMSE[sel] = computeRMSE(1 * (obsY[, sel, drop = FALSE] > 
                                     0), mPredO)
    mPredCY = mPredY[, sel, drop = FALSE]/mPredO
    CY = obsY[, sel, drop = FALSE]
    CY[CY == 0] = NA
    C.SR2[sel] = computeR2(CY, mPredCY, method = "spearman")
    C.RMSE[sel] = computeRMSE(CY, mPredCY)
  }
  MF = list(RMSE = RMSE)
  if (!is.null(R2)) {
    MF$R2 = R2
  }
  if (!is.null(AUC)) {
    MF$AUC = AUC
  }
  if (!is.null(TjurR2)) {
    MF$TjurR2 = TjurR2
  }
  if (!is.null(SR2)) {
    MF$SR2 = SR2
  }
  if (!is.null(O.AUC)) {
    MF$O.AUC = O.AUC
  }
  if (!is.null(O.TjurR2)) {
    MF$O.TjurR2 = O.TjurR2
  }
  if (!is.null(O.RMSE)) {
    MF$O.RMSE = O.RMSE
  }
  if (!is.null(C.SR2)) {
    MF$C.SR2 = C.SR2
  }
  if (!is.null(C.RMSE)) {
    MF$C.RMSE = C.RMSE
  }
  return(MF)
}

#modification of Hmsc function, making it more flexible for validation with external data

tuning_computePredictedValues <- function(newData, oldData, modelo, Yc=NULL) {
  require(dplyr)
  new_index <- !rownames(newData) %in% rownames(oldData)
  old_index <- rownames(newData) %in% rownames(oldData)
  
  old_sites <- rownames(newData)[old_index]
  rL1 = modelo$ranLevels[[1]]
  
  if(sum(old_index) > 0) {
    StudyDesignOld <- filter(modelo$studyDesign, units %in% old_sites)
    xyOld = filter(rL1$s, rownames(rL1$s) %in% old_sites)
  }
  
  if(sum(new_index) > 0) {
    com <- newData[new_index,]
    grad <- prepareGradient(modelo, XDataNew = com[, c("x1", "x2")], sDataNew = list(units = com[, c("x_coord", "y_coord")]))
    xyNew = com[,1:2] 
    rownames(xyNew) = grad$StudyDesignNew[,1]
    colnames(xyNew) = colnames(xyOld)
  } 

  if(sum(new_index) > 0 & sum(old_index) > 0){
    StudyDesignAll <- rbind(StudyDesignOld, grad$studyDesignNew)
    xy <- rbind(xyOld, xyNew)
  }
  if(sum(new_index) > 0 & sum(old_index) == 0){
    StudyDesignAll <- grad$studyDesignNew
    xy <- xyNew
  }
  if(sum(new_index) == 0 & sum(old_index) > 0){
    StudyDesignAll <- StudyDesignOld
    xy <- xyOld
  }

  rL1$pi = StudyDesignAll[,1]
  rL1$s = xy
  rL1$N = nrow(rL1$s)
  
  if(is.null(Yc)){
    predYR1 = predict(modelo, studyDesign = StudyDesignAll, XData = newData[, c("x1", "x2")], ranLevels = list("units"=rL1),
                    expected = TRUE, predictEtaMean = TRUE)
  }else{ 
   Yc[,1] <- NA
     predYR1 = predict.Hmsc(modelo, studyDesign = StudyDesignAll, XData = newData[, c("x1", "x2")], ranLevels = list("units"=rL1),
                           expected = TRUE, predictEtaMean = TRUE, Yc=as.matrix(Yc))
    
  }
  predArray = abind(predYR1, along = 3)
  predArray
}

#sensitivity specificty analysis for testing model performance

sensibilidad_especifidad <- function(sp_names, pred, obs, tr) {

  pred <- as.data.frame(pred)
  pred$joinID <- rownames(pred)
  obs$joinID <- rownames(obs)
  
  original_cut<- semi_join(obs,pred, by= c("joinID"="joinID"))
  
  presencias_originalsp <- original_cut[which(original_cut[,sp_names]==1),]
  ausencias_originalsp <- original_cut[which(original_cut[,sp_names]==0),]
  
  presenciaspredsp<- pred
  presenciaspredsp<- presenciaspredsp[presencias_originalsp$joinID, ]
  presenciaspredsp<- presenciaspredsp[,sp_names] 
  
  ausenciaspredsp <- pred
  ausenciaspredsp <- ausenciaspredsp[ausencias_originalsp$joinID, ]
  ausenciaspredsp <- ausenciaspredsp[,sp_names] 
  
  e <- dismo::evaluate(p=presenciaspredsp, a=ausenciaspredsp, tr=tr)
  t <- threshold(e, "spec_sens")
  e <- dismo::evaluate(p=presenciaspredsp, a=ausenciaspredsp, tr=t)
  return(e)
}

extraccion <- function(e){
  e <- data.frame(e@TPR, e@TNR, e@kappa, e@t)
   # e[,"TSS"]<- e@TPR + e@TNR - 1
  colnames(e) <- c("TPR", "TNR", "kappa", "t")
  return(e)
}

#Gives performance metrics (threshold dependent and independent) for explanatory power

tuning.metrics<- function(com_testing, com_training, model, sp_names, p, Yc=NULL, i) {
 
  predY<- tuning_computePredictedValues(com_testing, com_training, model, Yc=Yc)
  mpredY<- mean_prediction(predY, model)
  metricas<- data.frame(custom_evaluateModelFit(obsY = com_testing[, sp_names], predY = predY, hM = model), species = sp_names, p = p)
   
  sensitivity <- lapply(sp_names, FUN= sensibilidad_especifidad, mpredY, com_testing)
  sensitivity <- do.call("rbind", lapply(sensitivity, FUN= extraccion))
  
  metricas <- cbind(metricas, sensitivity)
  metricas[,"TSS"]<- metricas$TPR + metricas$TNR - 1
  metricas[,"i"] <- i
  metricas <- reshape2::melt(metricas, measure.vars = c("RMSE", "AUC", "TjurR2", "TPR", "TNR", "kappa", "TSS"))
  metricas
}

#Gives performance metrics (threshold dependent and independent) for predictive power

tuning.predictive <- function(com_testing, model, sp_names, p, i) {
  predY = computePredictedValues(hM= model, partition= com_testing$partition , nParallel = 5)
  mpredY<- data.frame(mean_prediction(predY, model))
  colnames(mpredY)<- sp_names
  rownames(mpredY)<- rownames(com_testing)
  
  PP <- data.frame(evaluateModelFit(hM=model, predY=predY), species=sp_names, p=p)
  
  sensitivity <- lapply(sp_names, FUN= sensibilidad_especifidad, mpredY, com_testing)
  sensitivity <- do.call("rbind", lapply(sensitivity, FUN= extraccion))
  
  PP <- cbind(PP, sensitivity)
  PP[,"TSS"]<- PP$TPR + PP$TNR - 1
  PP[,"i"] <- i
  PP <- reshape2::melt(PP, measure.vars = c("RMSE", "AUC", "TjurR2", "TPR", "TNR", "kappa", "TSS"))
  PP
}


#####original cut####### Functions to get communities (without false absences) with the same size and the same sample as the subsampled ones

  original.same.sample <- function(predictions, originals) {
    predictions[,"rownames"] <- rownames(predictions)
    originals[,"rownames"] <- rownames(originals)
    
    original_cut<- semi_join(originals, predictions, by = c("rownames"="rownames"))
  }


original_same_size <- function(predictions, originals) {
  comunidad_original_mismo_tamano<- originals %>%
    sample_n(size=(nrow(predictions)),replace=FALSE)
}


#########PLOTS##########

#Plot to check the diferences between the different metrics (i.e., explanatory and predicitive power, same sample)

comparative.plot <- function(explanatory, predictive, same_sample) {
  
  same_sample <- same_sample[which(same_sample[,"variable"] == c("AUC", "TjurR2", "TSS", "kappa")),]
  explanatory <- explanatory[which(explanatory[,"variable"] == c("AUC", "TjurR2", "TSS", "kappa")),]
  predictive <- predictive[which(predictive[,"variable"] == c("AUC", "TjurR2", "TSS", "kappa")),]

  explanatory[,"validation"]<- "EXP"
  predictive[,"validation"]<- "PRED"
  same_sample[,"validation"]<- "EXPwfn"
  
  mix <- rbind(explanatory, predictive, same_sample) #bind here the metrics you want to compare
  mix <- mix %>% mutate(variable = recode(variable, "kappa" = "Kappa"))
  mix$validation <- factor(mix$validation, levels = c("EXP", "PRED", "EXPwfn"))
  mix<- rename(mix, Validation = validation)
  cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9")
                 
  ggplot(mix, aes(x=factor(p, level = c(0.9, 0.75, 0.5, 0.25, 0.1), labels = c("90%", "75%", "50%", "25%", "10%")), y=value)) +
    # scale_x_reverse() +
    geom_boxplot(aes(fill=Validation)) +
    scale_fill_manual(values=cbPalette) +
    scale_y_continuous(limits= c(0,1)) +
    facet_wrap(.~variable) +
    xlab("Retention percentage") +
    ylab("") +
    theme_bw() 
  
  ggplot2::ggsave("./figures/performance_figure.pdf", device = "pdf", width = 7, height = 4.5)
}

#Plot to draw the intercept and the slope of the beta parameters along retention percentage.

betas.plot <- function(combined_rmas) {
  rmas <- combined_rmas %>% pivot_longer(cols = c("intercept", "slope"))
  rmas_sum <- rmas %>% group_by (p, variable, name) %>% summarise(mean = mean(value), sd = sd(value))
  
  complete_labels <- c("Intercept", "Beta x1", "Beta x2")
  names(complete_labels) <- c("int", "x1", "x2")
  
  ggplot(rmas, aes(x = p*100, y = value, color = name)) +
    scale_x_reverse() +
    geom_jitter(width = 3, alpha = 0.3, size = 1.5) +
    geom_line(data = rmas_sum, aes(x = p*100, y = mean), size = 0.8) +
    geom_point(data = rmas_sum, aes(x = p*100, y = mean), size = 2.5) +
    facet_grid(~variable, labeller = as_labeller(complete_labels)) +
    scale_color_manual(values=c("#762a83", "#1b7837"),
                       labels = c("Intercepts", "Slopes"),
                       guide = "legend",
                       name = "RMAs coeffs") +
    xlab("Retention percentage (%)") +
    ylab("") +
    theme_bw() +
    theme(legend.position="right")
  ggplot2::ggsave("./figures/betas_figure.pdf", device = "pdf", width = 7.6, height = 2.1)
}

#Representation of the slope and intercept of the beta parameters, grouped by retention percentage

supl.betas.plot <- function(betas, rmas) {
  
  complete_labels <- c(int = "Intercept", x1 = "Beta x1", x2 = "Beta x2")
  names(complete_labels) <- c("int", "x1", "x2")
  
  p_labels <- c("10%", "25%", "50%", "75%", "90%")
  names(p_labels) <- c(0.1, 0.25, 0.5, 0.75, 0.9)
 
  ggplot(betas, aes_string(x="original", y="model")) +
    geom_hex() +
    scale_fill_gradient2(low = "#E69F00" , mid = "#56B4E9" , high = "#0072B2", midpoint = 3.5) +
    theme(panel.grid=element_blank()) +
    theme_bw() +
    geom_abline(linetype="dashed", colour="#999999", size=0.75) +
    scale_x_continuous(breaks=c(-1,0,1), name="Betas (original matrix)") +
    scale_y_continuous(breaks=c(-2,-1,0,1), name="Betas (model)") +
    geom_abline(aes(intercept=intercept, slope=slope), data=rmas, color="#CC79A7", size=0.75) +
    facet_grid(p ~ variable, labeller = labeller(variable = complete_labels, p = p_labels)) +
    geom_text(aes(-0.5, -1.85, label=paste("b=", formatC(round(intercept, 2), format="f", digits=2), "; m=", formatC(round(slope, 2), format="f", digits=2), "\nRsquare=", formatC(round(Rsquare, 2), format="f", digits=2), sep=""), group=NULL), size=3, data=rmas, hjust=0) +
    geom_text(aes(-0.7, 0.6, label=paste("RMSE=", formatC(round(RMSE, 2), format="f", digits=2), sep=""), group=NULL), size=3, data=rmas, hjust=1)
  
  ggplot2::ggsave("./figures/supl_betas_figure.pdf", device = "pdf", width = 7.5, height = 8)
}

#Correlation plot of the residuals given by the model for a retention percentage.

correlation.plot <- function(model, p, plotOrder) {
  OmegaCor = computeAssociations(model)
  supportLevel = 0.95
  
  toPlot = OmegaCor[[1]]$mean
  pdf(paste0("./figures/residuals_plot", p, ".pdf"))
  corrplot(toPlot[plotOrder,plotOrder], method = "number", tl.cex=0.5,
           col=colorRampPalette(c("blue", "yellow", "red"))(255))
  dev.off()
}

#Correlation plot of the coefficients of the original community

original.correlation.plot <- function(comm, plotOrder) {
  toPlot = comm$param$spCor
  pdf("./figures/original_residuals_plot.pdf")
  corrplot(toPlot[plotOrder,plotOrder], method = "number", tl.cex=0.5,
           col=colorRampPalette(c("blue", "yellow", "red"))(255))
  dev.off()
}

to_plot <- function(model, p) {
  OmegaCor = computeAssociations(model)
  supportLevel = 0.95
  
  toPlot = OmegaCor[[1]]$mean
  toPlot <- melt(toPlot)
  toPlot$p_target <- p
  toPlot <- toPlot
  toPlot
}

individual.corr.plot <- function(p, array_corr, plotOrder) {
  pdf(paste0("residuals_plot", p, ".pdf"))
  p <- as.character(p)
  array_corr <- array_corr[,,p]
  corrplot(array_corr[plotOrder,plotOrder], method = "number", tl.cex=0.5,
           col=colorRampPalette(c("blue", "yellow", "red"))(255))
  dev.off()
}


####################

####RMAS BETAS######

BETAS <- function(data, model, p, i) {
  #original
  beta.original <- data.frame(data$param$paramX)
  colnames(beta.original) <- c("int", "x1", "x2")
  beta.original <- melt(beta.original, measure.vars= c("int", "x1", "x2"))
  beta.original$species <- paste0("sp", 1:10)
  colnames(beta.original) <- c("variable", "original", "species")
  beta.original$p <- p
  
  #model
  postBeta <- getPostEstimate(model, parName="Beta")
  
  beta.mod <- data.frame(t(postBeta$mean))
  colnames(beta.mod) <- c("int", "x1", "x2")
  beta.mod <- melt(beta.mod, measure.vars= c("int", "x1", "x2"))
  beta.mod$species <- paste0("sp", 1:10)
  colnames(beta.mod) <- c("variable", "model", "species")
  beta.mod$p <- p
  
  betas <- merge(beta.original, beta.mod, by = c("variable", "species", "p"))
  betas[,"i"] <- i
  return(betas)
}

#calibration of rma models 
rmaFit <- function(dat, a, b, var, p_tar, i_tar){
  dat<- filter(dat, variable == var & p == p_tar & i == i_tar)
  require(lmodel2)
  form <- as.formula(paste(b, "~", a, sep=""))
  mod <- lmodel2(form, data=dat, "interval", "interval", 1000)
  reg <- mod$regression.results
  names(reg) <- c("method", "intercept", "slope", "angle", "p-value")
  reg$Rsquare <- mod$rsquare
  reg$Pparam <- mod$P.param
  # Calculate RMSE. I keep this into this loop to get one estimate for each model fitted.
  e <- dat[,a] - dat[,b]
  reg$RMSE <- sqrt(mean(e^2))
  reg$p <- p_tar
  reg$i<- i_tar
  return(reg[4,-1]) #row 4 , rma method
}

RMAS <- function(betas, p_tar, i_tar) {
  int.rma<- as.data.frame(rmaFit(betas, "original", "model", "int", p_tar, i_tar))
  x1.rma<- as.data.frame(rmaFit(betas, "original", "model", "x1", p_tar, i_tar))
  x2.rma<- as.data.frame(rmaFit(betas, "original", "model", "x2", p_tar, i_tar))
  rmas<- rbind(int.rma, x1.rma, x2.rma)
  rmas$variable<- c("int", "x1", "x2")
  return(rmas)
}


supl_rmaFit <- function(dat, a, b, var, p_tar){
  dat<- filter(dat, variable == var & p == p_tar)
  require(lmodel2)
  form <- as.formula(paste(b, "~", a, sep=""))
  mod <- lmodel2(form, data=dat, "interval", "interval", 1000)
  reg <- mod$regression.results
  names(reg) <- c("method", "intercept", "slope", "angle", "p-value")
  reg$Rsquare <- mod$rsquare
  reg$Pparam <- mod$P.param
  # Calculate RMSE. I keep this into this loop to get one estimate for each model fitted.
  e <- dat[,a] - dat[,b]
  reg$RMSE <- sqrt(mean(e^2))
  reg$p <- p_tar
  return(reg[4,-1]) #row 4 , rma method
}

supl.RMAS <- function(betas, p_tar) {
  int.rma<- as.data.frame(supl_rmaFit(betas, "original", "model", "int", p_tar))
  x1.rma<- as.data.frame(supl_rmaFit(betas, "original", "model", "x1", p_tar))
  x2.rma<- as.data.frame(supl_rmaFit(betas, "original", "model", "x2", p_tar))
  rmas<- rbind(int.rma, x1.rma, x2.rma)
  rmas$variable<- c("int", "x1", "x2")
  return(rmas)
}

#########################

coo_pattern <- function(community_data, sp_names, sample, i, real_matrix = FALSE) {
  if(real_matrix == TRUE){
    prob <- cooccur(mat = t(community_data[,sp_names]), type = "spp_site", thresh = TRUE,
                    spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE, eff_matrix = TRUE)
  }
  else{prob <- community_data}
  prob_matrix <- as.matrix(prob)
  prob_matrix[upper.tri(prob_matrix, diag = TRUE)] <- NA
  prob_matrix <- na.omit(melt(prob_matrix))
  prob_matrix <- cbind(prob_matrix, paste0(prob_matrix[,1], "-", prob_matrix[,2]))
  prob_matrix[,5]<-  sample
  prob_matrix[,"i"] <- i
  colnames(prob_matrix) <- c("var1", "var2", "prob", "sp-sp", "sample", "i")
  prob_matrix
}


pattern.plot <- function(coo_plot, coo_levels) {
  coo_plot$sample <- factor(coo_plot$sample)
  coo_plot$sample <-factor(coo_plot$sample, levels = sort(levels(coo_plot$sample), decreasing= TRUE))
  
  coo_plot$`sp-sp` <- factor(coo_plot$`sp-sp`)
  coo_plot$`sp-sp` <-factor(coo_plot$`sp-sp`, levels = coo_levels)
  
  ggplot(coo_plot) + 
    geom_tile(aes(x= `sp-sp`, y= sample, fill = prob)) +
    scale_fill_gradientn(colors= c("#8e0152", "#c51b7d", "#f5f5f5", "#4d9221", "#276419"), values = c(0, 0.4, 0.5, 0.6, 1), limits= c(-0.11, 0.11)) +
    xlab("") +
    ylab("") +
    scale_y_discrete(labels=c(bquote(CM[O2]), bquote(CM[S90]), bquote(CM[S75]), bquote(CM[S50]), bquote(CM[S25]), bquote(CM[S10]))) +
    theme_bw() +
    theme(legend.title= element_blank(), axis.text.x=element_text(angle = 45, hjust = 1))
  ggplot2::ggsave("./figures/coo_figure.pdf", device = "pdf", width = 10, height = 4)
}


############################

corr.extraction <- function(model) {
  OmegaCor <- computeAssociations(model)
  corr <- OmegaCor[[1]]$mean
  corr
}

corr.pattern.plot <- function(corr_plot, corr_levels, filename, sign = FALSE) {
  corr_plot$sample <- factor(corr_plot$sample)
  corr_plot$sample <-factor(corr_plot$sample, levels = sort(levels(corr_plot$sample), decreasing= TRUE))
  levels(corr_plot$sample)   
  
  # corr_plot$`sp-sp` <- factor(corr_plot$`sp-sp`)
  corr_plot$`sp-sp` <- factor(corr_plot$`sp-sp`, levels = corr_levels)
  if(sign == TRUE){ corr_plot$prob <- sign(corr_plot$prob) }
  ggplot(corr_plot) + 
    geom_raster(aes(x = `sp-sp`, y = sample, fill = prob)) +
    scale_fill_gradientn(colors= c("#7f3b08", "#e08214", "#f5f5f5", "#8073ac", "#2d004b"), values = c(0, 0.35, 0.5, 0.65, 1), limits= c(-1, 1)) +
    ylab("") + 
    xlab("") +
    scale_y_discrete(labels=c(expression(Omega), bquote(CM[O2]), bquote(CM[S90]), bquote(CM[S75]), bquote(CM[S50]), bquote(CM[S25]), bquote(CM[S10]))) +
    geom_rect(xmin = -Inf, xmax = +Inf,   ymin = -Inf, ymax = 1.5,   color = "black", alpha = 0, linetype = "dashed") +
    theme_bw() +
    theme(legend.title= element_blank(), axis.text.x=element_text(angle = 45, hjust = 1))
  # filename <- "./figures/omegas.pdf"
  ggplot2::ggsave(filename, device = "pdf", width = 10, height = 4)
}


######false negatives count

compare.matrices <- function(sp_names, sample, original, i, p) {

  sample <- sample[rownames(original),]
  
  false_negatives <- apply(sample[,sp_names] != original[,sp_names], 2, sum) / nrow(sample)
  false_negatives <- data.frame(species = sp_names, FN = false_negatives, i = i, p= p)

  return(false_negatives)
}


####procrustes analysis

procrustes.analysis <- function(i_tar, p_tar, pattern, data) {
  sampled_coor <- reshape2::acast(pattern, var1 ~ var2 ~ sample ~ i, value.var = "prob" )
  sampled_coor<- drop(sampled_coor)
  sampled_coor <- rbind(sp1=NA, sampled_coor)
  sampled_coor <- cbind(sampled_coor, sp10=NA)
  sampled_coor[upper.tri(sampled_coor)] <- t(sampled_coor)[upper.tri(sampled_coor)]
  diag(sampled_coor) <- 1
  
  procrus_sign <- procrustes(sign(data$param$spCor), sign(sampled_coor), scale = F)
  procrus <- procrustes(data$param$spCor, sampled_coor, scale = F)
  return(data.frame(rmse= summary(procrus)$rmse, rmse_sign= summary(procrus_sign)$rmse, i= i_tar, p= p_tar))
}


#PREDICTIONS OVER SPACE

spatial.prediction <- function(p, i, env_xy_data, model, sp_names) {
  
  xy.grid <- env_xy_data %>% dplyr::select(c(x_coord, y_coord))
  colnames(xy.grid) <- c("x", "y")
  
  XData.grid <- env_xy_data %>% dplyr::select(c(x1, x2))
  
  # XData.grid <- as.data.frame(getValues(env_data))
  # xy.grid <- as.data.frame(coordinates(env_data))
  # 
  complete.index <- which(complete.cases(XData.grid))
  Gradient <- prepareGradient(model, XDataNew = XData.grid[complete.index,], sDataNew = list(units=xy.grid[complete.index,]))
  
  predY <- predict(model, Gradient=Gradient, expected = TRUE, predictEtaMean = TRUE, nParallel=2)
  
  length(predY)
  dim(predY[[1]])
  EpredY <- matrix(NA, nrow(XData.grid), 10)
  library(abind)
  EpredY[complete.index,] <- apply(abind(predY,along=3),c(1,2),mean)
  length(EpredY)
  
  mapData=data.frame(x=xy.grid[,1],y=xy.grid[,2], EpredY)
  colnames(mapData) <- c("x", "y", sp_names)
  library(reshape2)
  mapData <- melt(mapData, id.vars = c("x", "y")) %>% mutate(sample = p, iteration = i)
  
  # library(RColorBrewer)
  # library(ggplot2)
  # 
  # ggplot(data = mapData, aes(x=x, y=y, fill=value)) +
  #   geom_raster() +
  #   scale_fill_distiller(palette="RdYlGn", na.value="royalblue") +
  #   theme(legend.position="bottom", aspect.ratio= 1) +
  #   facet_wrap(~variable)
}

# spat.plot <- function(Data) {
#   distributions <- ggplot(data = Data, aes(x=x, y=y, fill=value)) +
#     geom_raster() +
#     scale_fill_distiller(palette="RdYlGn", na.value="royalblue") +
#     theme(legend.position="right", aspect.ratio= 1) +
#     facet_grid(variable~sample) +
#     theme(axis.text.x=element_blank(),
#           axis.ticks.x=element_blank()) +
#     labs(x = "")
#   
#   distributions
#   
#   ggplot2::ggsave("./figures/spatial_figure.pdf", device = "pdf", width = 7.5, height = 8)
#   
#   distributions
# }



spatial.clust.plot <- function(Data) {
  
  wide_data <- Data %>% pivot_wider(names_from = variable, values_from = value) 

  # sp_data <- wide_data %>% ungroup() %>% dplyr::select(paste0("sp", 1:10))
  # cls <- kmeans(x = sp_data, centers = 5)
  
  # set.seed(123)
  # fviz_nbclust(wide_data %>% filter(sample == 0.25) %>% ungroup() %>% dplyr::select(paste0("sp", 1:10)),
  #             kmeans, method = "wss")
  # fviz_nbclust(wide_data %>% filter(sample == 0.25) %>% ungroup() %>% dplyr::select(paste0("sp", 1:10)),
  #              kmeans, method = "silhouette")
  # fviz_nbclust(wide_data %>% filter(sample == 0.25) %>% ungroup() %>% dplyr::select(paste0("sp", 1:10)),
  #              kmeans, method = "gap_stat")
  # 
  calculate_clusters <- function(sample) {
    set.seed(1)
    cls <- kmeans(x = wide_data %>% 
                       filter(sample == sample) %>% 
                       ungroup() %>% 
                       dplyr::select(paste0("sp", 1:10)), centers = 3) %>% 
      .$cluster %>% 
      as.numeric()
    data.frame(value = cls, variable = "cluster", sample = sample) %>%
      bind_cols(wide_data %>% 
                  filter(sample == sample) %>% 
                  ungroup() %>% 
                  dplyr::select(c("x", "y")))    
  }             

  cls_001 <- calculate_clusters(0.1)
  cls_025 <- calculate_clusters(0.25)
  cls_050 <- calculate_clusters(0.5)
  cls_075 <- calculate_clusters(0.75)
  cls_090 <- calculate_clusters(0.9)
  cls_100 <- calculate_clusters(1)
  
  cls <- bind_rows(cls_001, cls_025, cls_050, cls_075, cls_090, cls_100)
  
  Data <- bind_rows(Data, cls)

  distributions <- Data %>% filter(variable != "cluster") %>% 
    ggplot(aes(x=x, y=y, fill=value)) +
    geom_raster() +
    # scale_fill_distiller(palette="RdYlGn", na.value="royalblue") +
    scale_fill_viridis_c(na.value="royalblue") +
    theme(legend.position="right", aspect.ratio= 1) +
    facet_grid(variable~sample) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y = element_text(size = 8)) +
    labs(x = "") 
  
  cluster <- Data %>% filter(variable == "cluster") %>%  
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = as.factor(value))) +
    scale_fill_manual("value", values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                 "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")) +
    theme(legend.position="right", aspect.ratio= 1) +
    facet_grid(variable~sample) +
    theme(strip.text.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
          axis.text.y = element_text(size = 8)) 
    
  
  library(ggplot2)
  library(ggpubr)
  theme_set(theme_pubr())
 
  ggarrange(distributions, cluster, ncol = 1, nrow = 2, heights = c(6.2,1), align = "v") %>%
    ggexport(filename = "./figures/spatial_figure.pdf", width = 7,5, height = 8)

}
