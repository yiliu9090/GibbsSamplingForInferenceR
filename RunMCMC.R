library(DirichletReg)
library(jsonlite)
library(matrixStats)
library(Rcpp)
sourceCpp("gsl.cpp")

#ptm <- proc.time()

arg = commandArgs(trailingOnly=TRUE)[1]
config = read_json(arg)

set.seed(config$SEED)

#This is 

AlphaSize = length(config$ALPHA)
print(AlphaSize)

for(k in 1:length(config$DATA_LOCATION)){

  iter = config$MCMCITER
  burnin = config$BURNIN

  gamma.prior.a = config$GAMMAPRIORA
  gamma.prior.b = config$GAMMAPRIORB



  maxN  = config$MAXN

  SepFac = config$SEPARATION_FACTOR[[k]][[1]]

  waiting_times1 <- read.table(config$DATA_LOCATION[[k]], quote="\"", comment.char="")
  t = waiting_times1$V1
  m = length(t)
  ## Computing the Posterior 

  if(length(config$ALPHA)>1){
    alpha = config$ALPHA[[k]][[1]]
  }else{
    alpha = config$ALPHA[[1]]
  }



  #Initialize the value 
  N.Posterior = 1 
  lambda.Posterior = rgamma(1, gamma.prior.a , rate= gamma.prior.b)
  A.Posterior = c(1)
  I.j.Posterior = which(rmultinom(m, 1, A.Posterior) ==1, arr.ind =T)[,1] #classes

  #Initialize Prior 

  Posterior.sampes.N = rep(0, (iter - burnin))
  Posterior.samples.A = matrix(0,(iter - burnin), maxN)
  Posterior.samples.lambda = matrix(0,(iter - burnin), maxN)

  NewLambdaI = NULL

  #Likelihood data 
  likelihood_changes = NULL

  ## Gibbs sampling
  for( i in 1:iter){
  
    ## sample lambda give I(j) and t
    for(lam in 1:N.Posterior){
      gamma.posterior.a = (gamma.prior.a+sum(I.j.Posterior==lam))
      gamma.posterior.b = gamma.prior.b + sum(t*(I.j.Posterior==lam))
    
      if(lam >=2){

      ApproveI = 0 
      Trials = 0 
    
      while(ApproveI == 0 & Trials < 50){
        Slambda = rgamma(1, gamma.posterior.a , rate= gamma.posterior.b)# suggested lambda 
        Trials = Trials + 1
        ApproveI = 1
        for(u in 1:(lam-1)){
          if(Slambda < lambda.Posterior[u]*SepFac & Slambda >lambda.Posterior[u]/SepFac){
            ApproveI = 0
          }
        }
      }
      if(Trials < 50){
        lambda.Posterior[lam] = Slambda
      }else{
        lambda.Posterior[lam] = 0 #sample(c(min(lambda.Posterior[1:(lam-1)])/SepFac,max(lambda.Posterior[1:(lam-1)])*SepFac),1,c(0.5,0.5))
      }
    
      }else{
        #first lambda
        lambda.Posterior[lam] = rgamma(1, gamma.posterior.a , rate= gamma.posterior.b)
      }
    }  

    nonZeroIndexLambda <- which( lambda.Posterior>0)
    N.Posterior <- length(nonZeroIndexLambda) 
    A.Posterior <- A.Posterior[nonZeroIndexLambda]
    A.Posterior <- A.Posterior/(sum(A.Posterior))
    lambda.Posterior <-  lambda.Posterior[nonZeroIndexLambda]

    ## Sample I(j) for each data 
    Aj = matrix(0,m,N.Posterior)


    for(lam in 1:N.Posterior){
      Aj[,lam] = A.Posterior[lam]*lambda.Posterior[lam]*exp(-t*lambda.Posterior[lam])
    }
  

    Aj = Aj/rowMaxs(Aj)

    AjrowSums = rowSums(Aj)
 
    Aj = Aj/AjrowSums
  
    #print(Aj)
    D = gsl_mmm(Aj)

    I.j.Posterior = which(D==1,arr.ind = TRUE)[,2]
  
  
    #for(j in 1:m){
    #  I.j.Posterior[j] = which(rmultinom(1, 1, Aj[j,]) ==1, arr.ind =T)[,1]
    #}

    #print(max(I.j.Posterior))
    ## Sample A 
    A.Dirichlet.Posterior = NULL

    for(lam in 1:N.Posterior){
      if(sum(I.j.Posterior == lam)>0){
        A.Dirichlet.Posterior = c(A.Dirichlet.Posterior, sum(I.j.Posterior == lam)+1)
      }
      else{
        A.Dirichlet.Posterior = c(A.Dirichlet.Posterior, 0)
      }
    
    }

    N.Posterior = sum(A.Dirichlet.Posterior>0)
    ix = which(A.Dirichlet.Posterior>0)
    A.Dirichlet.Posterior = A.Dirichlet.Posterior[ix]
    lambda.Posterior = lambda.Posterior[ix]

    AddGamma = 0

    if(N.Posterior< maxN){
      ApproveI = 0 
      Trials = 0 
      while(ApproveI == 0 & Trials < 50){
        Slambda = rgamma(1, gamma.prior.a , rate= gamma.prior.b)# suggested lambda 
        Trials = Trials + 1
        ApproveI = 1
        for(u in 1:N.Posterior){
          if(Slambda <= lambda.Posterior[u]*SepFac & Slambda >= lambda.Posterior[u]/SepFac){
            ApproveI = 0
          }
        }
      }
    
      if(Trials < 50){
        AddGamma = 1 #Indicating that one lambda is added.
        A.Dirichlet.Posterior = c(A.Dirichlet.Posterior,alpha)
        N.Posterior = N.Posterior + 1
        lambda.Posterior = c(lambda.Posterior, Slambda)
      }
      else{
        A.Dirichlet.Posterior = A.Dirichlet.Posterior #a hard cut off so that all the results can be saved
      } 
    }else{
      A.Dirichlet.Posterior = A.Dirichlet.Posterior #a hard cut off so that all the results can be saved
    }

    #Reorder lambda so we do not have switch spaces problem 
    A.Posterior = rdirichlet(1,A.Dirichlet.Posterior)

    if(AddGamma == 1){
      #Add a gamma
      Nnow = sum(lambda.Posterior>0)
      o = order(lambda.Posterior[1:(Nnow-1)], decreasing= TRUE)
      A.PosteriorSort = A.Posterior[1:(Nnow-1)]
      A.Posterior[1:(Nnow-1)] = A.PosteriorSort[o]
      lambda.PosteriorSort = lambda.Posterior[1:(Nnow-1)]
      lambda.Posterior[1:(Nnow-1)] = lambda.PosteriorSort[o]
    }else{
  
      #No lambda added
      o = order(lambda.Posterior, decreasing= TRUE)
      A.Posterior = A.Posterior[o]
      lambda.Posterior = lambda.Posterior[o]

    }
  
  
  

    if(i%%100 ==0){
      L = A.Posterior[1]*lambda.Posterior[1]*exp(-t*lambda.Posterior[1])
      if(N.Posterior>=2){
        for(lam in 2:N.Posterior){
          L = L + A.Posterior[lam]*lambda.Posterior[lam]*exp(-t*lambda.Posterior[lam])
        }
      }
      likelihood_changes =c(likelihood_changes, sum(log(L)))
    }

    #collect Samples
    if(i > burnin){
    
      if(Trials < 50){
        NewLambdaI = c(NewLambdaI,1)
      }
      else{
        NewLambdaI = c(NewLambdaI,0)
      }

      Posterior.sampes.N[i - burnin] = N.Posterior
      if(N.Posterior< maxN){
        Posterior.samples.A[(i-burnin),] = c(A.Posterior,rep(0,maxN-N.Posterior))
        Posterior.samples.lambda[(i-burnin),] = c(lambda.Posterior,rep(0,maxN-length(lambda.Posterior)))
      }else{
        Posterior.samples.A[(i-burnin),] = A.Posterior
        Posterior.samples.lambda[(i-burnin),] = lambda.Posterior
      }
    
    }

    }#done with MCMC

    probcomput = rep(0,maxN)


    pdf(config$LIKELIHOOD_PLOT[[k]])
    plot((1:length(likelihood_changes))*100,likelihood_changes)
    dev.off()

    for(i in 1:maxN){
      probcomput[i] = mean(Posterior.sampes.N==i)
    }

    N.est = which.max(probcomput) - 1
    var.name = c('N')
    var.est  = c(N.est)
    var.stat = c(max(probcomput))

    for(i in 1:N.est){
      var.name = c(var.name, paste0("Alpha", as.character(i)))
      var.est  = c(var.est , mean(Posterior.samples.A[,i]))
      var.stat = c(var.stat, sd(Posterior.samples.A[,i]))
    }
    for(i in 1:N.est){
      var.name = c(var.name, paste0("Lambda", as.character(i)))
      var.est  = c(var.est , mean(Posterior.samples.lambda[,i]))
      var.stat = c(var.stat, sd(Posterior.samples.lambda[,i]))
    }
    pdf('Data/Output/Histplot.pdf')
    ix = which(Posterior.sampes.N==(N.est+1))
    hist(Posterior.samples.lambda[ix], breaks = 100)
    dev.off()
   



    output.data = cbind(var.name,var.est, var.stat )
    output.data = data.frame(output.data)
    output.data$var.est = as.numeric(output.data$var.est)
    output.data$var.stat = as.numeric(output.data$var.stat)
    write.csv(output.data, config$SUMMARY_DATA_LOCATION[[k]])
  
    if(config$DUMP){
      save(Posterior.sampes.N, Posterior.samples.A,Posterior.samples.lambda,  file = config$DUMP_LOCATION[[k]])
    }


} 


