library(DirichletReg)
library(json.lite)

arg = commandArgs(trailingOnly=TRUE)[1]
config = read_json(arg)




for(j in 1:length(config$DATA_LOCATION)){

iter = config$MCMCITER
burnin = config$BURNIN
gamma.prior.a = config$GAMMAPRIORA
gamma.prior.b = config$GAMMAPRIORB
alpha = config$ALPHA
maxN  = config$MAXN

waiting_times1 <- read.table(config$DATA_LOCATION[[j]], quote="\"", comment.char="")
t = waiting_times1$V1
m = length(t)
## Computing the Posterior 


#Initialize the value 
N.Posterior = 2 
lambda.Posterior = c(2,2)
A.Posterior = c(0.5, 0.5)
I.j.Posterior = which(rmultinom(m, 1, A.Posterior) ==1, arr.ind =T)[,1] #classes

#Initialize Prior 


#set seed

Posterior.sampes.N = rep(0, (iter - burnin))
Posterior.samples.A = matrix(0,(iter - burnin), 50)
Posterior.samples.lambda = matrix(0,(iter - burnin), 50)


## Gibbs sampling
for( i in 1:iter){
  
  
  ## sample lambda give I(j) and t
  for(lam in 1:N.Posterior){
    gamma.posterior.a = (gamma.prior.a+sum(I.j.Posterior==lam))
    gamma.posterior.b = gamma.prior.b + sum(t*(I.j.Posterior==lam))
    lambda.Posterior[lam] = rgamma(1, gamma.posterior.a , rate= gamma.posterior.b)
  }  


  ## Sample I(j) for each data 
  Aj = matrix(0,m,N.Posterior)
  for(lam in 1:N.Posterior){
    Aj[,lam] = A.Posterior[lam]*lambda.Posterior[lam]*exp(-t*lambda.Posterior[lam])
  }
  
  
  AjrowSums = rowSums(Aj)
  
  if(sum(AjrowSums==0)>0){
    error_indx = which(AjrowSums==0)
  }
  
  
  Aj = Aj/AjrowSums
  
  if(sum(AjrowSums==0)>0){
    Aj[error_indx,] = 1/nrow(Aj)
  }
  
  
  
  for(j in 1:m){
    I.j.Posterior[j] = which(rmultinom(1, 1, Aj[j,]) ==1, arr.ind =T)[,1]
  }
  
  ## Sample A 
  A.Dirichlet.Posterior = NULL
  N.Posterior.update = N.Posterior
  for(lam in 1:N.Posterior){
    A.Dirichlet.Posterior = c(A.Dirichlet.Posterior, sum(I.j.Posterior == lam))
  }

  N.Posterior = sum(A.Dirichlet.Posterior>0)

  A.Dirichlet.Posterior = A.Dirichlet.Posterior[A.Dirichlet.Posterior>0]
  lambda.Posterior = lambda.Posterior[A.Dirichlet.Posterior>0]


  if(N.Posterior< 50){
    A.Dirichlet.Posterior = c(A.Dirichlet.Posterior,alpha)
    N.Posterior = N.Posterior + 1
    lambda.Posterior = c(lambda.Posterior, rgamma(1, gamma.prior.a , rate= gamma.prior.b))
  }else{
    A.Dirichlet.Posterior = c(A.Dirichlet.Posterior) #a hard cut off so that all the results can be saved
  }
  
  
  A.Posterior = rdirichlet(1,A.Dirichlet.Posterior)
  
  o = order(A.Posterior, decreasing= TRUE)
  A.Posterior = A.Posterior[o]
  lambda.Posterior = lambda.Posterior[o]
  
  #collect Samples
  if(i > burnin){
    
    Posterior.sampes.N[i - burnin] = N.Posterior
    if(N.Posterior< maxN){
      Posterior.samples.A[(i-burnin),] = c(A.Posterior,rep(0,maxN-N.Posterior))
      Posterior.samples.lambda[(i-burnin),] = c(lambda.Posterior,rep(0,maxN-length(lambda.Posterior)))
    }else{
      Posterior.samples.A[(i-burnin),] = A.Posterior
      Posterior.samples.lambda[(i-burnin),] = lambda.Posterior
    }
    
  }
}


}