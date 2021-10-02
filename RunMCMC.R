library(DirichletReg)
library(jsonlite)
arg = commandArgs(trailingOnly=TRUE)[1]
config = read_json(arg)




for(k in 1:length(config$DATA_LOCATION)){

iter = config$MCMCITER
burnin = config$BURNIN
gamma.prior.a = config$GAMMAPRIORA
gamma.prior.b = config$GAMMAPRIORB
alpha = config$ALPHA
maxN  = config$MAXN

waiting_times1 <- read.table(config$DATA_LOCATION[[k]], quote="\"", comment.char="")
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
Posterior.samples.A = matrix(0,(iter - burnin), maxN)
Posterior.samples.lambda = matrix(0,(iter - burnin), maxN)


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
  
  Aj = Aj/AjrowSums
  
  
  for(j in 1:m){
    I.j.Posterior[j] = which(rmultinom(1, 1, Aj[j,]) ==1, arr.ind =T)[,1]
  }
  
  ## Sample A 
  A.Dirichlet.Posterior = NULL

  for(lam in 1:N.Posterior){
    A.Dirichlet.Posterior = c(A.Dirichlet.Posterior, sum(I.j.Posterior == lam))
  }

  N.Posterior = sum(A.Dirichlet.Posterior>0)

  A.Dirichlet.Posterior = A.Dirichlet.Posterior[A.Dirichlet.Posterior>0]
  lambda.Posterior = lambda.Posterior[A.Dirichlet.Posterior>0]


  if(N.Posterior< maxN){
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
}#done with MCMC
probcomput = rep(0,maxN)
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
output.data = cbind(var.name,var.est, var.stat )
write.csv(data.frame(output.data), config$SUMMARY_DATA_LOCATION[[k]])

}


