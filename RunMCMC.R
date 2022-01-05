library(DirichletReg)
library(jsonlite)
library(matrixStats)
library(Rcpp)
sourceCpp("gsl.cpp")

#V2: moving across



arg = commandArgs(trailingOnly=TRUE)[1]
config = read_json(arg)

set.seed(config$SEED)

#This is 

nsd = config$NSD

AlphaSize = length(config$ALPHA)

if(length(config$NAME)!=length(config$DATA_LOCATION)){
  stop("This is an error message")
}




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

  #if(length(config$ALPHA)>1){
  #  alpha = config$ALPHA[[k]][[1]]
  #}else{
  #  alpha = config$ALPHA[[1]]
  #}

  PosteriorSampleListN = list()
  PosteriorSampleListA = list()
  PosteriorSampleListL = list()
  PosteriorLikelihood  = list()
  MaxLikelihood = -Inf
  MinLikelihood = Inf
  Accepted = rep(0,AlphaSize)

  
  #
  NEST     = NULL #This is to keep an record of all the N.estimates
  var.name = NULL
  var.est  = NULL
  var.stat = NULL
  var.up   = NULL
  var.low  = NULL
  var.alpha= NULL
  varl     = NULL
  var.est.l= NULL
  var.est.N= NULL
  var.est.a= NULL
  var.est.s= NULL #standard deviation
  for(alphanumber in 1:AlphaSize){
    #Likelihood data 
    likelihood_changes = NULL

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

    alpha = config$ALPHA[[alphanumber]][[1]]
  
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
      
      }#done with recording the MCMC samples



    }#done with MCMC

    PosteriorLikelihood[[alphanumber]] = likelihood_changes
    MaxLikelihood = max(c(max(likelihood_changes, MaxLikelihood)))
    MinLikelihood = min(c(min(likelihood_changes, MinLikelihood)))

    #find the right set of values
    #Summarize the data for that particular alpha
    probcomput = rep(0,maxN)

    for(i in 1:maxN){
      probcomput[i] = mean(Posterior.sampes.N==i)
    }


    N.est = which.max(probcomput) - 1

    NEST = c(NEST, N.est)

    var.name = c(var.name, 'N')
    var.est  = c(var.est , N.est)
    var.stat = c(var.stat, max(probcomput))
    var.up   = c(var.up , 0)
    var.low  = c(var.low, 0)
    var.alpha = c(var.alpha, alpha)

    ix = which(Posterior.sampes.N == (N.est+1))

    for(i in 1:N.est){
      var.name = c(var.name, paste0("A", as.character(i)))
      var.est  = c(var.est , mean(Posterior.samples.A[ix,i]))
      var.stat = c(var.stat, sd(Posterior.samples.A[ix,i]))
      var.up   = c(var.up,   mean(Posterior.samples.A[ix,i]) + nsd*sd(Posterior.samples.A[ix,i]))
      var.low  = c(var.low,  mean(Posterior.samples.A[ix,i]) - nsd*sd(Posterior.samples.A[ix,i]))
      var.alpha = c(var.alpha, alpha)
    }
    lambda_samples = NULL 
    for(i in 1:N.est){
      varl     = c(varl, Posterior.samples.lambda[ix,i])
      lambda_samples = c(lambda_samples, Posterior.samples.lambda[ix,i]) #lambda samples
      var.name = c(var.name, paste0("Lambda", as.character(i)))
      var.est  = c(var.est , mean(Posterior.samples.lambda[ix,i]))
      var.stat = c(var.stat, sd(Posterior.samples.lambda[ix,i]))
      var.up   = c(var.up,   mean(Posterior.samples.lambda[ix,i]) + nsd*sd(Posterior.samples.lambda[ix,i]))
      var.low  = c(var.low,  mean(Posterior.samples.lambda[ix,i]) - nsd*sd(Posterior.samples.lambda[ix,i]))
      var.alpha = c(var.alpha, alpha)
      var.est.l = c(var.est.l, mean(Posterior.samples.lambda[ix,i]))
      var.est.N = c(var.est.N, N.est)
      var.est.a = c(var.est.a, alpha)
      var.est.s = c(var.est.s, sd(Posterior.samples.lambda[ix,i]))
    }

    PosteriorSampleListN[[alphanumber]] = Posterior.sampes.N
    PosteriorSampleListA[[alphanumber]] = Posterior.samples.A
    PosteriorSampleListL[[alphanumber]] = Posterior.samples.lambda

    AcceptInd = 1
    checkupper = mean(Posterior.samples.lambda[ix,1]) + nsd*sd(Posterior.samples.lambda[ix,1])
    checklower = mean(Posterior.samples.lambda[ix,1]) - nsd*sd(Posterior.samples.lambda[ix,1]) 

    if(N.est >= 2){
      #compute upper and lower bound 
      for(i in 2:N.est){
        checkupper = c(checkupper, mean(Posterior.samples.lambda[ix,i]) + nsd*sd(Posterior.samples.lambda[ix,i]))
        checklower = c(checklower, mean(Posterior.samples.lambda[ix,i]) - nsd*sd(Posterior.samples.lambda[ix,i]))
      }
      #Check upper and lower bound of other values 
      for(i in 1:N.est){
        for(j in 1:N.est){
          M = mean(Posterior.samples.lambda[ix,i])
          if(j != i){
            if( M < checkupper[j] & M > checklower[j]){
              AcceptInd = 0
            }
          }
        }
      }
    }

    Accepted[alphanumber] = AcceptInd

    SettingsChar = paste0("GAMMAAB",as.character(gamma.prior.a),"_",as.character(gamma.prior.b),"MCMC",as.character(iter),"BURN",as.character(burnin))
    logfileName = paste0(config$DUMP_LOCATION[[k]],config$NAME[[k]],SettingsChar,'progress.txt')
    fileConn    = file(logfileName)
    progressline = paste("Complete ", as.character(alphanumber/AlphaSize*100), "%\n")
    writeLines(progressline, fileConn)
    close(fileConn)

  } #Done with all alpha's
  #Best estimates of N
  
  

  if(sum(Accepted) == 0){
    AlphaTooLarge = 1
    Bestix = 1
    Best.N = NEST[1]
  }else{
    AlphaTooLarge = 0
    BestCandidate = which(NEST*Accepted == max(NEST*Accepted))
    Bestix = max(BestCandidate)
    Best.N = NEST[Bestix]
  }

  #HISTOGRAM
  maxlim = max(varl)

  #all plotting
  AlphaCollection = NULL
  for(a in 1:AlphaSize){
    print(a)
    New.N = NEST[a]
    alpha = config$ALPHA[[a]][[1]]
    AlphaCollection = c(AlphaCollection, alpha)
    lambda_s = NULL 
    ix = which(PosteriorSampleListN[[a]]== (New.N+1))
    for(l in 1:New.N){
      lambda_s = c(lambda_s, PosteriorSampleListL[[a]][ix,l])
    }
    SettingsChar = paste0("GAMMAAB",as.character(gamma.prior.a),"_",as.character(gamma.prior.b),"ALPHA",as.character(alpha),"MCMC",as.character(iter),"BURN",as.character(burnin),"HIST")
    hist_name = paste0(config$DUMP_LOCATION[[k]],config$NAME[[k]],SettingsChar,'.pdf' )
    pdf(hist_name)
    hist_main = paste('Histogram, Alpha = ', as.character(alpha))
    hist(lambda_s, breaks= 100 , main = hist_main, xlab = "Lambdas", xlim= c(0, maxlim) ,sub = SettingsChar)
    dev.off()

    SettingsChar = paste0("GAMMAAB",as.character(gamma.prior.a),"_",as.character(gamma.prior.b),"ALPHA",as.character(alpha),"MCMC",as.character(iter),"BURN",as.character(burnin),"LIKELIHOOD")
    likelihoodplot_name = paste0(config$DUMP_LOCATION[[k]],config$NAME[[k]],SettingsChar,'.pdf' )
    #Likelihood plots
    pdf(likelihoodplot_name)
    likelihood_name = paste("Likelihood plot","Alpha = ", as.character(alpha))
    likelihood_changes = PosteriorLikelihood[[a]]
    plot(((1:length(likelihood_changes))*100),likelihood_changes,main = likelihood_name, xlab = 'Number of iterations', ylab = "Likelihood", ylim = c(MinLikelihood,MaxLikelihood), sub = SettingsChar)
    dev.off()
  }


  #log N against log alpha plot 
  SettingsChar = paste0("GAMMAAB",as.character(gamma.prior.a),"_",as.character(gamma.prior.b),"ALPHA",as.character(alpha),"MCMC",as.character(iter),"BURN",as.character(burnin),"NVSAL")
  NVSAplot_name = paste0(config$DUMP_LOCATION[[k]],config$NAME[[k]],SettingsChar,'.pdf' )
  pdf(NVSAplot_name)
  NVSA_name = "log(N) Against log(alpha)"
  plot(log10(AlphaCollection),log10(NEST), main = NVSA_name, xlab = 'log(alpha)', ylab = "log(N)", sub = SettingsChar)
  dev.off()
  #JITTERING PLOTS
  #var.est.a
  #var.est.l
  #var.est.N
  #var.est.s
  SettingsChar = paste0("GAMMAAB",as.character(gamma.prior.a),"_",as.character(gamma.prior.b),"MCMC",as.character(iter),"BURN",as.character(burnin),"LVSN")
  LVNplot_name = paste0(config$DUMP_LOCATION[[k]],config$NAME[[k]],SettingsChar,'.pdf' )
  pdf(LVNplot_name)
  LVN_name = "Lambda Estimates Against N"
  jittN = jitter(var.est.N, 1)
  plot(var.est.l ~ jittN, pch = 15, main = LVN_name, xlab = "N", ylab = "Lambdas", sub = SettingsChar)
  arrows(x0=jittN, y0=var.est.l-nsd*var.est.s, x1=jittN, y1=var.est.l+nsd*var.est.s, code=3, angle=90, length=0.1)
  dev.off()

  dis.a = min(abs(diff(log10(var.est.a))))

  SettingsChar = paste0("GAMMAAB",as.character(gamma.prior.a),"_",as.character(gamma.prior.b),"MCMC",as.character(iter),"BURN",as.character(burnin),"LVSLOGA")
  LVNplot_name = paste0(config$DUMP_LOCATION[[k]],config$NAME[[k]],SettingsChar,'.pdf' )
  pdf(LVNplot_name)
  LVN_name = "Lambda Estimates Against LOG(A)"
  jittA = jitter(log10(var.est.a), 0.01*dis.a)
  plot(var.est.l ~ jittA, pch = 15, main = LVN_name, xlab = "log(alpha)", ylab = "Lambdas", sub = SettingsChar)
  arrows(x0=jittA, y0=var.est.l-nsd*var.est.s, x1=jittA, y1=var.est.l+nsd*var.est.s, code=3, angle=90, length=dis.a/2)
  dev.off()

  dis.a = min(abs(diff(var.est.a)))

  SettingsChar = paste0("GAMMAAB",as.character(gamma.prior.a),"_",as.character(gamma.prior.b),"MCMC",as.character(iter),"BURN",as.character(burnin),"LVSA")
  LVNplot_name = paste0(config$DUMP_LOCATION[[k]],config$NAME[[k]],SettingsChar,'.pdf' )
  pdf(LVNplot_name)
  LVN_name = "Lambda Estimates Against A"
  jittA = jitter(var.est.a, 0.01*dis.a)
  plot(var.est.l ~ jittA, pch = 15, main = LVN_name, xlab = "alpha", ylab = "Lambdas", sub = SettingsChar)
  arrows(x0=jittA, y0=var.est.l-nsd*var.est.s, x1=jittA, y1=var.est.l+nsd*var.est.s, code=3, angle=90, length=dis.a/2)
  dev.off()


  SettingsChar = paste0("GAMMAAB",as.character(gamma.prior.a),"_",as.character(gamma.prior.b),"MCMC",as.character(iter),"BURN",as.character(burnin))
  logfileName = paste0(config$DUMP_LOCATION[[k]],config$NAME[[k]],SettingsChar,'.txt')

  ix = which(PosteriorSampleListN[[Bestix]]== (Best.N+1))
  BestA = PosteriorSampleListA[[Bestix]][ix,]
  BestL = PosteriorSampleListL[[Bestix]][ix,]

  fileConn    = file(logfileName)
  SettingLine = paste('We run the following experiment setting:')
  Gammaline   = paste('Our Gamma prior is given as Gamma', as.character(gamma.prior.a),'and', as.character(gamma.prior.b))
  MCMCline    = paste('We run MCMC for',as.character(iter),'iterations, with a burn in of',as.character(burnin))
  Alphaline   = paste('The best alpha', as.character(config$ALPHA[Bestix]))
  Alphaline   = paste(Alphaline,"\n We tried the following alpha \n")
  
  for(a in 1:AlphaSize){
    Alphaline   = paste(Alphaline,as.character( config$ALPHA[[a]][[1]]),"," )
  }

  if(AlphaTooLarge){
    Alphaline   = paste(Alphaline,"\n However, alpha proposal is too large so only the smallest alpha value taken")
  }
  Nline       = paste('The best N is', as.character(Best.N),"(",as.character(mean(PosteriorSampleListN[[Bestix]]== (Best.N+1))),")")
  Aline       = 'The best A estimates are'
  Lline       = 'The best lambda estimates are '
  for(i in 1:Best.N){
    Aline = paste(Aline, as.character(mean(BestA[,i])), "(",as.character(sd(BestA[,i])),")",",")
    Lline = paste(Lline, as.character(mean(BestL[,i])), "(",as.character(sd(BestL[,i])),")",",")
  }

  writeLines(c(SettingLine, Gammaline, MCMCline,Alphaline,Nline,Aline, Lline), fileConn)
  close(fileConn)

  

  CSVFileName = paste0(config$DUMP_LOCATION[[k]],config$NAME[[k]],SettingsChar,'FullEst.csv')
  output.data = cbind(var.name,var.est, var.stat,var.up,var.low,var.alpha)
  output.data = data.frame(output.data)
  output.data$var.est  = as.numeric(output.data$var.est)
  output.data$var.stat = as.numeric(output.data$var.stat)
  output.data$var.up   = as.numeric(output.data$var.up )
  output.data$var.low  = as.numeric(output.data$var.low)
  output.data$var.alpha= as.numeric(output.data$var.alpha)

  write.csv(output.data, CSVFileName)

  BestCSVFileName = paste0(config$DUMP_LOCATION[[k]],config$NAME[[k]],SettingsChar,'BestEst.csv')
  ixforbest = which(output.data$var.alpha==config$ALPHA[Bestix])
  write.csv(output.data[ixforbest,], BestCSVFileName)

  #Full Data 
  Alpha = config$ALPHA
  FullDataName = paste0(config$DUMP_LOCATION[[k]],config$NAME[[k]],SettingsChar,'FullData.RData')
  save(PosteriorSampleListN, PosteriorSampleListA, PosteriorSampleListL,Alpha,Bestix,file = FullDataName)

} 


