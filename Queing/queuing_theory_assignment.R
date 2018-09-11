install.packages("ggplot2")
library(ggplot2)



########-------- Instantiate Parameters --------########

# Calcualte service time input data
load_factor <- c(0.4,0.7,0.85,0.925) ## loading factors
exp_tau <- 79 ## Expected value of interarrivals
a <- 0.3135 ## Shape Parameter 
exp_service <- c(load_factor*exp_tau) ## Expected value of service time
b <- c(exp_service/gamma((a+1)/a)) ## Scale Parameter 



########-------- Analyse Service times --------########

serviceTimes <- function(loadfactor,seeds = c(7000,4455,9988,3222)){
  
    # Data structures for housing statistics
    sampleStats <- data.frame(matrix(NA, nrow = 1, ncol = 4)) # statistics of samples
    colnames(sampleStats) <- c("Sample mean","Sample variance","Sample median","Sample cv")  
    theoreticalStats <- data.frame(matrix(NA, nrow = 1, ncol = 4)) # statistics of theoretical
    colnames(theoreticalStats) <- c("Thoretical mean","Thoretical variance","Thoretical median","Thoretical cv")  
    compareStates <- data.frame(matrix(NA, nrow = 1, ncol = 4)) # comparison of statistics
    colnames(compareStates) <- c("Compared mean","Compared variance","Compared median","Compared cv")  
  
    # Generate sample service times for varying loading factors
    set.seed(seeds[loadfactor])
    service <- rweibull(10000,shape = a,scale = b[loadfactor]) ## Service times for loading factor specified
    
    # Plot service times with a histogram for visual inspection
    qplot(service,binwidth = 10,main = "Service Times", xlab = "Service Time", ylab = "Count",fill=I("blue"), col=I("black"),alpha=I(.2),xlim=(c(0,400)),ylim=(c(0,1500)))+theme(plot.title = element_text(hjust = 0.5))

    # Calcuate statistics of service times
    meanSample <- mean(service)
    varSample <- var(service)
    stdSample <- sd(service)
    medSample <- median(service)
    cvSample <- stdSample/meanSample
    sampleStats[1,] <- c(meanSample,varSample,medSample,cvSample)
    
    # Calculate the theoretical statistics of service times
    meanTheoretical <- b[loadfactor]*gamma((a+1)/a)
    varTheoretical <- b[loadfactor]*b[loadfactor]*(gamma((a+2)/a)-(gamma((a+1)/a)^2))
    stdTheoretical <- sqrt(varTheoretical)
    medTheoretical <- b[loadfactor]*(log(2)^(1/a))
    cvTheoretical <- stdTheoretical/meanTheoretical
    theoreticalStats[1,] <- c(meanTheoretical,varTheoretical,medTheoretical,cvTheoretical)
    
    # Compare the theoretical and sample results
    meanComp <- meanTheoretical/meanSample
    varComp <- varTheoretical/varSample
    stdComp <- stdTheoretical/stdSample
    medComp <- medTheoretical/medSample
    cvComp <- cvTheoretical/cvSample
    compareStates[1,] <- c(meanComp,varComp,medComp,cvComp)
    
    returnValues <- list(sampleStats,theoreticalStats,compareStates)
    return(returnValues)
}


########-------- G/G/1 Simulation --------########

# **Instructions - Run simulate(loadfactor = n) where loadfactor is the chosen load factor 1=0.4, 2=0.7, 3=0.85, 4=0.925

simulate <- function(loadfactor,n = 100000,seeds = c(2000,2220,3784,2120,2500,2600,2293,2384,2948,3094)){
  
    # Data structures to house simulation data 
    simStats <- data.frame(matrix(NA, nrow = 10, ncol = 4)) # statistics
    colnames(simStats) <- c("L","LQ","W","WQ")  
    confintDF <- data.frame(matrix(NA,nrow = 1,ncol=4)) # confidence intervals
    colnames(confintDF) <- c("Lq lower","Lq upper","Wq lower","Wq upper")
    simStatsMeans <- data.frame(matrix(NA,nrow=1,ncol=4)) # statistics averages
    colnames(simStatsMeans) <- c("mean L", "mean LQ", "mean W","mean WQ" )
    
    
    for(j in 1:10){
          # initialize all variables
          set.seed(seeds[j])
          x <- rweibull(n,a,b[loadfactor]) # service time
          tau <- rgamma(n,scale = 19.75, shape =4) # interarrival time
    
          theta <- numeric(length = n) # exit time instance
          ts <- numeric(length = n) # arrival time instance to the service system
          t <- numeric(length = n) # arrival time instance to the waiting system
          w <- numeric(length = n) # sojourn time in waiting system
          wq <- numeric(length = n) # sojourn time in queue
          lq <- numeric(length = n) # contribution to the occupancy of the queue
          l <- numeric(length = n) # contribution to the occupancy of the waiting system
          lt <- numeric(length = n) # average occupancy at time i
          
          L <- 0 # average occupancy of the waiting system
          LQ <- 0 # Average occupancy in the queue per client
          W <- 0 # average time in the waiting system per client
          WQ <- 0 # average time in the queue per client
    
          # Generate queueing system data
          for(i in 1:n){
                if(i==1){ 
                    ts[i] <- 0
                }else{ ts[i] <- max(theta[i-1],t[i]) } # arrival time intance to the waiting system for client i
                theta[i] <- ts[i] + x[i] # exit time instance for client i from the W.S.
                if(i<n){ t[i+1] <- t[i] + tau[i] } # entrance time instnace to W.S for for client i+1 
            
                # Compute step statistics
                l[i] <- w[i] <- theta[i] - t[i]
                L <- L + l[i]
                lt[i] <- L/(t[i]-t[1])
                W <- W + w[i]
                lq[i] <- wq[i] <- ts[i]-t[i]
                LQ <- LQ + lq[i]
                WQ <- WQ + wq[i]
                
                # Add average occupancy to data frame
                #stable[j,i] <- lt[i]
                #stable[j+10,i] <- t[i]
          }
          
          # Calculate simulation statistics
          W <- W/n
          WQ <- WQ/n
          L <- (L/(t[n]-t[1]))
          LQ <- (LQ/(t[n]-t[1]))
          
          # Add statistics to data frame
          simStats[j,] <- c(L,LQ,W,WQ)
    }
    
    #Calculate the mean of the 10 simulations
    simStatsMeans[1,] <- c(mean(simStats$L),mean(simStats$LQ),mean(simStats$W),mean(simStats$WQ))
    
    #Calculate confidence intervals
    tStat <- abs(qt(.05/2,9)) # 10 simulations with a 5% confidence level
    sdWQ <- sd(simStats$WQ) # sd of WQ
    sdLQ <- sd(simStats$LQ) # sd of LQ
    ciWQUpper <- mean(simStats$WQ)+(sdWQ*tStat)/sqrt(9) # upper ci of WQ
    ciWQLower <- mean(simStats$WQ)-(sdWQ*tStat)/sqrt(9) # lower ci of WQ
    ciLQUpper <- mean(simStats$LQ)+(sdLQ*tStat)/sqrt(9) # upper ci of LQ
    ciLQLower <- mean(simStats$LQ)-(sdLQ*tStat)/sqrt(9) # lower ci of LQ
    confintDF[1,] <- c(ciLQLower,ciLQUpper,ciWQLower,ciWQUpper)
  
  # Plot steady state of 10th simulation
  plot(t,lt, main ="Steady Steady",xlab="Time Instance",ylab="Average Occupancy",col="blue")
  returnValues <- list(simStats,simStatsMeans,confintDF)
  return(returnValues)
}

########-------- Allen Cuneen Formula --------########

# **Instructions - Siply run allenCuneen(). The output will be a data frame with the Wq and Lq for all loading factors

allenCuneen <- function(){
  
      meanTheoretical <- b*gamma((a+1)/a)
      varTheoretical <- b*b*(gamma((a+2)/a)-(gamma((a+1)/a)^2))
      stdTheoretical <- sqrt(varTheoretical)
      medTheoretical <- b*(log(2)^(1/a))

      lambda <- 1/exp_tau # Lambda
      mu <- 1/meanTheoretical # Mu
      theta <- lambda/mu # Theta
      rho <- lambda/mu # Loading factors
      varService <- varTheoretical # variance of service times
      varArrivals <- 4/((4/79)^2) # variance of inter-arrival times
      
      allenC <- data.frame(matrix(NA,nrow = 4,ncol=3)) # data frame to house data
      colnames(allenC) <- c("Load Factor","WQ","LQ")
      
      C <- (theta/(1-rho))/(1+(theta/(1-rho))) # Simplified due to only 1 server
      WQApprox <- C*((lambda*lambda*varArrivals)+(mu*mu*varService))/(2*mu*(1-rho)) # Simplified due to only 1 server
      LQApprox <- WQApprox*lambda
      allenC[,1] <- c("0.4","0.7","0.85","0.925")
      allenC[,2] <- WQApprox
      allenC[,3] <- LQApprox
      return(allenC)
}




