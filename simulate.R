rm(list=ls())

# params
N <- 500
mu0 <- 2
mu1 <- 2.5
sd0 <- 0.7
sd1 <- 0.9
lam <- 0.07

# age 
age <- floor(runif(N, 0,90))
ageG <- cut(age, c(-Inf,9,19,29,39,49,59,69,79,Inf))


# simulate
simdf <- data.frame(age=age, age_group=ageG, status=NA, titer=NA)
for(i in 1:nrow(simdf)){
  
  simdf$status[i] <- rbinom(1, 1, prob=1-exp(-lam*simdf$age[i]))
  
  if(simdf$status[i]==0){
    
    simdf$titer[i] <- rnorm(1, mu0, sd0)
    
  }else{
    
    simdf$titer[i] <- rnorm(1, mu0+mu1, sd1)
    
  }
  
  
}

hist(simdf$titer)
ggplot(simdf, aes(titer))+ geom_histogram()
ggplot(simdf, aes(age, titer))+ geom_point()

# save output
output <- list(simdf=simdf[,c("age_group","titer")],
               truepars=data.frame(lam=lam,
                                   mu0=mu0,
                                   mu1=mu1,
                                   sd0=sd0,
                                   sd1=sd1))

setwd("C:/Users/megan/Documents/GitHub/MixCat/data")
saveRDS(output, "SimulatedData.RDS")
