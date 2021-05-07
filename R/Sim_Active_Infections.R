#' Simulate Active Infections
#'
#' @param dat1 Standard infection pattern (vector) for black barplot
#' @param dat2 Comparison infection pattern (vector) for blue barplot
#' @param dat3 Comparison infectin pattern (vector) for red barplot
#' @param MU1 expected time  (integer) to getting infected for 1st column
#' @param MU2 expected time (integer) to getting infected for 2nd column
#' @param size1 Dispersion parameter (integer) for negative binomial distribution so that var = MU + MU^2/size (first column)
#' @param size2 Dispersion parameter (integer) for negative binomial distribution so that var = MU + MU^2/size (second column)
#' @param parallel Option to parallelize simulation replicates
#' @param max_iter number if iterations for replicates
#'
#'
#' @return 2 x 3 barplots comparing daily active infections for different infection rates and different parameters (MU and size).
#'
#' @importFrom stats rnbinom
#' @importFrom stats rpois
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import fame
#'
#' @export
#'
#' @examples
#' Sim_Active_Infections(dat1 = StandardData,dat2 = EarlyIntervention,dat3 = LowerRate,MU1=4,MU2 = 3.6,size1=1,size2=0.9,parallel = FALSE,max_iter=5 )
Sim_Active_Infections = function(dat1 = StandardData,dat2 = EarlyIntervention,dat3 = LowerRate,MU1=4,MU2 = 3.6,size1=1,size2=0.9,parallel = FALSE,max_iter=5 ){

  Simul = function(days=300, period=30, rt=rr, muT=MU1, size=size1,limit= 1000000, pi=0.001,n0=1){
    #days = # of days in simulation
    #nd = the simulation period
    #rt = infection rate
    #muT = mean time transmission time for individual
    #size = dispersion parameter so var = mu + mu^2/size.
    #limit = study population size
    #pi = proportion of people with immunity
    #n0 = number of initial infectious people

    daily_new = rep(0,days) #daily new cases
    active_day = daily_new #number of active cases each day
    total_cases = 0
    nn = length(daily_new)

    if (period > length(rt)) {
      stop("length of rt should not be smaller than the period.")
    }
    stoplimit = limit*(1-pi)


    if (max_iter <= 1){
      stop("not enough iterations")
    }

    nk = n0
    for (i in 1:nk){
      if (total_cases > stoplimit) {
        rt[1] = 0.001
      }
      ni = rpois(1,rt[1]) #tells us number of people infected by virus carrier
      immune = sample(c(0,1),1,prob=c((1-pi),pi))
      if (immune == 1){
        ni = 0
      }
      total_cases = total_cases + ni
      if (ni > 0){
        nthdaycase= rep(0,ni)
        for (j in 1:ni){
          nthdaycase[j] = rnbinom(1,size=size,mu=muT)+ 1 #nth day on which new case occurs
          daily_new[nthdaycase[j]] = daily_new[nthdaycase[j]]+1
        }


        past = c(rep(1,(max(nthdaycase)-1)),rep(0,(days-max(nthdaycase)+1)))
        active_day = active_day+past
      }
    }

    for (i in 2:period){
      nk = daily_new[i-1] # num people newly infected on (i-1)th day
      if (nk > 0){
        for (k in 1:nk){
          if (total_cases > stoplimit){
            rt[i] = 0.001
          }
          ni = rpois(1,rt[i])
          immune = sample(c(0,1),1,prob=c(1-pi,pi))
          if (immune == 1){
            ni = 0
          }
          total_cases = total_cases + ni
          if (ni > 0){
            nthdaycase = rep(0,ni)
            for (j in 1:ni){
              nthdaycase[j] = rnbinom(1,size=size,mu=muT)+1+i
              daily_new[nthdaycase[j]] = daily_new[nthdaycase[j]] + 1
            }
            past = c(rep(0,(i-1)),rep(1,(max(nthdaycase)+1-i)),rep(0,(days-max(nthdaycase))))
            active_day = active_day+past

          }
        }
      }
    }

    list(riskpop = active_day, dailynew = daily_new,total=total_cases)



  }

  if(parallel == TRUE){

    allseeds = seq(1,10000,1)
    seeds = sample(allseeds,max_iter)

    cores=detectCores()
    cl <- makeCluster(7,type="SOCK")#not to overload your computer
    registerDoParallel(cl)

    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {

      run1.11 = (Simul(period=100, muT=MU1,rt = dat1,size=size1))$riskpop
      run1.11
    }


    run1.11 = rowMeans(final)

    stopCluster(cl)


    cores=detectCores()
    cl <- makeCluster(cores-1,type="SOCK")#not to overload your computer
    registerDoParallel(cl)

    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {
      set.seed(seeds[i])
      run1.12 = (Simul(period=100, muT=MU1,rt = dat2,size=size1))$riskpop
      run1.12
    }



    run1.12 = rowMeans(final)

    stopCluster(cl)

    cores=detectCores()
    cl <- makeCluster(cores-1,type="SOCK")#not to overload your computer
    registerDoParallel(cl)



    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {
      set.seed(seeds[i])
      run1.21 = (Simul(period=100, muT=MU2,rt = dat1,size=size1))$riskpop
      run1.21
    }


    run1.21 = rowMeans(final)

    stopCluster(cl)


    cores=detectCores()
    cl <- makeCluster(cores-1,type="SOCK")#not to overload your computer
    registerDoParallel(cl)

    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {
      set.seed(seeds[i])
      run1.22 = (Simul(period=100, rt=dat2,muT=MU2,size=size1))$riskpop
      run1.22
    }


    run1.22 = rowMeans(final)
    stopCluster(cl)


    cores=detectCores()
    cl <- makeCluster(cores-1,type="SOCK")#not to overload your computer
    registerDoParallel(cl)

    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {
      set.seed(seeds[i])
      run1.31 = (Simul(period=100,rt=dat1, muT=MU1,size=size2))$riskpop
      run1.31
    }


    run1.31 = rowMeans(final)
    stopCluster(cl)


    cores=detectCores()
    cl <- makeCluster(cores-1,type="SOCK")#not to overload your computer
    registerDoParallel(cl)

    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {
      set.seed(seeds[i])
      run1.32 = (Simul(period=100,rt=dat2,muT=MU1,size=size2))$riskpop
      run1.32
    }


    run1.32 = rowMeans(final)
    stopCluster(cl)


    cores=detectCores()
    cl <- makeCluster(cores-1,type="SOCK")#not to overload your computer
    registerDoParallel(cl)

    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {
      set.seed(seeds[i])
      run2.11 = (Simul(period=100,rt=dat1 ,muT=MU1,size=size1))$riskpop
      run2.11
    }


    run2.11 = rowMeans(final)

    stopCluster(cl)


    cores=detectCores()
    cl <- makeCluster(cores-1,type="SOCK")#not to overload your computer
    registerDoParallel(cl)

    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {
      set.seed(seeds[i])
      run2.12 = (Simul(period=100, rt=dat3,muT=MU1,size=size1))$riskpop
      run2.12
    }


    run2.12 = rowMeans(final)

    stopCluster(cl)


    cores=detectCores()
    cl <- makeCluster(cores-1,type="SOCK")#not to overload your computer
    registerDoParallel(cl)

    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {
      set.seed(seeds[i])
      run2.21 = (Simul(period=100, rt=dat1,muT=MU2,size=size1))$riskpop
      run2.21
    }



    run2.21 = rowMeans(final)
    stopCluster(cl)


    cores=detectCores()
    cl <- makeCluster(cores-1,type="SOCK")#not to overload your computer
    registerDoParallel(cl)


    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {
      set.seed(seeds[i])
      run2.22 = (Simul(period=100, rt=dat3,muT=MU2,size=size1))$riskpop
      run2.22
    }


    run2.22 = rowMeans(final)
    stopCluster(cl)



    cores=detectCores()
    cl <- makeCluster(cores-1,type="SOCK")#not to overload your computer
    registerDoParallel(cl)

    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {
      set.seed(seeds[i])
      run2.31 = (Simul(period=100, rt=dat1,muT=MU1,size=size2))$riskpop
      run2.31
    }


    run2.31 = rowMeans(final)
    stopCluster(cl)


    cores=detectCores()
    cl <- makeCluster(cores-1,type="SOCK")#not to overload your computer
    registerDoParallel(cl)


    final <- foreach(i=1:max_iter, .combine=cbind) %dopar% {
      set.seed(seeds[i])
      run2.32 = (Simul(period=100, rt=dat3,muT=MU1,size=size2))$riskpop
      run2.32
    }


    run2.32 = rowMeans(final)
    stopCluster(cl)



    par(mfrow = c(2, 3))

    barplot(run1.11,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab = "Active Infections")
    barplot(run1.12,xlim=c(0,79),ylim=c(0,5000),col='blue', add=T)

    barplot(run1.21,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab = "Active Infections")
    barplot(run1.22,xlim=c(0,79), ylim=c(0,5000),col='blue', add=T)

    barplot(run1.31,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab = "Active Infections")
    barplot(run1.32,xlim=c(0,79),ylim=c(0,5000), col='blue', add=T)

    barplot(run2.11,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab ="Active Infections")
    barplot(run2.12,xlim=c(0,79),ylim=c(0,5000),col='red', add=T)

    barplot(run2.21,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab = "Active Infections")
    barplot(run2.22,xlim=c(0,79),ylim=c(0,5000), col='red', add=T)

    barplot(run2.31,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab = "Active Infections")
    barplot(run2.32,xlim=c(0,79), ylim=c(0,5000),col='red', add=T)

  }else{

    allseeds = seq(1,10000,1)
    seeds = sample(allseeds,max_iter)

    run1.11 = list()
    run1.12 = list()
    run1.21 = list()
    run1.22 = list()
    run1.31 = list()
    run1.32 = list()
    run2.11 = list()
    run2.12 = list()
    run2.21 = list()
    run2.22 = list()
    run2.31 = list()
    run2.32 = list()


    for (it in 1:max_iter){



      rr = dat1
      set.seed(seeds[it])
      run1.11[[it]] = (Simul(period=100, muT=MU1,size=size1))$riskpop
      rr=dat2
      set.seed(seeds[it])
      run1.12[[it]] = (Simul(period=100, muT=MU1,size=size1))$riskpop


      rr=dat1
      set.seed(seeds[it])
      run1.21[[it]] = (Simul(period=100, muT=MU2,size=size1))$riskpop
      rr=dat2
      set.seed(seeds[it])
      run1.22[[it]] = (Simul(period=100, muT=MU2,size=size1))$riskpop


      rr=dat1
      set.seed(seeds[it])
      run1.31[[it]] = (Simul(period=100, muT=MU1,size=size2))$riskpop
      rr=dat2
      set.seed(seeds[it])
      run1.32[[it]] =  (Simul(period=100, muT=MU1,size=size2))$riskpop

      rr=dat1
      set.seed(seeds[it])
      run2.11[[it]] = (Simul(period=100, muT=MU1,size=size1))$riskpop
      rr = dat3
      set.seed(seeds[it])
      run2.12[[it]] = (Simul(period=100, muT=MU1,size=size1))$riskpop

      rr=dat1
      set.seed(seeds[it])
      run2.21[[it]] = (Simul(period=100, muT=MU2,size=size1))$riskpop
      rr = dat3
      set.seed(seeds[it])
      run2.22[[it]] = (Simul(period=100, muT=MU2,size=size1))$riskpop

      rr=dat1
      set.seed(seeds[it])
      run2.31[[it]] =  (Simul(period=100, muT=MU1,size=size2))$riskpop
      rr = dat3
      set.seed(seeds[it])
      run2.32[[it]] = (Simul(period=100, muT=MU1,size=size2))$riskpop
    }

    run1.11 = Reduce("+", run1.11) / length(run1.11)
    run1.12 =  Reduce("+", run1.12) / length(run1.12)
    run1.21 =  Reduce("+", run1.21) / length(run1.21)
    run1.22 =  Reduce("+", run1.22) / length(run1.22)
    run1.31 =  Reduce("+", run1.31) / length(run1.31)
    run1.32 =  Reduce("+", run1.32) / length(run1.32)
    run2.11 =  Reduce("+", run2.11) / length(run2.11)
    run2.12 =  Reduce("+", run2.12) / length(run2.12)
    run2.21 = Reduce("+", run2.21) / length(run2.21)
    run2.22 =  Reduce("+", run2.22) / length(run2.22)
    run2.31 =  Reduce("+", run2.31) / length(run2.31)
    run2.32 =  Reduce("+", run2.32) / length(run2.32)

    par(mfrow = c(2, 3))

    barplot(run1.11,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab = "Hospitalizations")
    barplot(run1.12,xlim=c(0,79),ylim=c(0,5000),col='blue', add=T)

    barplot(run1.21,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab = "Hospitalizations")
    barplot(run1.22,xlim=c(0,79), ylim=c(0,5000),col='blue', add=T)

    barplot(run1.31,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab = "Hospitalizations")
    barplot(run1.32,xlim=c(0,79),ylim=c(0,5000), col='blue', add=T)

    barplot(run2.11,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab = "Hospitalizations")
    barplot(run2.12,xlim=c(0,79),ylim=c(0,5000),col='red', add=T)

    barplot(run2.21,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab = "Hospitalizations")
    barplot(run2.22,xlim=c(0,79),ylim=c(0,5000), col='red', add=T)

    barplot(run2.31,xlim=c(0,79),ylim=c(0,5000),xlab="Days",ylab = "Hospitalizations")
    barplot(run2.32,xlim=c(0,79), ylim=c(0,5000),col='red', add=T)



  }

}

