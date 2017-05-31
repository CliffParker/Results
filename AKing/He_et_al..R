
#Packages to install
install.packages("pomp")
install.packages("plyr")
install.packages("reshape2")
install.packages("magrittr")
install.packages("foreach")
install.packages("doParallel")
install.packages("ggplot2")
###############################################################################################################'
##################     Library and things to call                       #######################################'
###############################################################################################################'
###############################################################################################################'
###############################################################################################################'


library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
library(ggplot2)
library(foreach)
library(doParallel)

registerDoParallel()
set.seed(998468235L,kind="L'Ecuyer")
stopifnot(packageVersion("pomp")>="1.4.8")


Biweekly=function(Data){
  n=nrow(Data)/2
  m=ncol(Data)
  mat = matrix( ,n,m)
  for (i in 0:n - 1 ){
    mat[i + 1,]= rep(0, m)
    for (j in 1:2) {
      x = (2*i)+j
      mat[i+1,] = c(mat[i+1,]) + c(Data[x,])
    }
  }
  return(mat)
}

Cumulative = function(Data){
  n=nrow(Data)
  m=ncol(Data)
  Dta=matrix( ,n,m)
  Dta[1,] = Data[1,]
  for (i in 2:n){
    Dta[i,] = Dta[i-1,] + Data[i,]
  }
  
  return(Dta)
}


daturl <- "http://kingaa.github.io/pomp/vignettes/twentycities.rda"
datfile <- file.path(tempdir(),"twentycities.rda")
download.file(daturl,destfile=datfile,mode="wb")
load(datfile)

demog$town = factor(demog$town)
measles$town = factor(measles$town)

"creating City datasets"
for (names in "London") {
  tmp<- subset(demog, town == names)
  tmp<-tmp[,-1]
  tmp %>% subset(year>=1944 & year<1964) %>%
    summarize(
      time=seq(from=min(year),to=max(year),by=1/12),
      pop=predict(smooth.spline(x=year,y=pop),x=time)$y,
      birthrate=predict(smooth.spline(x=year+0.5,y=births),x=time-4)$y
    ) -> covar

  assign( paste0(names,"_covar"),covar)
}


"Cases"
for (names in c("London")) {
  tmp <- subset(measles, town == names)
  tmp %>%
    dcast(date~"cases", fun.aggregate = sum) %>%
    mutate(year=as.integer(format(date,"%Y"))) %>%
    subset(year>=1944 & year<1965) %>%
    mutate(time=(julian(date,origin=as.Date("1944-01-01")))/365.25+1944) %>%
    subset(time>1944 & time<1965, select=c(time,cases)) -> tmp #Weekly incidence data

  ################## Removing the column with town names
  covar<- subset(demog, town == names)
  covar<-covar[,-1]
  ##################

  covar %>% subset(year>=1944 & year<1964) %>%
    summarize(
      time= tmp$time,

      birthrate=(predict(smooth.spline(x=year,y=births),x=time-4)$y)/52,

      pop=(predict(smooth.spline(x=year,y=pop),x=time)$y)/52 ) -> covar #weekly birth and pop data, birth adjusted for marternal immunity

  merge(tmp,covar, by = "time")->x

  #Converting to Bi-weeks
  x<- as.matrix(x)
  Bdat = Biweekly(x)# Biweekly data

  #here
  result = cbind(time=seq(from = 1944 , by = (14/365.25), length.out = 547),cases=Bdat[,2])
  result = as.data.frame(result)



  assign( paste0(names,"_BiData"),result)

}




###############################################################################################################'
###############################################################################################################'
########################       Creating pomp Body parts      ##################################################'
###############################################################################################################'
###############################################################################################################'

rproc <- Csnippet("
                  double beta, br, seas, va, foi, dw, births;
                  double rate[6], trans[6];

                  //Vacination uptake
                  if ( t< 1968)
                  va = 0;
                  else if (t>=1968 && t<=1978)
                  va = 0.4  + 0.4 * (t-1968)/10;
                  else
                  va = 0.8;

                  // cohort effect
                  if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
                  br = cohort*birthrate/dt + (1-cohort)*birthrate;
                  else
                  br = (1.0-cohort)*birthrate;

                  // term-time seasonality
                  t = (t-floor(t))*365.25;
                  if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
                  seas = 1.0+amplitude*0.2411/0.7589;
                  else
                  seas = 1.0-amplitude;

                  // transmission rate
                  beta = R0*(gamma+mu)*seas;
                  // expected force of infection
                  foi = beta*pow(I+iota,alpha)/pop;
                  // white noise (extrademographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);

                  rate[0] = foi*dw/dt;  // stochastic force of infection
                  rate[1] = mu;             // natural S death
                  rate[2] = sigma;        // rate of ending of latent stage
                  rate[3] = mu;             // natural E death
                  rate[4] = gamma;        // recovery
                  rate[5] = mu;             // natural I death

                  // Poisson births
                  births = rpois(br*(1-va)*dt);

                  // transitions between classes
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,E,&rate[2],dt,&trans[2]);
                  reulermultinom(2,I,&rate[4],dt,&trans[4]);

                  S += births   - trans[0] - trans[1];
                  E += trans[0] - trans[2] - trans[3];
                  I += trans[2] - trans[4] - trans[5];
                  R = pop - S - E - I;
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[4];           // true incidence
                  ")

initlz <- Csnippet("
                   double m = pop/(S_0+E_0+I_0+R_0);
                   S = nearbyint(m*S_0);
                   E = nearbyint(m*E_0);
                   I = nearbyint(m*I_0);
                   R = nearbyint(m*R_0);
                   W = 0;
                   C = 0;
                   ")
dmeas <- Csnippet("
                  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  ")
rmeas <- Csnippet("
                  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")


toEst <- Csnippet("
                  Tmu = log(mu);
                  Tsigma = log(sigma);
                  Tgamma = log(gamma);
                  Talpha = log(alpha);
                  Tiota = log(iota);
                  Trho = logit(rho);
                  Tcohort = logit(cohort);
                  Tamplitude = logit(amplitude);
                  TsigmaSE = log(sigmaSE);
                  Tpsi = log(psi);
                  TR0 = log(R0);
                  to_log_barycentric (&TS_0, &S_0, 4);
                  ")

fromEst <- Csnippet("
                    Tmu = exp(mu);
                    Tsigma = exp(sigma);
                    Tgamma = exp(gamma);
                    Talpha = exp(alpha);
                    Tiota = exp(iota);
                    Trho = expit(rho);
                    Tcohort = expit(cohort);
                    Tamplitude = expit(amplitude);
                    TsigmaSE = exp(sigmaSE);
                    Tpsi = exp(psi);
                    TR0 = exp(R0);
                    from_log_barycentric (&TS_0, &S_0, 4);
                    ")#



################################################################################################'
#################################################################################################'
#################################################################################################'
##### Main ######################################################################################'
######################################## Object #################################################'
#################################################################################################'
########################### to ##################################################################'
#################################################################################################'
#################################################################################################'
#################################################################### run ########################'
#################################################################################################'
#################################################################################################'

  name = "London"
    #Assigning to parest the right name associated parameters
    #######################################################################'
    #                            POMP BODY
    #######################################################################'
    get(paste0(name,"_BiData")) %>%
      pomp(t0=1940,
           time="time",
           rprocess=euler.sim(rproc,delta.t=1/365.25),
           initializer=initlz,
           dmeasure=dmeas,
           rmeasure=rmeas,
           covar=get(paste0(name,"_covar")),
           tcovar="time",
           zeronames=c("C","W"),
           toEstimationScale=toEst,
           fromEstimationScale=fromEst,
           statenames=c("S","E","I","R","C","W"),
           paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                        "rho","sigmaSE","psi","cohort","amplitude",
                        "S_0","E_0","I_0","R_0")) -> m1


         coef(m1)<- coef(secondFit) #secondFit has the best parameter estimates



          secondFit -> mf

        ## Runs 10 particle filters to assess Monte Carlo error in likelihood
        ###################################################
                          #logLik
        ###################################################'
        foreach(i=1:5,
                .packages="pomp",
                .options.multicore=list(set.seed=TRUE)
        ) %dopar% {
          pfilter(mf, Np = 5000)
        } -> pf

        ll <- sapply(pf,logLik)
        ll <- logmeanexp(ll, se = TRUE)
        (ll) #logLik here!!!




        ####################################################################################
                                          #residuals
        ####################################################################################
        tt = time(m1) + 20
        tt = tt[tt>1964.99]
        
        ttt = time(m1) - 4
        ttt = ttt[ttt<1944]
        
        tt = c(ttt,time(m1),tt)
        
        
        
        set.seed(998468235L,kind="L'Ecuyer")
        int = function(x,k){
          x = as.double(x)
          th = (10^k)
          x= x * th
          x = as.integer(x)
          (x/th)
        }

        f = function(){
         
          #creating dataset by simulating 500 times
          m1%>%
            simulate(params=coef(secondFit),nsim=200,as.data.frame=TRUE , times= tt ,include.data=TRUE) %>%
            subset(time>=1944 & time< 1965,select=c(time,sim,cases)) %>%
            mutate(data=sim=="data") ->dta
          
          #Creating variance dataset
          dta %>%
            subset(sim!="data") %>%
            ddply( .(time), summarize, sim = "Var", cases=var(cases), data = "Var") -> VAR
          VAR[1:547,]->VAR
          
          
          #creating mean dataset
          dta %>%
            subset(sim!="data") %>%
            ddply( .(time), summarize, sim = "mean", cases=mean(cases), data = "mean") -> MEAN
          MEAN[1:547,]->MEAN
          
          
          #creating cases dataset
          subset(dta,data == TRUE)->DATA
          
          #creating an extened dataset with cases, the mean and variance 
          rbind(dta,VAR,MEAN) ->NEW
          #subsetting dataset to data cases, mean and variance 
          subset(NEW,data!= FALSE)->NEW
          
          #Data managemenet
          NEW$sim = factor(NEW$sim)
          NEW$data = factor(NEW$data)
          
          
          #creating wide dataset
          dcast(NEW,time~sim,value.var = "cases") -> NEW
          
          # 
          PR = (NEW$data-NEW$mean ) / NEW$data        
          chi2=  sum(PR^2); chi2
        }
        
        
        
        ll=replicate(30,f())
        mean(ll); sqrt(var(ll))
        
        


        ######################################################################################################'
        #############################          plotting       ################################################'
        ######################################################################################################'
        ######################################################################################################'
        ######################################################################################################'
        # time extention into vaccine era
        tt = time(m1) + 20
        tt = tt[tt>1964.99]
        tt = c(time(m1),tt)

        # Model variability
        m1 %>%
          simulate(params=coef(m1),nsim=200,as.data.frame=TRUE , times= tt ,include.data=TRUE) %>%
          subset(select=c(time,sim,cases)) %>%
          mutate(data=sim=="data") %>%
          ddply(~time+data,summarize,
                p=c(0.05,0.5,0.95),q=quantile(cases,prob=p,names=FALSE)) %>%
          mutate(p=mapvalues(p,from=c(0.05,0.50,0.95),to=c("lo","med","hi")),
                 data=mapvalues(data,from=c(TRUE,FALSE),to=c("data","simulation"))) %>%
          dcast(time+data~p,value.var='q') %>%
          ggplot(aes(x=time,y=med,color=data,fill=data,ymin=lo,ymax=hi))+
          geom_ribbon(alpha=0.2)+ ylab("cases") + ggtitle("Variability in the He et al. model")+theme_bw()

        # Model simulations into the vaccine era
        m1 %>%
          simulate(params=coef(m1),nsim=6,as.data.frame=TRUE, times= tt ,include.data=F) %>%
          ggplot(aes(x=time,y=cases,group=sim,color=(sim=="data")))+
          guides(color=FALSE)+
          geom_line()+facet_wrap(~sim,ncol=2)+ ggtitle("He et al. model simulations into vaccine era")+ theme_bw()


        #"Mean Plotting"
        m1%>%
          simulate(params=coef(firstFit),nsim=200,as.data.frame=TRUE ,include.data=TRUE) %>%
          subset(select=c(time,sim,cases)) %>%
          mutate(data=sim=="data") ->dta

        dta %>%
          subset(sim!="data") %>%
          ddply( .(time), summarize, sim = "mean", cases=mean(cases), data = "mean") %>%
          rbind(dta)->dta

        ggplot(subset(dta,data !=FALSE), mapping=aes(x=time, y=cases, color=data)) +
          geom_line() + ggtitle("He et al. average cases with data") +
          xlab("time") + ylab("cases") + theme_bw()
