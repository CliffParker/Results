
library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
library(ggplot2)
library(foreach)
library(doParallel)

registerDoParallel()
set.seed(998468235L,kind="L'Ecuyer")

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




"Loading Datasets"
daturl <- "http://kingaa.github.io/pomp/vignettes/twentycities.rda"
datfile <- file.path(tempdir(),"twentycities.rda")
download.file(daturl,destfile=datfile,mode="wb")
load(datfile)

demog$town = factor(demog$town)
measles$town = factor(measles$town)



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

      birthrate=(predict(smooth.spline(x=year,y=births),x=time-4)$y),

      pop=(predict(smooth.spline(x=year,y=pop),x=time)$y) ) -> covar #weekly birth and pop data, birth adjusted for marternal immunity

  merge(tmp,covar, by = "time")->x

  #Converting to Bi-weeks
  x<- as.matrix(x)
  Bdat = Biweekly(x)# Biweekly data

  #here
  result = cbind(time=seq(from = 1944 , by = (14/365.25), length.out = 547),cases=Bdat[,2])
  result = as.data.frame(result)



  Data<- result

}
######################################### For log-Likelihood


dtat$loglik ->x
x=x[!is.na(x)]
max(x)
p=subset(dtat, dtat$loglik == max(x))
p <- unlist(p)
coef(m1)=p




Fit <- traj.match(m1 , start=coef(m1), est = c(),
                  method ="subplex",transform = T)

#log-likelihood assosiated with best parameter
logLik(Fit)


################################################ Residuals



#creating time (starting 4 years before 1944 because of initial adjustments)
tt = time(m1) + 20
tt = tt[tt>1964.99]

ttt = time(m1) - 4
ttt = ttt[ttt<1944]

tt = c(ttt,time(m1),tt)


#plot to make sure its working
plot(trajectory(m1,times= tt,as=T)$H~trajectory(m1,times= tt,as=T)$time,type = "l")
#Data of trajectories
trajectory(m1, times= tt, as = T)[-1,]->z
subset(z, time> 1944 & time<1965)->z

H<- z$H
t<- z$time

#trajectory(m1, params=coef(m1),as=T)$H ->H
#trajectory(m1, params=coef(m1),as=T)$time->t
plot(H~t,type = "l")

rho = p["rho"];rho
psi = p["psi"];psi
casemean= rho*H
#sum of square relative error. relative to data
chi2 = (Data$cases - casemean)/sqrt(casevar); sum(chi2^2)

#Sum of square error
R = (Data$cases[-1] - casemean[-1])
SSR = sum(R^2); SSR


####################################################################################'
####################################################################################'
####################################################################################
                                  #Plots
####################################################################################
#Dataset of trajectoriies
trajectory(m1, times= tt, as = T)[-1,]->zz
subset(zz, time> 1944 & time<1965)->zz


Cases = zz$H * rho

Time =zz$time

Z = cbind(Time, Data, Cases )
Z = as.data.frame(Z)
#############################################################################
ggplot(Z, mapping=aes(x=Time, y=Cases)) +
  geom_line() + ggtitle("Deterministic model with school-term forcing") +
  xlab("Time") + ylab("Cases") + theme_bw()



ZZ = data.frame(Time=Time, Value=Data$cases,  Trajectory = "Data")
ZZZ=  data.frame(Time=Time,  Value=Cases, Trajectory = "Cases")
ZZ = rbind(ZZ,ZZZ)

###################################################################
ggplot(ZZ, mapping=aes(x=Time, y=Value, color=Trajectory)) +
    geom_line() + ggtitle("Deterministic model (school-term forcing)  with data") +
    xlab("Time") + ylab("Value") + theme_bw()





trajectory(m1, times= tt, as = T)[-1,]->ww
subset(ww, time>= 1944 )->ww

Cases = ww$H * rho

Time =ww$time

Z = cbind(Time, Cases )
Z = as.data.frame(Z)

####################################################################################
ggplot(Z, mapping=aes(x=Time, y=Cases)) +
  geom_line() + ggtitle("Deterministic model's (school-term forcing) trajectory into the vaccine era") +
  xlab("Time") + ylab("Cases") + theme_bw()
