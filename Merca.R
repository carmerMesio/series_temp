library(tsfa)
library(nFactors)
#library(sem)
library(vars)
library(readr)
library(foreach)
library(doParallel)
library(lubridate)
library(forecast)

superm<-read_delim("Quarter Hourly Electrical Energy_01-01-2017_31-12-2017.csv", 
                   ";", escape_double = FALSE, trim_ws = TRUE)

##cargamos los festivos en España 2017
load("tabla_festivos.RData")

serie=ts(superm$`General [kWh]`,freq=96)

plot(serie)

param = 310
dif = 364-param

## average by day and by variance.
m=apply(matrix(serie,nr=96),2,mean)
v=apply(matrix(serie,nr=96),2,var)
plot(m,v,xlab="Medias anuales",ylab="Varianzas anuales",main="serie")
abline(lm(v~m),col=2,lty=3,lwd=2)

boxplot(serie~floor(time(serie)))

##USaremos 300 dias para predecir los otros 65.

##Miramos el numero maximo de factores que puede tener la serie

datos=matrix(superm$`General [kWh]`,nrow=param,ncol=96, byrow = T)
predicc=matrix(superm$`General [kWh]`[-c(1:(param*96))],nrow=dif,ncol=96, byrow = T)

##domingos y festivos
valdom=ifelse(weekdays(dmy(superm$Date))=="domingo",1,0)
datdom=matrix(valdom,nrow = param,ncol = 96,byrow = T)[,1]
preddom=matrix(valdom[-c(1:(param*96))],nrow=dif,ncol=96, byrow = T)[,1]

##festivos y fiestas navidad y sem santa
festiusdat=matrix(news$festivo,nrow = param,ncol = 96,byrow = T)[,1]
festiuspred=matrix(news$festivo[-c(1:(param*96))],nrow=dif,ncol=96, byrow = T)[,1]
  
hwfiesdat=matrix(news$holw,nrow=param,ncol=96, byrow = T)[,1]
hwfiespred=matrix(news$holw[-c(1:(param*96))],nrow=dif,ncol=96, byrow = T)[,1]


#numero máximo de factores que puede tener la série.
maxnf=LedermannBound(datos)
#set.seed(12dif)
#idx=sample(1:9,365,replace=T)
#fac=rnorm(365*96,1,0.01)
#dades=ts(datos[idx,]*fac,freq=7)


##Miraremos la correlación al tratarse de analisis factorial.
#eig=eigen(var(diff(datos,lag=1)))
##los eigen dan maximo de 6.
zz <- eigen(cor(datos), symmetric=TRUE)[["values"]] 
print(zz)
par(omi=c(0.1,0.1,0.1,0.1),mar=c(4.1,4.1,0.6,0.1)) 
plot(zz, ylab="Value", xlab="Eigenvalue Number", pch=20:20,cex=1,type="o")

## cogeremos los eigenvalues mayores a 1.
#Contamos los autovalores mayores de 1 con count (metodo antiguo)
count=0
for (i in 1:length(zz)){
  if (zz[i]>1){
    count=count+1
  }
}
## Nos devuelve 6 valores con un eigenvalue mayor a 1.


#Calculamos el numero optimo de autovalores por el Barlett ?QUEHACE?
#Numero optimo de factores es 5 segun Bartlett. Tiene sentido con el plot de eigenvalues.
optnf=nBartlett(cor(datos),N = ncol(datos), alpha=0.05, cor=TRUE, details=TRUE, correction=TRUE)
numfatfit <- FAfitStats(datos, diff.=TRUE, maxit=30000) ##se estima mal nose si faltan iteraciones para hacerlo bien 
print(numfatfit)

nf=optnf$nFactors[1]
##normalize TRUE o FALSE.
factores<-estTSF.ML(datos,nf,normalize=F,maxit=30000)
summary(factores)

##Barplot por factores que se usan en el análisis.
##Barplot factor loadings by day.
par(mfrow=c(2,3))
for (i in 1:nf){
  barplot(factores$loadings[,i])
}

print(summary(1 - factores$stats$uniquenesses))

fmodel<-factors(factores) #modelo dado por los factores TSFA
plot(fmodel)
FAfitStats(factores)
exp=explained(factores, fmodel, names=factores$names) 



seriefactor <- ts(fmodel)
varf1 <- VAR(seriefactor,p=1, type="const")
##selection of variables by AIC.
VARselect(seriefactor,lag.max = 10,type = "const")
VARselect(seriefactor,lag.max = 7,type = "const",exogen = cbind(dom=datdom,s1=fourier_value[1:300,1],c1=fourier_value[1:300,2]))
varf1 <- MTS::VARX(seriefactor,p = 1,xt = datdom)


#MTS::VARXpred(refVARX(varf1),hstep = 1,preddom)
##Usar mas cores del pc.
cl <- makeCluster(detectCores()) #demano els 8 cores del meu pc posa 8 ?
registerDoParallel(cl) #paralelitzo els seguents models el primer tarda i despres els altres els fa rapid.

##Sin valores exogenos en los modelos VAR y AR + Estacionalidad

{
## function bucle
j<-1
datos2 <- datos
datdom2<-datdom
fest2 <- festiusdat
nad2 <- hwfiesdat
p<-vector()
factorspred <- matrix(nrow = nf,ncol = nrow(predicc))
resppred <- matrix(nrow = nrow(predicc),ncol = 96)
resppredarima <- matrix(nrow = nrow(predicc),ncol = 96)
for(i in nrow(datos):363){

factores<-estTSF.ML(datos2,nf,normalize=F,maxit=30000)
fmodel<-factors(factores)
serief1ts <- ts(fmodel)
##Mejor modelo con const.
varf1 <- VAR(serief1ts,p=1, type="const")

##Autoarima with also three covariables auxiliary. And ets with trend of 0.2
prediccarima <- vector()
for(t in 1:ncol(serief1ts)){
vararima <- auto.arima(serief1ts[,t],max.p = 7,max.q = 0,stationary = TRUE)
varets <- as.vector(forecast(ets(serief1ts[,t],model='AAN'),h = 1)[2]$mean)
ab<-predict(vararima)
#prediccarima <- c(prediccarima,ab$pred)
prediccarima <- c(prediccarima,c(.5*ab$pred+.5*varets))
}

resppredarima[j,]<-t(factores$loadings%*%prediccarima)


aa<-predict(varf1,n.ahead = 1)

f1<-aa$fcst$Factor.1[,1]
f2<-aa$fcst$Factor.2[,1]
f3<-aa$fcst$Factor.3[,1]
f4<-aa$fcst$Factor.4[,1]
f5<-aa$fcst$Factor.5[,1]


fac <- t(matrix(c(f1,f2,f3,f4,f5),nrow=length(f1),ncol=5))
factorspred[,j]<-fac
p<-c(p,t(factores$loadings%*%factorspred[,j]))
resppred[j,]<-t(factores$loadings%*%factorspred[,j])
datos2 <- rbind(datos2,predicc[j,])
datdom2 <- c(datdom2,preddom[j])
fest2 <- c(fest2,festiuspred[j])
nad2 <- c(nad2,hwfiespred[j])

par(mfrow=c(2,2))
##plot individual
plot(ts(as.vector(t(datos2[-c(1:(nrow(datos2)-1)),])),freq=96),type="l",ylab='',
     xlab=paste0('Time | RMSE=',round(Metrics::rmse(as.vector(t(datos2[-c(1:(nrow(datos2)-1)),])),
                                                        as.vector(t(resppred[j,]))),2)),
     lty=2,col="blue",
     ylim=c(0,max(c(as.vector(t(datos2[-c(1:(nrow(datos2)-1)),])),as.vector(t(resppred[j,]))))))
title("Individual VARX Model")
lines(ts(as.vector(t(resppred[j,])),freq=96))

##plot colectivo
plot(ts(as.vector(t(datos2[-c(1:nrow(datos)),])),freq=96),type="l",lty=2,col="blue",ylab='')
title("VARX Collective Model")
lines(ts(as.vector(t(resppred)),freq=96))

##plot individual arima
plot(ts(as.vector(t(datos2[-c(1:(nrow(datos2)-1)),])),freq=96),type="l",ylab='',
     xlab=paste0('Time | RMSE=',round(Metrics::rmse(as.vector(t(datos2[-c(1:(nrow(datos2)-1)),])),
                                                    as.vector(t(resppredarima[j,]))),2)),
     lty=2,col="red",ylim=c(0,max(c(as.vector(t(datos2[-c(1:(nrow(datos2)-1)),])),as.vector(t(resppredarima[j,]))))))
title("Individual Autoarima Model")
lines(ts(as.vector(t(resppredarima[j,])),freq=96))

##plot colectivo arima
plot(ts(as.vector(t(datos2[-c(1:nrow(datos)),])),freq=96),type="l",lty=2,col="red",ylab='')
title("Autoarima Model")
lines(ts(as.vector(t(resppredarima)),freq=96))


j<-j+1

}
}
#parado paralelizado.
stopCluster(cl)

## Metrica VARX
Metrics::rmse(as.vector(t(predicc)),as.vector(t(resppred)))

#taula
rmse_v <- numeric(54)
mape_v <- numeric(54)
median_ape_v <- numeric(54)
rmse_vv <- numeric(54)
mape_vv <- numeric(54)
median_ape_vv <- numeric(54)

for (i in 1:dim(resppred)[1]) {

  rmse_v[i] <- Metrics::rmse(predicc[i,],resppred[i,])
  mape_v[i] <- Metrics::mape(predicc[i,],resppred[i,])
  median_ape_v[i] <- median(Metrics::ape(predicc[i,],resppred[i,]))
  
  rmse_vv[i] <- Metrics::rmse(predicc[i,],resppredarima[i,])
  mape_vv[i] <- Metrics::mape(predicc[i,],resppredarima[i,])
  median_ape_vv[i] <- median(Metrics::ape(predicc[i,],resppredarima[i,]))
  
  
}

options(digits=3)
taula1 <- data.frame(VAR_rmse = rmse_v, Ens_rmse = rmse_vv,
                     VAR_mape=mape_v, Ens_mape=mape_vv,
                     VAR_median_ape = median_ape_v, Ens_median_ape = median_ape_vv)



Metrics::mape(as.vector(t(predicc)),as.vector(t(resppred)))
median(ape(as.vector(t(predicc)),as.vector(t(resppred))))

# metric ARIMAX
Metrics::rmse(as.vector(t(resppredarima)),as.vector(t(resppred)))
Metrics::mape(as.vector(t(resppredarima)),as.vector(t(resppred)))
median(ape(as.vector(t(resppredarima)),as.vector(t(resppred))))

par(mfrow=c(2,1))
plot(ts(as.vector(t(datos2[-c(1:nrow(datos)),])),freq=96),type="p",lty=2,col="blue")
lines(ts(as.vector(t(resppred)),freq=96))

ts.plot(ts(as.vector(t(resppred)),freq=96),ts(as.vector(t(datos2[-c(1:nrow(datos)),])),freq=96),lty=c(1,2),col=c(1,4))

summary(aa)

#aa$model$varresult$Factor.1

### Con valores exogenos.
{
  ## function bucle
  j<-1
  datos2 <- datos
  datdom2<-datdom
  fest2 <- festiusdat
  nad2 <- hwfiesdat
  p<-vector()
  factorspred <- matrix(nrow = nf,ncol = nrow(predicc))
  resppred <- matrix(nrow = nrow(predicc),ncol = 96)
  resppredarima <- matrix(nrow = nrow(predicc),ncol = 96)
  for(i in nrow(datos):363){
    
    factores<-estTSF.ML(datos2,nf,normalize=F,maxit=30000)
    fmodel<-factors(factores)
    serief1ts <- ts(fmodel)
    ##Mejor modelo con const.
    varf1 <- VAR(serief1ts,p=1, type="const",exogen = as.matrix(cbind(dat=datdom2,fest=fest2,nad=nad2)))
    
    ##Autoarima with also three covariables auxiliary. And ets with trend of 0.2
    prediccarima <- vector()
    for(t in 1:ncol(serief1ts)){
      vararima <- auto.arima(serief1ts[,t],max.p = 7,max.q = 0,stationary = TRUE,xreg = cbind(dat=datdom2,fest=fest2,nad=nad2))
      varets <- as.vector(forecast(ets(serief1ts[,t],model='AAN'),h = 1)[2]$mean)
      ab<-predict(vararima,newxreg = cbind(preddom[j],festiuspred[j],hwfiespred[j]))
      #prediccarima <- c(prediccarima,ab$pred)
      prediccarima <- c(prediccarima,c(.8*ab$pred+.2*varets))
    }
    
    resppredarima[j,]<-t(factores$loadings%*%prediccarima)
    
    
    aa<-predict(varf1,dumvar = as.matrix(cbind(preddom[j],festiuspred[j],hwfiespred[j])),n.ahead = 1)
    
    f1<-aa$fcst$Factor.1[,1]
    f2<-aa$fcst$Factor.2[,1]
    f3<-aa$fcst$Factor.3[,1]
    f4<-aa$fcst$Factor.4[,1]
    f5<-aa$fcst$Factor.5[,1]
    
    
    fac <- t(matrix(c(f1,f2,f3,f4,f5),nrow=length(f1),ncol=5))
    factorspred[,j]<-fac
    p<-c(p,t(factores$loadings%*%factorspred[,j]))
    resppred[j,]<-t(factores$loadings%*%factorspred[,j])
    print(j)
    datos2 <- rbind(datos2,predicc[j,])
    datdom2 <- c(datdom2,preddom[j])
    fest2 <- c(fest2,festiuspred[j])
    nad2 <- c(nad2,hwfiespred[j])
    
    par(mfrow=c(2,2))
    ##plot individual
    plot(ts(as.vector(t(datos2[-c(1:(nrow(datos2)-1)),])),freq=96),type="l",ylab='',
         xlab=paste0('Time | RMSE=',round(Metrics::rmse(as.vector(t(datos2[-c(1:(nrow(datos2)-1)),])),
                                                        as.vector(t(resppred[j,]))),2)),
         lty=2,col="blue",
         ylim=c(0,max(c(as.vector(t(datos2[-c(1:(nrow(datos2)-1)),])),as.vector(t(resppred[j,]))))))
    title("Individual VARX Model")
    lines(ts(as.vector(t(resppred[j,])),freq=96))
    
    ##plot colectivo
    plot(ts(as.vector(t(datos2[-c(1:nrow(datos)),])),freq=96),type="l",lty=2,col="blue",ylab='')
    title("VARX Collective Model")
    lines(ts(as.vector(t(resppred)),freq=96))
    
    ##plot individual arima
    plot(ts(as.vector(t(datos2[-c(1:(nrow(datos2)-1)),])),freq=96),type="l",ylab='',
         lty=2,col="red",ylim=c(0,max(c(as.vector(t(datos2[-c(1:(nrow(datos2)-1)),])),as.vector(t(resppredarima[j,]))))))
    title("Individual Autoarima Model")
    lines(ts(as.vector(t(resppredarima[j,])),freq=96))
    
    ##plot colectivo arima
    plot(ts(as.vector(t(datos2[-c(1:nrow(datos)),])),freq=96),type="l",lty=2,col="red",ylab='')
    title("Autoarima Model")
    lines(ts(as.vector(t(resppredarima)),freq=96))
    
    
    j<-j+1
    
  }
}

#####ToDo#########
##Coger ahora y mirar efectos calendario. Vacaciones o fines de semana.
##Provar random forest + ARIMA por cada variable (provar mejor ARIMAX). Mirar VarpX.
##################

##Multiplicamos cada elemento de los loadings para saber como sera el consumo del dia 181 en 96 bloques.
data = matrix(0,nrow = 96, ncol= length(f1))

for(i in 1:ncol(factors)){data[,i]=factores$loadings%*%factors[,i]}

Metrics::rmse(datos[301:nrow(datos),],t(data))



