library(quantmod)
library(ggplot2)
library(qrmtools)
library(gridExtra)
library(QRM)
library(moments)
library(EnvStats)
library(copula)


### si se quieren los datos en internet ####

# IBEX<-getSymbols("^IBEX", auto.assign = F, from="2004-10-27")
# head(IBEX)
# tail(IBEX)
# 
# DAX<-getSymbols("^GDAXI", auto.assign = F, from="2004-10-27")
# head(DAX)
# tail(DAX)
# 
# X<-na.omit(merge(IBEX$IBEX.Adjusted,DAX$GDAXI.Adjusted))
# head(X)
# log_ret_IBEXyDAX<-diff(log(X))[-1]["2004/2020-04-21"]
# colnames(log_ret_IBEXyDAX)<-c("IBEX", "DAX")
# 
# log_ret_IBEXyDAX<-na.omit(merge(log_ret_IBEXyDAX$IBEX[-which(log_ret_IBEXyDAX$IBEX==0)],log_ret_IBEXyDAX$DAX))
# 
# dim(log_ret_IBEXyDAX)
# 
# anyNA(log_ret_IBEXyDAX)
# 
# saveRDS(log_ret_IBEXyDAX, file = "datosIBEXyDAX.Rda")


log_ret_IBEXyDAX<-readRDS("datosIBEXyDAX.rda")



head(log_ret_IBEXyDAX)
tail(log_ret_IBEXyDAX)
dim(log_ret_IBEXyDAX)
par(mfrow=c(2,1))
plot.xts(log_ret_IBEXyDAX$IBEX, multi.panel =T,
     col=c(1,2), main="a) Log-retornos del IBEX 35",grid.ticks.lwd = 0,grid.ticks.lty = 0,xaxis.same = F, box=T)
plot.xts(log_ret_IBEXyDAX$DAX, multi.panel =T,
         col="red", main="b) Log-retornos del DAX",grid.ticks.lwd = 0,grid.ticks.lty = 0,xaxis.same = F, box=T)




graficoIBEXyDAX<- ggplot(as.data.frame(log_ret_IBEXyDAX),aes(x=IBEX,y=DAX))+geom_point(size=1, shape=16)+xlim(-0.16, 0.16)+ylim(-0.16, 0.16)+theme_bw()+theme(
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),panel.border = element_rect(color="black"))+labs(title.position="bottom",title="a) Log-retornos IBEX 35 y DAX", x="IBEX 35", y="DAX")+theme(plot.title=element_text(size=15, face="bold", hjust = 0.5))+theme(axis.text=element_text(size=8),
                                                                                                                                                                                                                                                                   axis.title=element_text(size=14))
graficoIBEXyDAX

### SIMULACIONES DISTRIBUCI?N NORMAL ####

mu_n <- colMeans(log_ret_IBEXyDAX) #vector de medias
Sigma_n <-var(log_ret_IBEXyDAX)    #matriz varianza-covarianza
A_n <- t(chol(Sigma_n))             #c?lculo del factor Cholesky (funci?n chol())
n <-dim(log_ret_IBEXyDAX)[1]     #tama?o de muestra simulada
d <- dim(log_ret_IBEXyDAX)[2]    #dimensi?n
Z <- matrix(rnorm(n*d), ncol=d) #muestra de simulaciones de i.i.d. N(0,1)
Rsimnor<-rep(mu_n,each=n)+t(A_n%*%t(Z)) #muestra de simulaciones iid N(mu,Sigma)
head(Rsimnor)

graficoIBEXyDAXSimN<- ggplot(as.data.frame(Rsimnor),aes(x=IBEX,y=DAX))+geom_point(size=1, shape=16, col="red")+xlim(-0.16, 0.16)+ylim(-0.16, 0.16)+theme_bw()+theme(
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color="black"))+labs(title.position="bottom",title=expression(paste("b) Log-retornos simulados (","N(",mu,",",Sigma,") ajustada)")), x="IBEX 35", y="DAX")+theme(plot.title=element_text(size=15, face="bold", hjust = 0.5))+theme(axis.text=element_text(size=8),
                                                                                                                                                                                                                                                                         axis.title=element_text(size=14))
graficoIBEXyDAXSimN

grid.arrange(graficoIBEXyDAX,graficoIBEXyDAXSimN, ncol=2, as.table=T)

### SIMULACIONES DISTRIBUCI?N T #####

parametros<-fit.mst(log_ret_IBEXyDAX, method="BFGS")

(Sigma<-as.matrix(parametros$Sigma))           # Par?metro de dispersi?n estimado
(mu<-parametros$mu)                            # Par?metro de  estimado
(df<-parametros$df)                            # Grados de libertad es estimado


A <- t(chol(Sigma))             #calculo del factor Cholesky (funci?n chol())
n <-dim(log_ret_IBEXyDAX)[1]     #tama?o de muestra simulada
d <- dim(log_ret_IBEXyDAX)[2]    #dimensi?n
Z <- matrix(rnorm(n*d), ncol=d) #muestra de simulaciones de i.i.d. N(0,1)
W<-1/rgamma(n,shape=df/2,rate= df/2)  #simulaciones de i.i.d. gamma inversa
Rsimt<-rep(mu,each=n)+sqrt(W)*t(A%*%t(Z))
plot(rbind(coredata(log_ret_IBEXyDAX),Rsimt),col = rep(c("black", "royalblue3"), each = n))
colnames(Rsimt)<-colnames(Rsimnor)

graficoIBEXyDAXSimT<- ggplot(as.data.frame(Rsimt),aes(x=IBEX,y=DAX))+xlim(-0.2, 0.2)+ylim(-0.2, 0.2)+geom_point(size=1, shape=16, col="red")+theme_bw()+theme(
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color="black"))+labs(title.position="bottom",title=expression(paste("b) Log-retornos simulados (","t"[nu],"(",mu,",",Sigma,") ajustada)")), x="IBEX 35", y="DAX")+theme(plot.title=element_text(size=15, face="bold", hjust = 0.5))+theme(axis.text=element_text(size=8),
                                                                                                                                                                                                                                                                                axis.title=element_text(size=14))

graficoIBEXyDAXSimT

grid.arrange(graficoIBEXyDAX+xlim(-0.2, 0.2)+ylim(-0.2, 0.2),graficoIBEXyDAXSimT, ncol=2, as.table=T)

#### Probabilidades valores extremos ####


(pt<-pmvt(lower = -Inf,upper = c(-0.07,-0.07), delta = mu,sigma = Sigma,df=3))[1]

(1-pbinom(0,nrow(log_ret_IBEXyDAX),pt[1]))*100 

(pn<-pmvnorm(lower =-Inf,upper= c(-0.07,-0.07), mean  = mu_n,sigma = Sigma_n))[1]*100

(1-pbinom(0,nrow(log_ret_IBEXyDAX),pn[1]))*100 
(dbinom(5,nrow(log_ret_IBEXyDAX),pn[1]))*100

log_ret_IBEXyDAX[which(log_ret_IBEXyDAX$IBEX <= -0.07 & log_ret_IBEXyDAX$DAX <= -0.07)]

1-pbinom(0,nrow(log_ret_IBEXyDAX),pt[1]) 

### test de multinormalidad ####
maha2_test(log_ret_IBEXyDAX, type = c("ks.test"))

mardia_test(log_ret_IBEXyDAX, type = "skewness")

mardia_test(log_ret_IBEXyDAX, type = "kurtosis")
MardiaTest(Rsimnor)

### histogramas DIST NORM #####

par(mfrow=c(1,2))  


hist(log_ret_IBEXyDAX$IBEX, breaks = 45, probability = T, main="a) Histograma log-retornos IBEX", xlab="Log-retornos IBEX 35", ylim = c(0,47) )

lines(seq(-0.15,0.15,length=100),dnorm(seq(-0.15,0.15,length.out =100),mean = mu_n[1], sd=sqrt(Sigma_n[1,1])), col="red", lwd=2)
lines(density(log_ret_IBEXyDAX$IBEX), col="green", lwd=2)

legend(x=0.01,y=40,legend=c(expression(paste("N(",mu["IBEX"],",",sigma["IBEX"],") ajustada")),
                             "Densidad Kernel"),lty=c(1,1),lwd = c(3,3) ,col=c("red","green"),inset=0.0,pt.lwd = 1,pt.cex = 2,cex = 1,bty="n",seg.len=0.5)


hist(log_ret_IBEXyDAX$DAX, breaks = 45, probability = T, main="b) Histograma log-retornos DAX", xlab="Log-retornos DAX", ylim = c(0,47))
lines(seq(-0.15,0.15,length=100),dnorm(seq(-0.15,0.15,length.out =100),mean = mu_n[2], sd=sqrt(Sigma_n[2,2])), col="red", lwd=2)
lines(density(log_ret_IBEXyDAX$DAX), col="green", lwd=2)

legend(x=0.01,y=40,legend=c(expression(paste("N(",mu["DAX"],",",sigma["DAX"],") ajustada")),
                            "Densidad Kernel"),lty=c(1,1),lwd = c(3,3) ,col=c("red","green"),inset=0.0,pt.lwd = 1,pt.cex = 2,cex = 1,bty="n",seg.len=0.5)

###Normalqqplot ####

par(mfrow=c(1,2))
qqPlot(as.vector(log_ret_IBEXyDAX$IBEX), 
       y = NULL, distribution = "norm",
       ylim = c(-0.16,0.16),xlim = c(-0.05,0.05), 
       param.list = list(mean = mu_n[1], sd = sqrt(Sigma_n[1,1])),
       add.line =T, line.col="red",line.lwd =3, ylab="Cuantiles log-retornos IBEX 35",xlab=expression(paste("Cuantiles teóricos distribución normal ajustada ","N(",mu["IBEX"],",",Sigma["IBEX"],")")),
       main=expression(paste("a) Q-Q Plot: log-ret. IBEX35 vs. distribución normal ajustada ","N(",mu["IBEX"],",",Sigma["IBEX"],")")) )


qqPlot(as.vector(log_ret_IBEXyDAX$DAX), 
       y = NULL, distribution = "norm",
       ylim = c(-0.16,0.16), xlim = c(-0.05,0.05),
       param.list = list(mean = mu_n[2], sd = sqrt(Sigma_n[2,2])),
       add.line =T, line.col="red",line.lwd =3, ylab="Cuantiles log-retornos DAX",xlab=expression(paste("Cuantiles teóricos distribución normal ajustada ","N(",mu["DAX"],",",Sigma["DAX"],")")),
       main=expression(paste("b) Q-Q PLot : log-ret. DAX vs. distribución normal ajustada ","N(",mu["DAX"],",",Sigma["DAX"],")")) )




####### test normalidad de distribuciones marginales ####

ks.test((log_ret_IBEXyDAX$IBEX-mu_n[1])/sqrt(Sigma_n[1,1]),"pnorm")
ks.test((log_ret_IBEXyDAX$DAX-mu_n[2])/sqrt(Sigma_n[2,2]),"pnorm")
ks.test((Rsimnor[,1]-mu[1])/sqrt(Sigma_n[1,1]),"pnorm")

jarque.test(as.vector(coredata(log_ret_IBEXyDAX$IBEX)))
skewness(as.vector(coredata(log_ret_IBEXyDAX$IBEX)))
kurtosis(as.vector(coredata(log_ret_IBEXyDAX$IBEX)), na.rm = FALSE)


jarque.test(as.vector(coredata(log_ret_IBEXyDAX$DAX)))
skewness(as.vector(coredata(log_ret_IBEXyDAX$DAX)))
kurtosis(as.vector(coredata(log_ret_IBEXyDAX$DAX)))
#calculo de coeficientes de asimetria de IBEX y DAX
coef.asimetría <-c(skewness(as.vector(coredata(log_ret_IBEXyDAX$IBEX))),skewness(as.vector(coredata(log_ret_IBEXyDAX$DAX))))
#c?lculo de curtosis de IBEX y DAX 
curtosis<-c(kurtosis(as.vector(coredata(log_ret_IBEXyDAX$IBEX))),kurtosis(as.vector(coredata(log_ret_IBEXyDAX$DAX))))
#Uni?n de los resultados en una matriz
rbind(coef.asimetría,curtosis)


###### Histogrmas con t de Student ####

par(mfrow=c(1,2))  

hist(log_ret_IBEXyDAX$IBEX, breaks = 47, probability = T, main="a) Histograma log-retornos IBEX 35", xlab="Log-retornos IBEX 35", ylim = c(0,47) )

lines(density(log_ret_IBEXyDAX$IBEX), col="green", lwd=3, lty=1)

lines(seq(-0.15,0.15,length=100),dt((seq(-0.15,0.15,length.out =100)-mu[1])/sqrt(Sigma[1,1]), df=df)/sqrt(Sigma[1,1]), col="darkblue", lwd=3, lty=4)

lines(seq(-0.15,0.15,length=100),dnorm(seq(-0.15,0.15,length.out =100),mean = mu_n[1], sd=sqrt(Sigma_n[1,1])), col="red", lwd=2)

legend(x=0.01,y=40,legend=c(expression(paste("N(",mu["IBEX"],",",sigma["IBEX"],") ajustada")),
                             "Densidad Kernel",
                             expression(paste("t"[nu],"(",mu["IBEX"],",",sigma["IBEX"],") ajustada"))),lty=c(1,1,4),lwd = c(3,3,3) ,col=c("red","green","darkblue"),inset=0.0,pt.lwd = 1,pt.cex = 2,cex = 1,bty="n",seg.len=0.4)


hist(log_ret_IBEXyDAX$DAX, breaks = 45, probability = T, main="b) Histograma log-retornos DAX", xlab="Log-retornos DAX", ylim = c(0,47), xlim = c(-0.1,0.13))

lines(density(log_ret_IBEXyDAX$DAX), col="green", lwd=3)

lines(seq(-0.15,0.15,length=100),dnorm(seq(-0.15,0.15,length.out =100),mean = mu_n[2], sd=sqrt(Sigma_n[2,2])), col="red", lwd=2)

lines(seq(-0.15,0.15,length=100),dt((seq(-0.15,0.15,length.out =100)-mu[2])/sqrt(Sigma[2,2]), df=df)/sqrt(Sigma[2,2]), col="darkblue", lwd=3, lty=4)

legend(x=0.01,y=40,legend=c(expression(paste("N(",mu["DAX"],",",sigma["DAX"],") ajustada")),
                            "Densidad Kernel",
                            expression(paste("t"[nu],"(",mu["DAX"],",",sigma["DAX"],") ajustada")) ),lty=c(1,1,4),lwd = c(3,3,3) ,col=c("red","green", "darkblue"),inset=0.0,pt.lwd = 1,pt.cex = 2,cex = 1,bty="n",seg.len=0.4)



### QQplot t ####
par(mfrow=c(1,2))  
qqPlot(as.vector(log_ret_IBEXyDAX$IBEX), 
       y = NULL, distribution = "t",
       param.list = list( df=df),ylim = c(-0.16,0.16),
       add.line =T, line.col="red",line.lwd =3, ylab="Cuantiles log-retornos IBEX",xlab=expression(paste("Cuantiles teóricos distribución t de Student estándar")),
       main=expression(paste("a) Gráfico Q-Q:log-ret. IBEX vs. distribución t de Student estándar" )) )

qqPlot(as.vector(log_ret_IBEXyDAX$DAX), 
       y = NULL, distribution = "t",
       param.list = list( df=df),ylim = c(-0.16,0.16),
       add.line =T, line.col="red",line.lwd =3, ylab="Cuantiles log-retornos DAX",xlab=expression(paste("Cuantiles teóricos distribución t de Student estándar")),
       main=expression(paste("b) Gráfico Q-Q:log-ret. de DAX vs. distribución t de Student estándar")) )

### Kolmogorov-Smirnov distribucion t ####

ks.test((log_ret_IBEXyDAX$IBEX-mu[1])/sqrt(Sigma[1,1]), "pt", df=df)

ks.test((log_ret_IBEXyDAX$DAX-mu[2])/sqrt(Sigma[2,2]), "pt", df=df)

ks.test((Rsimt[,2]-mu[2])/sqrt(Sigma[2,2]), "pt", df=df)

### VaR cartera de ####

TauKendall=cor(log_ret_IBEXyDAX, method = "kendall")
Up=as.data.frame(pobs(log_ret_IBEXyDAX))

ggplot(Up,aes(x=IBEX,y=DAX))+xlim(0, 1)+ylim(0, 1)+geom_point(size=1, shape=16, col="black")+theme_bw()+theme(
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color="black"))+labs(title.position="bottom",title=expression(paste("Pseudo-observaciones de la cópula de los log-retornos del IBEX 35 y del DAX ")))+theme(plot.title=element_text(size=9, face="bold", hjust = 0.5))+theme(axis.text=element_text(size=8),
 axis.title=element_text(size=14))



U=rCopula(3917,tCopula(param = 0.8204, df=3))  
colnames(U)<-colnames(log_ret_IBEXyDAX)
f1=quantile(coredata(log_ret_IBEXyDAX$IBEX), U[,1], names=F)
head(f1)
f2=quantile(coredata(log_ret_IBEXyDAX$DAX), U[,2], names=F)
head(f2)
rsimcop=cbind(f1,f2)
head(rsimcop)
colnames(rsimcop)<-colnames(log_ret_IBEXyDAX)
quantile(exp(0.5*rsimcop[,1]+0.5*rsimcop[,2])-1, 0.01)

hist(exp(0.5*rsimcop[,1]+0.5*rsimcop[,2])-1, ylab="Densidad", xlab = expression(paste("r"["t sim"])),
     axes =F, probability = T, breaks = 40, main = expression(paste("Histograma de ","r"["t sim"])))

axis(1 , c(-0.1,-0.075,quantile(exp(0.5*rsimcop[,1]+0.5*rsimcop[,2])-1, 0.01),0,0.05,0.1), labels =c(-0.1,-0.075, round(quantile(exp(0.5*rsimcop[,1]+0.5*rsimcop[,2])-1, 0.01,names=F),3),0,0.05,0.10) , cex=0.4)

axis(2,seq(0,50,10))
axis(3,tick = F, labels = F)
axis(4,tick = F, labels = F)
lines(density(rsimcop),lwd=2, col="green")
lines(rep(quantile(exp(0.5*rsimcop[,1]+0.5*rsimcop[,2])-1, 0.01),2),c(0,30), lwd=3)
legend(x=quantile(exp(0.5*rsimcop[,1]+0.5*rsimcop[,2])-1, 0.01)-0.025, y =35,cex=1.3,legend=expression(paste("VaR"["1%"])),bty="n", text.width=5)
legend(x=0.04,y=40,legend=c(expression(paste("Densidad Kernel"))),lty=1,lwd=2,col="green",bty="n")

ggplot(as.data.frame(U),aes(x=IBEX,y=DAX))+xlim(0, 1)+ylim(0, 1)+geom_point(size=1, shape=16, col="black")+theme_bw()+theme(
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color="black"))+labs(title.position="bottom",title=expression(paste("Simulaciones de cópula t-Student","", rho,"=0.8204 y ",nu,"=3")))+theme(plot.title=element_text(size=7, face="bold", hjust = 0.5))+theme(axis.text=element_text(size=8), axis.title=element_text(size=14))

a=ggplot(as.data.frame(rsimcop),aes(x=IBEX,y=DAX))+xlim(-0.16, 0.16)+ylim(-0.16, 0.16)+geom_point(size=1, shape=16, col="red")+theme_bw()+theme(
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color="black"))+labs(title.position="bottom",title=expression(paste("b) Simulaciones de log-ret. con cópula ","C"["3, 0.8204"]^"t"," y  marginales", hat(F)[IBEX]," y ", hat(F)[DAX])))+theme(plot.title=element_text(size=12, face="bold"))+theme(axis.text=element_text(size=8), axis.title=element_text(size=14))

grid.arrange(graficoIBEXyDAX+theme(plot.title=element_text(size=16,face = "plain", hjust = 0.5)),a, ncol=2, as.table=T)
a
