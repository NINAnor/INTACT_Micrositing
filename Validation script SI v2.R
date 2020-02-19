#==============================
# INTACT GIS Validation
# Author: R. May
# Date: 09.02.2017
#=============================

library(XLConnect)
library(lme4)
library(MuMIn)
library(effects)
library(lattice)
#library(circular)
#library(stlplus)
library(LMERConvenienceFunctions)
library(AICcmodavg)
library(merTools)
library(lmerTest)
library(cvAUC)

auc.ci <- function(model, data, V=10){
  .cvFolds <- function(Y, V){ #Create CV folds (stratify by outcome)
    Y0 <- split(sample(which(Y==0)), rep(1:V, length=length(which(Y==0))))
    Y1 <- split(sample(which(Y==1)), rep(1:V, length=length(which(Y==1))))
    folds <- vector("list", length=V)
    for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}
    return(folds)
  }
  .doFit <- function(v, folds, data, model){ #Train/test glmer for each fold
    fit <- glmer(formula(model), data=data[-folds[[v]],], family=binomial)
    pred <- predict(fit, newdata=data[folds[[v]],], type="response",allow.new.levels=T)
    return(pred)
  }
  folds <- .cvFolds(Y=data[,"USED"], V=V) #Create folds
  predictions <- unlist(sapply(seq(V), .doFit, folds=folds, data=data, model=model)) #CV train/predict
  predictions[unlist(folds)] <- predictions #Re-order pred values
  # Get CV AUC and confidence interval
  out <- ci.cvAUC(predictions=predictions, labels=data$USED, folds=folds, confidence=0.95)
  return(out)
}

#==========================
#LOAD ALL DATA AND PREPARE DATASET:
setwd("C:/Users/roel.may/OneDrive - NINA/Werk/INTACT/GIS")

data1 <- readWorksheet(loadWorkbook("GPS_RESULTSx.xls"),sheet="GPS_RESULTS",header=TRUE)
Ocols <- c(59,60,62,64,seq(65,85,by=2))
Tcols <- c(61,63,seq(66,86,by=2))
data1$Orographic <- apply(data1[,Ocols],1,max,na.rm=T)
data1$Thermal <- apply(data1[,Tcols],1,max,na.rm=T)
data1 <- data1[,-c(59:86)]
data1 <- data1[,c(1:5,13,15,19:20,57:60)]
data1$USED <- 1

filenames <- "RANDOM_RESULTS"
for(a in 1:5){
  tmp <- readWorksheet(loadWorkbook(paste(filenames[1],"_",a,".xls",sep="")),sheet="RANDOM_RESULTS",header=TRUE)
  tmp$Orographic <- apply(tmp[,(Ocols-51)],1,max,na.rm=T)
  tmp$Thermal <- apply(tmp[,(Tcols-51)],1,max,na.rm=T)
  tmp <- tmp[,-c(1,8:35)]
  tmp$USED <- 0
  tmp <- cbind(tmp[,1],data1[,2],tmp[,c(2:4)],data1[,c(6:9)],tmp[,c(5:9)])
  
  if(a==1){data0 <- tmp}
  if(a>1){data0 <- rbind(data0,tmp)}

}

colnames(data1) <- colnames(data0) <- c("ID","IND","Date_Time","Lat","Lon","Speed","Altitude","Sex","Age","Elevation","Distance","Orographic","Thermal","USED")
data <- rbind(data0,data1)
data <- as.data.frame(data)
data[,3] <- as.POSIXct(data[,3])

rm(data1,data0,Ocols,Tcols,tmp)

data$Orographic[which(is.infinite(data$Orographic))] <- NA
data$Thermal[which(is.infinite(data$Thermal))] <- NA

for(cc in 6:7){
  data[,cc] <- as.numeric(data[,cc])
}

data$Season <- "Winter"
data$Season[which(as.POSIXlt(data$Date_Time)$yday>=73)] <- "Spring"
data$Season[which(as.POSIXlt(data$Date_Time)$yday>=165)] <- "Summer"
data$Season[which(as.POSIXlt(data$Date_Time)$yday>=257)] <- "Autumn"
data$Season <- factor(data$Season,levels=c("Winter","Spring","Summer","Autumn"))

data$Year <- as.POSIXlt(data$Date_Time)$year+1900
data$Year <- as.factor(data$Year)
data$Julian <- as.POSIXlt(data$Date_Time)$yday+1

data$Hour <- as.POSIXlt(data$Date_Time)$hour

data$Updraft <- data$Orographic + data$Thermal

data$Flight <- data$Altitude-data$Elevation

data$OUT <- 0
data$OUT[is.na(data$Speed)] <- 1
data$OUT[which((data$Speed<=0 | data$Speed>=1000))] <- 1
data$OUT[which(data$Distance<100)] <- 1
data$OUT[which(is.na(data$Orographic) | is.na(data$Thermal))] <- 1

# data$Bout <- NA
# dT <- diff(data$Date_Time)
# dT[diff(data$IND!=0)] <- NA
# data <- cbind(data,dT=c(NA,dT))
# interval <- 3600*24
# data$Bout[1] <- 1
# for(rr in 2:nrow(data)){
#   if(data$IND[rr]==data$IND[rr-1] & data$dT[rr]<=interval){
#     data$Bout[rr] <- data$Bout[rr-1]
#   }
#   if(data$IND[rr]!=data$IND[rr-1] | data$dT[rr]>interval){
#     data$Bout[rr] <- data$Bout[rr-1]+1
#   }
# }
# 
# for(rr in 1:nrow(data)){
#   if(data$USED[rr]==1){next}
#   data$Bout[rr] <- data$Bout[which(data$ID==data$ID[rr] & data$USED==1)]
# }

data <- data[which(data$OUT==0),]

data$OUT <- 1
data$OUT[which(data$USED==1)] <- 0
ii <- unique(subset(data,USED==1)$ID)
data$OUT[which(data$USED==0 & data$ID %in% ii)] <- 0

data$Flight2 <- data$Flight
data$Flight2[data$Flight<0 | is.na(data$Flight)] <- 0
data$Flight2[data$Flight>1000] <- 1000

xx<-boxplot(subset(data,OUT==0 & USED==1)$Flight2)
data$FlightF <- "Above RSZ"
data$FlightF[data$Flight2<=40] <- "Below RSZ"
data$FlightF[data$Flight2>40 & data$Flight2<=110] <- "Within RSZ"
#data$FlightF[data$Flight2>xx$stats[3] & data$Flight2<=xx$stats[4]] <- "High"
data$FlightF <- as.factor(data$FlightF)
data$FlightF <- factor(data$FlightF,levels=c("Below RSZ","Within RSZ","Above RSZ"))
data$FlightF <- relevel(data$FlightF,ref="Below RSZ")
data$Flight <- data$FlightF

#Contrasting updraft models:
fit0 <- glmer(USED~1+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit1 <- glmer(USED~Orographic+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit2 <- glmer(USED~Thermal+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit3 <- glmer(USED~Orographic+Thermal+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
aictab(list("null"=fit0,"O"=fit1,"T"=fit2,"O+T"=fit3))

#fit0f <- glmer(USED~1+Flight+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit1f <- glmer(USED~Orographic*Flight+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit2f <- glmer(USED~Thermal*Flight+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit3f <- glmer(USED~(Orographic+Thermal)*Flight+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit4f <- glmer(USED~Orographic+Thermal*Flight+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit5f <- glmer(USED~Orographic*Flight+Thermal+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
aictab(list("OF"=fit1f,"TF"=fit2f,"OF+TF"=fit3f,"O+TF"=fit4f,"OF+T"=fit5f))

#fit0s <- glmer(USED~1+Season+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit1s <- glmer(USED~Orographic*Season+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit2s <- glmer(USED~Thermal*Season+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit3s <- glmer(USED~(Orographic+Thermal)*Season+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit4s <- glmer(USED~Orographic+Thermal*Season+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit5s <- glmer(USED~Orographic*Season+Thermal+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
aictab(list("OS"=fit1s,"TS"=fit2s,"OS+TS"=fit3s,"O+TS"=fit4s,"OS+T"=fit5s))

fit1sf <- glmer(USED~Orographic*Flight+Thermal*Season+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
fit2sf <- glmer(USED~Orographic*Season+Thermal*Flight+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0),family=binomial)
aictab(list("OF+TS"=fit1sf,"OS+TF"=fit2sf))

aictab(list("null"=fit0,"O"=fit1,"T"=fit2,"O+T"=fit3,
            "OF"=fit1f,"TF"=fit2f,"OF+TF"=fit3f,"O+TF"=fit4f,"OF+T"=fit5f,
            "OS"=fit1s,"TS"=fit2s,"OS+TS"=fit3s,"O+TS"=fit4s,"OS+T"=fit5s,
            "OF+TS"=fit1sf,"OS+TF"=fit2sf))

anova(fit3f)
summary(fit3f)

tiff(file="Figure_3.tiff",width=4000,height=2000,units="px",res=400,compression="lzw")
plot(allEffects(fit3f),type="response",ylab="Probability of presence",
     cex.lab=1.5,main="",factor.names=FALSE)
dev.off()

#Assess model performance:
source('kxv_glmer.r')
kxvglmer(formula(fit3f),fit3f@frame,k=10,nbin=10,family=binomial)
(auc <- auc.ci(model=fit3f,data=fit3f@frame,V=10))

# Assess effects of updrafts and season on flight height:
fit.null <- lmer(log(Flight2+0.01)~1+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0 & USED==1))
fitO <- lmer(log(Flight2+0.01)~Orographic+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0 & USED==1))
fitT <- lmer(log(Flight2+0.01)~Thermal+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0 & USED==1))
fitOT <- lmer(log(Flight2+0.01)~Orographic+Thermal+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0 & USED==1))
fitS <- lmer(log(Flight2+0.01)~Season+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0 & USED==1))
fitSO <- lmer(log(Flight2+0.01)~Orographic*Season+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0 & USED==1))
fitST <- lmer(log(Flight2+0.01)~Thermal*Season+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0 & USED==1))
fitSOT <- lmer(log(Flight2+0.01)~(Orographic+Thermal)*Season+(1|IND/Hour)+(1|Year/Season),data=subset(data,OUT==0 & USED==1))
AIC(fit.null,fitO,fitT,fitOT,fitS,fitSO,fitST,fitSOT)
anova(fitOT)
summary(fitS)
plot(allEffects(fitS),type="response",ylab="Flight height (m)",cex.lab=1.5,main="")


#Prepare data to produce figure 2 of the main document:
tiff(file="Figure_2.tiff",width=6000,height=3000,units="px",res=400,compression="lzw")
par(mfrow=c(1,3))
boxplot(subset(data,OUT==0 & USED==1)$Orographic~subset(data,OUT==0 & USED==1)$FlightF,outline=F,xlab="Altitudinal category",ylab="Orographic updraft velocity (m/s)",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1.5))
abline(h=0.75,col="red",lty=2)
mtext("min. sink",col="red",side=2,line=1)
boxplot(subset(data,OUT==0 & USED==1)$Thermal~subset(data,OUT==0 & USED==1)$FlightF,outline=F,xlab="Altitudinal category",ylab="Thermal updraft velocity (m/s)",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1.5))
abline(h=0.75,col="red",lty=2)
mtext("min. sink",col="red",side=2,line=1)
boxplot(subset(data,OUT==0 & USED==1)$Flight2~subset(data,OUT==0 & USED==1)$Season,outline=F,xlab="Season",ylab="Flight altitude (m)",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1000))
polygon(c(0:5, rev(0:5)), c(c(110,110,110,110,110,110), rev(c(40,40,40,40,40,40))),
        col = rgb(green=0,blue=0,red=1,alpha=0.4), border = NA)
mtext("RSZ",col="red",side=2,line=1,at=70)
dev.off()

#Prepare data to produce figure 4 of the main document:
Tdata <- readWorksheet(loadWorkbook("HitraParametre.xlsx"),sheet="Turbinverdier",header=TRUE)[1:24,]
cols <- seq(1:((ncol(Tdata)-3)/2))*2+2
Tdata$maxO <- apply(Tdata[,cols],1,max)
Tdata$maxT <- apply(Tdata[,cols+1],1,max)
Tdata$mnO <- apply(Tdata[,cols],1,mean)
Tdata$mnT <- apply(Tdata[,cols+1],1,mean)

Ndata <- matrix(0,nrow=nrow(Tdata),ncol=ncol(Tdata),byrow=T)
for(cc in 4:31){Ndata[which(Tdata[,cc]>0.75),cc] <- 1}
Tdata$NO <- apply(Ndata[,cols],1,sum)/length(cols)
Tdata$NT <- apply(Ndata[,cols+1],1,sum)/length(cols)
Tdata$Type <- "Turbines 1"
rm(Ndata)

temp <- readWorksheet(loadWorkbook("NyeTurbinerVerdier.xls"),sheet="NyeTurbinerVerdier",header=TRUE)[1:26,-4]
colnames(temp)[1:3] <- colnames(Tdata)[1:3]
cols <- seq(1:((ncol(temp)-3)/2))*2+2
temp$maxO <- apply(temp[,cols],1,max)
temp$maxT <- apply(temp[,cols+1],1,max)
temp$mnO <- apply(temp[,cols],1,mean)
temp$mnT <- apply(temp[,cols+1],1,mean)

Ndata <- matrix(0,nrow=nrow(temp),ncol=ncol(temp),byrow=T)
for(cc in 4:31){Ndata[which(temp[,cc]>0.75),cc] <- 1}
temp$NO <- apply(Ndata[,cols],1,sum)/length(cols)
temp$NT <- apply(Ndata[,cols+1],1,sum)/length(cols)
temp$Type <- "Turbines 2"
rm(Ndata)
Tdata <- rbind(Tdata,temp)

#Collisions at turbines
collisions <- c(24,23,23,3,14,16,15,16,6,24)
temp <- Tdata[collisions,]
temp$Type <- "Collisions"
Tdata <- rbind(Tdata,temp)
Tdata$Crisk <- 0
Tdata$Crisk[which(Tdata$Type=="Turbines 1" & Tdata$Turbines %in% collisions)] <- 1


test <- subset(Tdata,Type=="Turbines 1")
test$NC <- 0
for(c in 1:length(collisions)){test$NC[which(test$Turbines==collisions[c])] <- test$NC[which(test$Turbines==collisions[c])]+1}
test$maxTot <- test$maxT+test$maxO
summary(fit <- glm(Crisk~maxT+maxO,data=test,family=binomial))
summary(fit <- glm(Crisk~mnT+mnO,data=test,family=binomial))
summary(fit <- glm(Crisk~NT+NO,data=test,family=binomial))
summary(fit <- glm(NC~maxT+maxO,data=test,family=poisson))
plot(allEffects(fit),type="response")
wilcox.test(maxT~Crisk,data=test)
wilcox.test(maxO~Crisk,data=test)
wilcox.test(mnT~Crisk,data=test)
wilcox.test(mnO~Crisk,data=test)

tiff(file="Figure_4.tiff",width=5000,height=2500,units="px",res=400,compression="lzw")
par(mfrow=c(1,2))
maxT <- as.data.frame(rbind(cbind(Type=Tdata$Type,Updraft="Orographic",Value=Tdata$maxO),cbind(Type=Tdata$Type,Updraft="Thermal",Value=Tdata$maxT)))
maxT$Value <- as.numeric(as.vector(maxT$Value))
maxT$Type <- factor(maxT$Type,levels=c("Turbines 1","Collisions", "Turbines 2"))
x1<-boxplot(Value~Type+Updraft,data=maxT,xlab="Orographic updraft      Thermal updraft",ylab="Maximum updraft velocity (m/s)",cex.lab=1.5,names=c("Hitra I","Collisions","Hitra II","Hitra I","Collisions","Hitra II"))
abline(h=0.75,lty=2,col="red")

NT <- as.data.frame(rbind(cbind(Type=Tdata$Type,Updraft="Orographic",Value=Tdata$NO),cbind(Type=Tdata$Type,Updraft="Thermal",Value=Tdata$NT)))
NT$Value <- as.numeric(as.vector(NT$Value))
NT$Type <- factor(NT$Type,levels=c("Turbines 1","Collisions","Turbines 2"))
x2<-boxplot(Value~Type+Updraft,data=NT,xlab="Orographic updraft       Thermal updraft",ylab="Proportion over minimum sink",cex.lab=1.5,names=c("Hitra I","Collisions","Hitra II","Hitra I","Collisions","Hitra II"))
dev.off()

fit1 <- lm(Value~Type*Updraft,data=maxT)
fit2 <- glm(cbind(Value*14,14-Value*14)~Type*Updraft,data=NT,family=binomial)

#Produce histogram plots (part of Fig. 1):
oro <- read.table("OROGRAFISK.txt",header=T,sep=",")
ter <- read.table("TERMISK.txt",header=T,sep=",")
h.o <- hist(oro[,3],breaks=seq(0,7,by=0.1),plot=F)
h.o$counts <- h.o$counts/1000
tiff(file="Figure_1o.tiff",width=2000,height=2000,units="px",res=400,compression="lzw")
plot(h.o,xlab="Orografic updraft speed (m/s)",ylab="Frequency (in thousands)",cex.lab=1.5,main="",ylim=c(0,140),xlim=c(0,7))
dev.off()
h.t<- hist(ter[,3],breaks=seq(0,7,by=0.1),plot=F)
h.t$counts <- h.t$counts/1000
tiff(file="Figure_1t.tiff",width=2000,height=2000,units="px",res=400,compression="lzw")
plot(h.t,xlab="Thermal updraft speed (m/s)",ylab="Frequency (in thousands)",cex.lab=1.5,main="",ylim=c(0,2600),xlim=c(0,7))
dev.off()


