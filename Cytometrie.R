# Chargement du package
library(rflowcyt)
library(mclust)

file1="Dioc6_IP_Ho_08_12_2008 .001 ne1"
file2="AnnexV_IP_12_12_2008.002 new 1"

dioc.data<-read.FCS(file1)@data
anx.data<-read.FCS(file2)@data

ndioc<-dim(dioc.data)[1]
nanx<-dim(anx.data)[1]

{
if (ndioc<nanx)
	{
	index<-sample(1:nanx,ndioc)
	sub1<-cbind.data.frame(dioc.data[,-c(4,6)],AnxV=anx.data[index,3])
	}
else
	{
	index<-sample(1:ndioc,nanx)
	sub1<-cbind.data.frame(dioc.data[index,-c(4,6)],AnxV=anx.data[,3])
	}
}

#sub1[,c(3,5,6)]<-log(sub1[,c(3,5,6)],base=10)
colnames(sub1)<-c("FSC","SSC","Dioc6","Hoechst","IP","Annex.V")
min.point<-c(mean(sub1$FSC),mean(sub1$SSC),0,mean(sub1$Hoechst),0,0)
max.point<-c(mean(sub1$FSC),mean(sub1$SSC),1000,mean(sub1$Hoechst),1000,1000)
sub1<-rbind.data.frame(sub1,min.point,max.point)

write.table(sub1,"data1.xls",sep="\t")

# Représentation des densités
par(mfrow=c(2,3))
for (i in 1:6)
	{
	if (colnames(sub1)[i] %in% c("FSC","SSC","Hoechst")) {lim.x=range(0,1000);lab.x="linear"}
	else {lim.x=range(0,4);lab.x="Log10"}

	plot(density(sub1[,i]),xlim=lim.x,xlab=lab.x,main=colnames(sub1)[i])
	}

# Représentation des données
par(mfrow=c(1,1))
pairs(sub1,col="orangered",pch=3,cex=0.2)

# Analyse du cycle
cycle<-sub1[,4]
cycleBIC<-mclustBIC(cycle,G=1:5)					# Hoegst
cycle.model<-mclustModel(cycle,cycleBIC)
means<-as.numeric(cycle.model$parameters$mean)
prop<-round(cycle.model$parameters$pro,2)
par(mfrow=c(1,1))
mclust1Dplot(cycle, parameters = cycle.model$parameters, z = cycle.model$z, 
             what = "density", identify = TRUE)
for (i in 1:length(means))
	legend(x=means[i],y=0.015/i,legend=prop[i],x.intersp=-1,bty="n",cex=1)



# Sélection sur les critères FSC et SSC
XY<- sub1[,1:2]		# Critères FCS, SSC
library(mclust)

subBIC<-mclustBIC(XY)
mclustXY<-mclustModel(XY,subBIC)
nclass<-dim(mclustXY$parameters$mean)[2]
par(mfrow=c(1,1))
mclust2Dplot(XY, parameters = mclustXY$parameters, z = mclustXY$z,
	 what = "classification", identify = TRUE,
	 symbols=c(1:nclass),colors=c(1:nclass),CEX=0.5)
legend("topleft",legend=paste("Grp",seq(1,nclass)),
	 pch=c(1:nclass),col=c(1:nclass),bty="n")

# Sélection des groupes et pairsplot
grps<-summary(mclustXY)$classification
select<-c(1,3,4)
sub2<-sub1[which(grps %in% select),-c(1,2,4)]			# Sélection des autres variables
sub2<-rbind.data.frame(sub2, rep(1,dim(sub2)[2]))		
pairs(sub2,pch=19,cex=0.5,col="orangered")

# Recherche du modèle de mélanges
subBIC.FL<-mclustBIC(sub2)
mclust.res<-mclustModel(sub2,subBIC.FL)
nclass<-dim(mclust.res$parameters$mean)[2]
grps.FL<-summary(mclust.res)$classification
nclass.FL<-dim(mclust.res$parameters$mean)[2]
pairs(sub2, pch=c(1:nclass.FL)[grps.FL], cex=0.2, col=c(1:nclass.FL)[grps.FL])


# Représentation par ACP
prcomp(sub2)->acp.FL

# 2D
biplot(acp.FL,col=c("white","black"),expand=1.2)
points(I(acp.FL$x[,2]*20)~I(acp.FL$x[,1]*10),
	col=c(1:nclass.FL)[grps.FL], pch=c(1:nclass.FL)[grps.FL], cex=0.5)
legend("topright",legend=paste("Grp",seq(1,nclass.FL)),pch=c(1:nclass.FL),
	col=c(1:nclass.FL),cex=1,bty="n")

# 3D
source("D:\\Stats\\Doc R\\Scripts R\\graph3D.2.R")
graph3D.2(sub2,class1=grps.FL,class2=grps.FL,bold=F)


# Analyse Dioc6/IP à partir de sub1 (changer la variable si sub2)
Dioc.IP<-cbind.data.frame(sub2$Dioc6,sub2$IP)
plot(Dioc.IP)

Dioc.IP.BIC<-mclustBIC(Dioc.IP,G=1:3)
Dioc.IP.model<-mclustModel(Dioc.IP,Dioc.IP.BIC)
prop<-Dioc.IP.model$parameters$pro
means<-Dioc.IP.model$parameters$mean
nclass<-length(prop)

mclust2Dplot(Dioc.IP, parameters = Dioc.IP.model$parameters, z = Dioc.IP.model$z,what = "classification", identify = TRUE)
text(x=means[1,]+100,y=means[2,]+100,labels=paste(round(prop*100,2),"%",sep=""),font=4,cex=1)


