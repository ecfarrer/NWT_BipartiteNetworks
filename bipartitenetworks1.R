setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Tulane/Lab/GradStudents/MonicaBrady/ITS_TAXA_RAREFIED_5238")
library(tidyr)
library(vegan)
library(dplyr)
library(bipartite)
library(igraph)
library(Hmisc)

#save.image("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Tulane/Lab/GradStudents/MonicaBrady/workspace1.Rdata")  # 
load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Tulane/Lab/GradStudents/MonicaBrady/workspace1.Rdata")  # 


#undrstanding how stucture (especially specialization) of networks changes acoss space is interesting and useul to management. If one community has high specialization, consequences of biodiveristy loss are different.


#### EMILY code #####
otu_table_rare_L6T<-read.csv("otu_table_rare_L6T.csv",  sep=",", header=TRUE,comment.char = "#")
genus <-otu_table_rare_L6T%>%
  separate(col=PLOT.ID, into=c("Plant.Species","Plot.Number"), sep="\\.")%>%
  arrange(Plant.Species)%>%
  filter(Plant.Species!="UNKOPP",Plant.Species!="CARSPP",Plant.Species!="POASPP")

#get rid of singletons (taxa only appearing in one sample,prior to taking the mean), doubletons, or rare taxa (taxa with asummed relative abudnance of <0.6%). Not sure if you want to do this (maybe not for the stats but more filtering for the figures?)
#ind<-which(colSums(genus[,3:384]>0)>1)
ind<-which((colSums(genus[,3:384]>0)>2)&colSums(genus[,3:384])>0.006)
length(ind)
genus2<-cbind(genus[,1:2],genus[,ind+2])

rowSums(genus2[3:170]>0)

genus3<-genus2%>%
  select(-Plot.Number)%>%
  group_by(Plant.Species)%>%
  summarise_all(mean)
genus3<-as.data.frame(genus3)
rownames(genus3)<-as.character(genus3[,1])
genus4 <- as.matrix(genus3[,-1])

#delete rare taxa at <0.001 abundance, I did this b/c I realized that plants were associating with 58 microbial taxa out of a total of 171 microbes, so there was so much overlap in who plants associate with because this is presence/absence. it is likly that in some areas microbes would be abundant in roots and in other areas they would be rare, so I'm penalizing the rares here and making them 0
genus4[genus4<0.001]<-0
length(which(genus4==0))

microbecolors.genus<-c("black","chartreuse4","red","purple","lightblue",
                       "darkblue","orange")
plantcolors.genus<-c("darkblue")

plotweb(web=genus4, ybig=1, bor.col.interaction="gray80", col.interaction="gray80",
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)
grouplevel(genus4,index="mean number of links")





###### Read in environmental data #####
env<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/biogeo8.csv",  sep=",", header=TRUE)

#env<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Tulane/Lab/GradStudents/MonicaBrady/env2.csv",  sep=",", header=TRUE)

ind<-as.numeric(unique(genus2$Plot.Number))
env2<-env%>%
  rename(Plot.Number=Sample_name)%>%
  filter(Plot.Number%in%ind)

###I must have done write.csv to write the env2 file as an output b/c it is in Monica's folder and that is the file I sent to her via email


#####  low and hi networks ######

dim(env2)
env3<-arrange(env2,snowdepth)
env3$snowdepth[36:38]
#162 is my cutoff
indlo<-env3$Plot.Number[which(env3$snowdepth<161.5)]
indhi<-env3$Plot.Number[which(env3$snowdepth>161.5)]

genus2lo<-genus2%>%
  filter(Plot.Number%in%indlo)
  
genus2hi<-genus2%>%
  filter(Plot.Number%in%indhi)

genus3lo<-genus2lo%>%
  select(-Plot.Number)%>%
  group_by(Plant.Species)%>%
  summarise_all(mean)

#maybe need to take out microbe zeros? or it looks like bipartite ignores zeros?
genus3lo<-as.data.frame(genus3lo)
rownames(genus3lo)<-as.character(genus3lo[,1])
genus4lo <- as.matrix(genus3lo[,-1])

#delete rare taxa at <0.001 abundance
genus4lo[genus4lo<0.001]<-0
length(which(genus4lo==0))

plotweb(web=genus4lo, ybig=1, bor.col.interaction="gray80", col.interaction=c("gray80"),
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)

genus3hi<-genus2hi%>%
  select(-Plot.Number)%>%
  group_by(Plant.Species)%>%
  summarise_all(mean)
#need to take out microbe zeros?
genus3hi<-as.data.frame(genus3hi)
rownames(genus3hi)<-as.character(genus3hi[,1])
genus4hi <- as.matrix(genus3hi[,-1])

#delete rare taxa at <0.001 abundance
genus4hi[genus4hi<0.001]<-0
length(which(genus4hi==0))

plotweb(web=genus4hi, ybig=1, bor.col.interaction="gray80",  col.interaction="gray80",
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)

networklevel(genus4lo,index="linkage density")
networklevel(genus4hi,index="linkage density")
grouplevel(genus4lo,index="mean number of links")
grouplevel(genus4hi,index="mean number of links")

#calculate mean microbial richness, should equal "mean number of links" above
mean(rowSums(genus4lo>0))
mean(rowSums(genus4hi>0))


#trying to plot as one figure but with lines colored based on whether the interaction is in lo hi or both
#making matrixes with 1's when ther is an association in lo and hi and 0 when there is no association
collo<-rbind(genus4lo[1:6,],CARPER=rep(0,168),genus4lo[7:21,],OXYDIG=rep(0,168),genus4lo[22:23,],STEUMB=rep(0,168),genus4lo[24:25,])
collo[collo>0]<-1
colhi<-rbind(AGRVAR=rep(0,168),ANGGRA=genus4hi[1,],ANTMED=rep(0,168),AQUCOE=rep(0,168),BESALP=rep(0,168),genus4hi[2:5,],CARRUP=rep(0,168),CARSCO=rep(0,168),CERARV=rep(0,168),genus4hi[6:8,],ERISIM=rep(0,168),genus4hi[9:11,],LUZSPI=rep(0,168),MINOBT=genus4hi[12,],OREALP=rep(0,168),genus4hi[13:18,])
colhi[colhi>0]<-2

colall=collo+colhi
colall<-colall[,which(colSums(colall)>0)]  #need to take out microbes that are all zeros or the colors won't work
graphall<-colall
graphall[graphall>0]<-1
colall[colall==0]<-NA
colall[colall==3]<-"gray50" #grey for interactions in both sites
colall[colall==1]<-"royalblue3" #blue for interactions in low snow
colall[colall==2]<-"tomato2" #red for interactions in high snow sites

plotweb(web=graphall, ybig=1, bor.col.interaction=t(colall), col.interaction=t(colall),
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)
#hmm kind of cool!

#tying to figure out colors on a simplified diagram, oxydig only in high, carpyr in both, carrup in lo
test<-rbind(OXYDIG=colall[23,],colall[9:10,])
test[is.na(test)==T]<-0
test2<-test[,which(colSums(test)>0)]
dim(test2)
testna<-test2
testna[testna==0]<-NA
test3<-test2
test3[test3>0]<-1

plotweb(web=test3, ybig=1, bor.col.interaction=t(testna), col.interaction=t(testna),
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)


###I think the next things I would do are:
#1. maybe take out plants that are only present at one site and only do the matrix on plants present in both?
#2. only focus on the extremes of the gradient, rather than just splitting down the middle
#3. and most importantly - circle back to thinking about hypotheses you want to test and what types of analyses you should use (networks or not!) to answer those questions.



######### snowdepth lo, me, hi breaks #########

dim(env2)
env3<-arrange(env2,snowdepth)
env3$snowdepth[23:25]
#148.5 is my lo/me cutoff
env3$snowdepth[47:49]
#188.5 is my me/hi cutoff
indlo<-env3$Plot.Number[which(env3$snowdepth<148.5)]
indme<-env3$Plot.Number[which(env3$snowdepth>148.5&env3$snowdepth<188.5)]
indhi<-env3$Plot.Number[which(env3$snowdepth>188.5)]

genus2lo<-genus2%>%
  filter(Plot.Number%in%indlo)

genus2me<-genus2%>%
  filter(Plot.Number%in%indme)

genus2hi<-genus2%>%
  filter(Plot.Number%in%indhi)


###lo graph
genus3lo<-genus2lo%>%
  select(-Plot.Number)%>%
  group_by(Plant.Species)%>%
  summarise_all(mean)

#maybe need to take out microbe zeros? or it looks like bipartite ignores zeros?
genus3lo<-as.data.frame(genus3lo)
rownames(genus3lo)<-as.character(genus3lo[,1])
genus4lo <- as.matrix(genus3lo[,-1])

#delete rare taxa at <0.001 abundance
genus4lo[genus4lo<0.001]<-0
length(which(genus4lo==0))

plotweb(web=genus4lo, ybig=1, bor.col.interaction="gray80", col.interaction=c("gray80"),
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)


####me graph
genus3me<-genus2me%>%
  select(-Plot.Number)%>%
  group_by(Plant.Species)%>%
  summarise_all(mean)

#maybe need to take out microbe zeros? or it looks like bipartite ignores zeros?
genus3me<-as.data.frame(genus3me)
rownames(genus3me)<-as.character(genus3me[,1])
genus4me <- as.matrix(genus3me[,-1])

#delete rare taxa at <0.001 abundance
genus4me[genus4me<0.001]<-0
length(which(genus4me==0))

plotweb(web=genus4me, ybig=1, bor.col.interaction="gray80", col.interaction=c("gray80"),
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)


####hi graph
genus3hi<-genus2hi%>%
  select(-Plot.Number)%>%
  group_by(Plant.Species)%>%
  summarise_all(mean)
#need to take out microbe zeros?
genus3hi<-as.data.frame(genus3hi)
rownames(genus3hi)<-as.character(genus3hi[,1])
genus4hi <- as.matrix(genus3hi[,-1])

#delete rare taxa at <0.001 abundance
genus4hi[genus4hi<0.001]<-0
length(which(genus4hi==0))

plotweb(web=genus4hi, ybig=1, bor.col.interaction="gray80",  col.interaction="gray80",
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)

networklevel(genus4lo,index="linkage density")
networklevel(genus4me,index="linkage density")
networklevel(genus4hi,index="linkage density")
grouplevel(genus4lo,index="mean number of links")
grouplevel(genus4me,index="mean number of links")
grouplevel(genus4hi,index="mean number of links")

#calculate mean microbial richness, should equal "mean number of links" above
mean(rowSums(genus4lo>0))
mean(rowSums(genus4me>0))
mean(rowSums(genus4hi>0))

#specialization
H2fun(genus4lo,H2_integer = F)
H2fun(genus4me,H2_integer = F)
H2fun(genus4hi,H2_integer = F)
dfun(genus4lo)$d
dfun(genus4me)$d
dfun(genus4hi)$d



######### productivity, lo, me, hi breaks #########

dim(env2)
sum(is.na(env2$TC))#16 NAs
sum(is.na(env2$TN))#16 NAs
sum(is.na(env2$moisture))#10 NAs
sum(is.na(env2$WHC))#11 NAs
env3<-arrange(env2,moisture)
env4<-filter(env3,is.na(moisture)==F)
#the below cutoffs mke 20 plots lo, 21 plots me, and 21 plots hi
env4$moisture[20:22]
#9.1 is my lo/me cutoff
env4$moisture[41:43]
#15.1 is my me/hi cutoff
indlo<-env4$Plot.Number[which(env4$moisture<9.1)]
indme<-env4$Plot.Number[which(env4$moisture>9.1&env4$moisture<15.1)]
indhi<-env4$Plot.Number[which(env4$moisture>15.1)]

genus2lo<-genus2%>%
  filter(Plot.Number%in%indlo)

genus2me<-genus2%>%
  filter(Plot.Number%in%indme)

genus2hi<-genus2%>%
  filter(Plot.Number%in%indhi)

#lo graph
genus3lo<-genus2lo%>%
  select(-Plot.Number)%>%
  group_by(Plant.Species)%>%
  summarise_all(mean)

#maybe need to take out microbe zeros? or it looks like bipartite ignores zeros?
genus3lo<-as.data.frame(genus3lo)
rownames(genus3lo)<-as.character(genus3lo[,1])
genus4lo <- as.matrix(genus3lo[,-1])

#delete rare taxa at <0.001 abundance
genus4lo[genus4lo<0.001]<-0
length(which(genus4lo==0))

plotweb(web=genus4lo, ybig=1, bor.col.interaction="gray80", col.interaction=c("gray80"),
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)

#me graph
genus3me<-genus2me%>%
  select(-Plot.Number)%>%
  group_by(Plant.Species)%>%
  summarise_all(mean)

#maybe need to take out microbe zeros? or it looks like bipartite ignores zeros?
genus3me<-as.data.frame(genus3me)
rownames(genus3me)<-as.character(genus3me[,1])
genus4me <- as.matrix(genus3me[,-1])

#delete rare taxa at <0.001 abundance
genus4me[genus4me<0.001]<-0
length(which(genus4me==0))

plotweb(web=genus4me, ybig=1, bor.col.interaction="gray80", col.interaction=c("gray80"),
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)

#hi graph
genus3hi<-genus2hi%>%
  select(-Plot.Number)%>%
  group_by(Plant.Species)%>%
  summarise_all(mean)
#need to take out microbe zeros?
genus3hi<-as.data.frame(genus3hi)
rownames(genus3hi)<-as.character(genus3hi[,1])
genus4hi <- as.matrix(genus3hi[,-1])

#delete rare taxa at <0.001 abundance
genus4hi[genus4hi<0.001]<-0
length(which(genus4hi==0))

plotweb(web=genus4hi, ybig=1, bor.col.interaction="gray80",  col.interaction="gray80",
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)

networklevel(genus4lo,index="linkage density")
networklevel(genus4me,index="linkage density")
networklevel(genus4hi,index="linkage density")
grouplevel(genus4lo,index="mean number of links")
grouplevel(genus4me,index="mean number of links")
grouplevel(genus4hi,index="mean number of links")

#calculate mean microbial richness, should equal "mean number of links" above
mean(rowSums(genus4lo>0))
mean(rowSums(genus4me>0))
mean(rowSums(genus4hi>0))

#specialization
H2fun(genus4lo,H2_integer = F)
H2fun(genus4me,H2_integer = F)
H2fun(genus4hi,H2_integer = F)
dfun(genus4lo)$d
dfun(genus4me)$d
dfun(genus4hi)$d






#if you want to put this into igraph for visualization do the following
#I'm not sure if this would be nicer or not, didn't get very far
#convert to long format
genus2a<-genus2%>%
  rename(Plant=Group.1)%>%
  filter(Plant!="UNKOPP")%>% #take out unkopp, an unknown plant with only one occurrence
  gather(Microbe,abun,k__Fungi.Other.Other.Other.Other.Other:k__Fungi.p__Zygomycota.c__unidentified.o__unidentified.f__unidentified.g__unidentified)
head(genus2a)

g<-graph_from_data_frame(genus2a, directed = F, vertices = NULL)
bipartite.mapping(g)
V(g)$type <- bipartite_mapping(g)$type
plot(g,vertex.size=4,vertex.color=1,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2






######### microbe comp ordintion #######
#use genus and env from above
genus[1:5,1:5]
dim(genus)

#take out microbial doubletons and singletons
ind<-which((colSums(genus[,3:384]>0)>2))#&colSums(genus[,3:384])>0.006
length(ind)
genusord<-cbind(genus[,1:2],genus[,ind+2])
genusord[1:5,1:5]

#take out plant doubletons and singletons
plantcount<-genusord%>%count(Plant.Species)
plantcount2<-plantcount$Plant.Species[which(plantcount$n>2)]
genusord2<-genusord%>%
  filter(Plant.Species%in%plantcount2)%>%
  arrange(Plot.Number)
genusord2[1:5,1:5]
unique(genusord2$Plant.Species)
genusord2$Plot.Number<-as.numeric(genusord2$Plot.Number)

ind<-as.numeric(unique(genusord2$Plot.Number))
envord<-env%>%
  rename(Plot.Number=Sample_name)%>%
  filter(Plot.Number%in%ind)
allord<-genusord2%>%
  full_join(envord)

allordmic<-allord[,3:215]
allordenv<-allord[,c(1,216:235)]

ind<-which(is.na(rowSums(allordenv[,2:21]))==F)
allordmicnona<-allordmic[ind,]
allordenvnona<-allordenv[ind,]
  
myrda<-capscale(allordmic~Plant.Species,data=allordenv)
myrda<-capscale(allordmic~moisture+Condition(Plant.Species),data=allordenv,na.action=na.omit)
myrda<-capscale(allordmic~pH+snowdepth+moisture+WHC+TN+TC+Plant.Species,data=allordenv,na.action=na.omit)
anova(myrda,by="margin",permutations=how(nperm=9999))#permutations=how(blocks=sample_data(datPlantsS)$Site)

plot(myrda,color=allordenv$Plant.Species)

myrda1<-capscale(allordmicnona~1+Plant.Species,data=allordenvnona,na.action=na.omit)
myrda<-capscale(allordmicnona~pH+snowdepth+moisture+WHC+TN+TC,data=allordenvnona,na.action=na.omit)
myrda<-capscale(allordmicnona~snowdepth+moisture,data=allordenvnona,na.action=na.omit)
myrda2<-ordistep(myrda1,myrda,direction="both",na.action=na.omit)
myrda2<-ordistep(myrda,direction="backward",na.action=na.omit,permutations=how(nperm=9999))



######### plant density vs snow depth and moisture comparisons #####
env10<-env2[which(!is.na(rowSums(env2))),]

#just using the full dataset to make sure i can more or less replicate the ordination in my network paper, yes i can. so variation in the ordination on the reduced dataset is due to the plots that were not included (not a data error)
env11<-env[which(!is.na(rowSums(env))),]
#this is the only way the ordintion will show labels on it. no idea why.
ind2<-as.numeric(seq(0:200))
env12<-env11%>%
  rename(Plot.Number=Sample_name)%>%
  filter(Plot.Number%in%ind2)

plot(env2$snowdepth,env2$pH)
abline(lm(env2$pH~env2$snowdepth))
plot(env2$moisture,env2$TC)
plot(sqrt(env2$snowdepth),sqrt(env2$VascPlantRichness))
rcorr(as.matrix(env2))
m0<-lm(VascPlantRichness~1,data=env10,na.action = na.omit)
m1<-lm(VascPlantRichness~moisture+snowdepth+pH,data=env10,na.action = na.omit)
summary(m1)
step(m0,scope=formula(m1),direction = "forward")

m1<-rda(env12[,c("WHC","moisture","snowdepth","TN","TC","pH","Elev","VascPlantRichness","VascPlantDensity","vascularplantcover","MicC","MicN","NH4","NO3")],scale=T)
m1<-rda(env10[,c("WHC","moisture","snowdepth","TN","TC","pH")],scale=T)
summary(m1)
plot(m1,choices=c(1,2))



######## Plant ordination #######

plant<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Tulane/Lab/GradStudents/MonicaBrady/plantcomp3.csv",  sep=",", header=TRUE)

head(plant)

#take out doubletons singletons, and plots with no plants
#col 89 is lichen, 88 is liverwort, 87 is moss, 86 is unktil (the last vascular plant)
ind<-which((colSums(plant[,2:86]>0)>2))#&colSums(genus[,3:384])>0.006
length(ind)
plantord<-cbind(Sample_name=plant[,1],plant[,ind+1])
plantord[1:5,1:5]

ind<-which((rowSums(plant[,2:86]>0)>0))
plantord2<-plantord[ind,]
  
ind<-as.numeric(unique(plantord2$Sample_name))
envplant<-env%>%
  #rename(Plot.Number=Sample_name)%>%
  filter(Sample_name%in%ind)

allplantord<-plantord2%>%
  full_join(envplant)

allordplant<-allplantord[,2:54]
allordplantenv<-allplantord[,c(1,55:74)]

ind<-which(is.na(rowSums(allordplantenv[,2:21]))==F)
allordplantnona<-allordplant[ind,]
allordplantenvnona<-allordplantenv[ind,]


myrda<-capscale(allordmic~Plant.Species,data=allordenv)
myrda<-capscale(allordmic~moisture+Condition(Plant.Species),data=allordenv,na.action=na.omit)
myrda<-capscale(allordmic~pH+snowdepth+moisture+WHC+TN+TC+Plant.Species,data=allordenv,na.action=na.omit)
anova(myrda2,by="margin",permutations=how(nperm=9999))#permutations=how(blocks=sample_data(datPlantsS)$Site)

plot(myrda2,scaling=3)

myrda1<-capscale(allordplantnona~1,data=allordplantenvnona,na.action=na.omit)
myrda<-capscale(allordplantnona~WHC+pH+snowdepth+moisture+TN+TC,data=allordplantenvnona,na.action=na.omit)#
#myrda<-capscale(allordmicnona~snowdepth+moisture,data=allordenvnona,na.action=na.omit)
myrda3<-ordistep(myrda1,myrda,direction="both",na.action=na.omit,permutations=how(nperm=9999))
myrda2<-ordistep(myrda,direction="backward",na.action=na.omit,permutations=how(nperm=9999))


##### env plus full plant dataset ######
envplantfull<-env%>%
  full_join(plant)

plot(envplantfull$pH,envplantfull$XXMOSS)











########## old ########
#### Order(?) level, not used #####
otu_table_rare_L2T<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Tulane/Lab/GradStudents/MonicaBrady/ITS_TAXA_RAREFIED_5238/otu_table_rare_L2T.csv",  sep=",", header=TRUE)
otu_table_rare_L2T <-separate(data=otu_table_rare_L2T, col=PLANT.PLOT, 
                              into=c("Plant.Species", "Plot.Number"), sep="\\.")
data <- arrange(otu_table_rare_L2T, Plant.Species)
data1 <- subset(data, select=-c(Plot.Number))
data2 <- as.matrix(data1)
data2 <- aggregate.data.frame(x=data1[],by=list(data1$Plant.Species),FUN=mean)
rownames(data2)<-as.character(data2[,1])
data2 <- subset(data2, select=-c(Group.1,Plant.Species))
colnames(data2)<-c("Other","Asco","Basidio","Chytridio","Glomero","Rozello","Zygo")
microbecolors<-c("black","chartreuse4","red","purple",
                 "lightblue","darkblue","orange")
plantcolors<-c("darkblue")
plotweb(web=data2, ybig=1, bor.col.interaction="gray80",  col.interaction=c("gray80"), 
        arrow="down.center", method="normal", labsize=0.9,y.width.high=0.05,
        y.width.low=0.05, col.low=plantcolors, col.high=microbecolors, text.rot=90, 
        low.lab.dis=NULL, bor.col.high="grey47",bor.col.low="gray47", 
        high.lablength=20, low.lablength = 20)

networklevel(data2)
grouplevel(data2,index="mean number of links")
H2fun(web=data2, H2_integer = FALSE)
mean(rowSums(data2>0))

#### Genus level ####
otu_table_rare_L6T<-read.csv("otu_table_rare_L6T.csv",  sep=",", 
                             header=TRUE,comment.char = "#")
genus <-separate(data=otu_table_rare_L6T, col=PLOT.ID, into=c("Plant.Species", 
                                                              "Plot.Number"), sep="\\.")
genus <- arrange(genus, Plant.Species)
genus1 <- subset(genus, select=-c(Plot.Number))
genus2 <- aggregate.data.frame(x=genus1[,-1],by=list(genus1$Plant.Species),FUN=mean)
genus3 <- subset(genus2, select=-c(Group.1))
rownames(genus3)<-c(as.character(genus2[,1]))
genus4 <- as.matrix(genus3)

microbecolors.genus<-c("black","chartreuse4","red","purple","lightblue",
                       "darkblue","orange")
plantcolors.genus<-c("darkblue")
plotweb(web=genus4, ybig=1, bor.col.interaction="gray80",  col.interaction="gray80",
        arrow="down.center", method="normal", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)
plotweb(web=genus4, ybig=1, bor.col.interaction="gray80",  col.interaction="gray80",
        arrow="down.center", method="cca", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low=plantcolors.genus, 
        col.high=microbecolors.genus, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)
networklevel(genus4)
H2fun(web=genus4, H2_integer = FALSE)


asco<-subset(genus3, select = OTU.1:OTU.238) 
  #subsets dataframe to only genera in ascomycota
asco1 = subset(asco, select = -c(OTU.12,OTU.14,OTU.19,OTU.20,OTU.35,OTU.39,OTU.48,
      OTU.53,OTU.57,OTU.58,OTU.60,OTU.61,OTU.62,OTU.72,OTU.73,OTU.78,OTU.103,OTU.113,
      OTU.115,OTU.117,OTU.122,OTU.127,OTU.139,OTU.147,OTU.157,OTU.159,OTU.164,
      OTU.173,OTU.174,OTU.180,OTU.181,OTU.182,OTU.188,OTU.190,OTU.198,OTU.203,
      OTU.207,OTU.214,OTU.219,OTU.222,OTU.224,OTU.229,OTU.231,OTU.235,OTU.237,
      OTU.238))#removes unidentified genuses in ascomycota
asco2 = subset(asco1, select = -c(OTU.1, OTU.2, OTU.10, OTU.13, OTU.16, OTU.26, 
      OTU.29, OTU.42, OTU.50, OTU.54, OTU.63, OTU.64, OTU.66, OTU.79, OTU.80, 
      OTU.101, OTU.105, OTU.123, OTU.128, OTU.131, OTU.140, OTU.148, OTU.162, 
      OTU.176, OTU.177, OTU.192, OTU.199, OTU.200, OTU.201, OTU.204, OTU.208, 
      OTU.218, OTU.223, OTU.230))
asco3<-as.matrix(asco2)
plotweb(web=asco3, ybig=1, bor.col.interaction="gray80",  col.interaction="gray80", 
        arrow="down.center", method="normal", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, col.low="darkblue", 
        col.high="chartreuse4",
        text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47")

asco.c1= subset(asco3, select = c(OTU.3:OTU.100,OTU.183:OTU.197))
microbecolors.c1<-c(rep("black" , 1), rep("chartreuse4",37), rep("red", 10), 
                 rep("purple", 20), rep("lightblue",6), rep("orange",6))
plantcolors.c1<-c("darkblue")

plotweb(web=asco.c1, ybig=1, bor.col.interaction="gray80",  col.interaction="gray80", 
        arrow="down.center", method="normal", labsize=0.9,y.width.high=0.05,
        y.width.low=0.05, col.low=plantcolors.c1, col.high=microbecolors.c1, text.rot=90, 
        low.lab.dis=NULL, bor.col.high="grey47",bor.col.low="gray47", 
        high.lablength=20, low.lablength = 20)


otu.abd<-colSums(genus4)
genus5<-t(genus3)
genus6<-as.data.frame(genus5)
genus7<-tibble::rownames_to_column(genus6, var="rowname")
genus8<-as.data.frame(arrange(genus7, desc(otu.abd)))
genus9<-tibble::column_to_rownames(genus8, var="rowname")
otu.abd1<-rowSums(genus9)
genus10<-t(genus9)
genus11<-genus10[,c(1:36)]
microbecolors.genus1<-c("black","chartreuse4","red","purple","lightblue",
                       "darkblue","orange")
plantcolors.genus1<-c("darkblue")
plotweb(web=genus11, ybig=1, bor.col.interaction="gray80",  col.interaction="gray80",
        arrow="down.center", method="normal", labsize=0.9,
        y.width.high=0.05,y.width.low=0.05, arrow="both.center", col.low=plantcolors.genus1, 
        col.high=microbecolors.genus1, text.rot=90, low.lab.dis=NULL,
        bor.col.high="grey47",bor.col.low="gray47", high.lab.dis = 0)
networklevel(genus11)
H2fun(web=genus11, H2_integer = FALSE)


