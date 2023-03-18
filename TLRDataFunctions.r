#load libraries 

library(phytools)
library(ape)
library(tidyverse)
library(geiger)
library(ggstream)

#data paths
datapath<-"~/Documents/TLR_6_Project/ComparativeData/TLR_gene_history_R/TLR_copy_number_data.csv"
treepath<-"~/Documents/TLR_6_Project/ComparativeData/TLR_gene_history_R/timetree.tre"

#directory for saving
setwd("~/Documents/TLR_6_Project/ComparativeData/Routput")

#read in data
tree<-read.tree(treepath)
data<-read.csv(datapath)

#data cleanup
#first replace space with underscore
data<-data %>% mutate(Genus.species= str_replace(Genus.species, " ", "_"))

#check for mismatches
sample<-data[,2]
names(sample)<-data[,1]
treedata(tree, sample)
#replace names
tree$tip.label[tree$tip.label=="Poecilia_reticulata"]<-"Poecilia_formosa"
tree$tip.label[tree$tip.label=="Carcharias_taurus"]<-"Carcharodon_carcharias"
tree$tip.label[tree$tip.label=="Etheostoma_caeruleum"]<-"Etheostoma_spectabile"
tree$tip.label[tree$tip.label=="Periophthalmus_argentilineatus"]<-"Periophthalmus_magnuspinnatus"
tree$tip.label[tree$tip.label=="Polypterus_endlicherii"]<-"Polypterus_endlicheri"
tree$tip.label[tree$tip.label=="Sinibrama_macrops"]<-"Anabarilius_grahami"
tree$tip.label[tree$tip.label=="Ara_ararauna"]<-"Strigops_habroptila"
tree$tip.label[tree$tip.label=="Rhea_americana"]<-"Struthio_australis"

#check again
treedata(tree, sample)

#ok, lets get all the origins
TLR_TT<-NULL
for(i in 2:24){
	geneTarget<-data[data[,i]>0,]
	locus<-geneTarget[,i]
	names(locus)<-geneTarget[,1]
	td<-treedata(tree, locus)
	index<-i-1
	TLR_TT[index]<-max(branching.times(td$phy))
}
gtt<-(table(TLR_TT))
gtt <-rev(gtt)
gtt<-cumsum(gtt)

#plot it
gtt<-as.data.frame(gtt) 
gtt$times=as.numeric(row.names(gtt))
gtt %>% ggplot(aes(times, log(gtt))) +
geom_line() + 
scale_x_reverse() +
theme_classic()

#now to define losses
#first fit the trait models
models<-matrix(nrow=23, ncol=3)
for(i in 2:24){
locus<-data[,i]
names(locus)<-data[,1]
td<-treedata(tree, locus)
index<-i-1
models[index,]<-c(AIC(ace(td$data,td$phy,model="ER",type="discrete")),
AIC(ace(td$data,td$phy,model="SYM",type="discrete")),
AIC(ace(td$data,td$phy,model="ARD",type="discrete")))
}
models2<-as_tibble(models)
names(models2)<-c("ER","SYM", "ARD")
models2<-models2 %>% rowwise() %>% mutate(bestfit=min(ER, SYM, ARD), dAIC=-1*sum(ER, SYM, ARD)+max(ER, SYM, ARD)+ min(ER, SYM, ARD) + min(ER, SYM, ARD))
modelsTested<-c("ER", "SYM", "ARD") 
ModelName<-NULL
for (i in 1:23){
  target<-models2[i,c(1:3)]
  integer<-which(target==models2$bestfit[i])
  ModelName[i]<-modelsTested[integer] 
}

#write out the models for supplemental
#write.csv(file="TLR_Model_fits", models2)

#use best-fit models with simmap 
for(i in 2:24){
# Simulate 200 stochastic character maps
locus<-data[,i]
names(locus)<-data[,1]
td<-treedata(tree, locus)
index<-i-1
chartrees <- make.simmap(td$phy, locus,
                         model=ModelName, nsim = 1000)
# Output summary information
obj<-summary(chartrees)
fn<-paste("ParalogTLR",names(data)[i])
pdf(file=fn, width=8, height=10)
plot(obj)
cols<-setNames(palette()[1:length(unique(getStates(chartrees,
    "tips")))],sort(unique(getStates(chartrees,"tips"))))
add.simmap.legend(prompt=FALSE, colors=cols, x=5, y=10, leg=0:length(unique(getStates(chartrees))))

## stop pdf saving
dev.off() }

#make a negate function
#`%!in%` <- Negate(`%in%`)
#gene loss detection, note this will make warnings with treedata since it reports what is dropped
placeholder<--1
nameholder<-"a"
for(i in 1:23){
#isolate each TLR tree from its origin
index<-i+1
locus<-data[,index]
names(locus)<-data[,1]
tdtemp<-treedata(tree, locus)
geneMRCA<-locus[locus>0]
origin<-findMRCA(tdtemp$phy,names(geneMRCA))
newtree<-extract.clade(tdtemp$phy,origin)
td<-treedata(newtree, locus)
#check if the there are losses and place these together
if(0 %in% td$data){
  data2<-td$data[td$data==0,]
#going through each taxon to see where losses occur
for (j in 1:length(names(data2))){
  print(j)
taxon<-which(td$phy$tip.label==names(data2[j]))
firstnode<-getParent(td$phy,taxon)
nodes<-getDescendants(td$phy, firstnode)
ingroup<-nodes[nodes<length(td$phy$tip.label)+1]		
#if this is true, keep the divergence time as a loss
	if (sum(td$data[rownames(td$data)%in%td$phy$tip.label[ingroup],])>0)
	{
			placeholderi<-max(branching.times(drop.tip(td$phy,td$phy$tip.label[-ingroup])))	
			} 
#if this is not true, find the loss
else {
	  test<-sum(td$data[rownames(td$data)%in%td$phy$tip.label[ingroup],])
	  tracker<-findMRCA(td$phy, ingroup)
	  while(test < 0.5) {
	    oldest<-findMRCA(td$phy, ingroup)
	    nextnode<-getParent(td$phy,oldest)
	    nodes2<-getDescendants(td$phy, nextnode)
	    ingroup<-nodes2[nodes2<length(td$phy$tip.label)+1]
	    test<-sum(td$data[rownames(td$data)%in%td$phy$tip.label[ingroup],])
	    print("while")
	  }
	  nodeX<-nextnode
	  finalnode<-getDescendants(td$phy, nodeX[length(nodeX)])
	  finalnode<-finalnode[finalnode<length(td$phy$tip.label)]
	  placeholderi<-max(branching.times(drop.tip(td$phy,td$phy$tip.label[-finalnode])))
		}
placeholder<-c(placeholder,placeholderi)
geneNames<-rep(names(data[index]),length(placeholderi))
nameholder<-c(nameholder,geneNames)
}
}
}

#Assemble losses through time for plotting
assembledLoss<-data.frame(time=placeholder[-1],gene=nameholder[-1]) %>% 
  group_by(gene) %>%
    summarize(lossTime=unique(time))
   
lossTibble<- assembledLoss %>%
   group_by(gene, lossTime) %>%
   summarize(n())

#order for plotting
lossTibble$gene<- factor(lossTibble$gene, 
levels = c("TLR1","TLR3","TLR4","TLR5","TLR7","TLR8","TLR9","TLR10","TLR11","TLR12","TLR13","TLR18","TLR19","TLR20","TLR21","TLR22","TLR25","TLR27","Novel.3..a."),
ordered = TRUE)
  
 
lossTibble %>% 
  ggplot(aes(x=lossTime, y=`n()`, fill=gene)) +
  geom_stream(type="ridge") + 
  ylab("Number of losses") +
  xlab("Time from present") +
  scale_fill_viridis(discrete=TRUE) +
  theme_classic() +
  guides(fill = guide_legend(title = "Gene"))

