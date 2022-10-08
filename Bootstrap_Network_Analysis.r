#load packages
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(igraph)

#data for analysis
args = commandArgs(T)

species = read.table(args[1], sep = ',',header = 1, row.names = 1)

#gdm 42 ngt 42
gdm = species[1:42,]
ngt = species[43：84，]

#bootstrap for 100 times
extract = function(data.gdm,data.ngt,save){
    random_1 = data.gdm[sample(nrow(data.gdm), 30,replace = F), ]
    random_2 = data.ngt[sample(nrow(data.ngt), 30,replace = F), ]
    
    #save tables in the path you like
    write.table(random_1,file=paste("/home/jiating/gdm_t",save,".txt",sep=""),sep = '\t',col.names =NA)
    write.table(random_2,file=paste("/home/jiating/ngt_t",save,".txt",sep=""),sep = '\t',col.names =NA)
}

for(n in 1:100){
    extract(gdm,ngt,n)
}

#the bootstrap table link
table.gdm.path = '/home/jiating/GDM_random'
table.ngt.path = '/home/jiating/NGT_random'

table.gdm = list.files(table.gdm.path, pattern = '*.txt$',full.name = TRUE)
table.ngt = list.files(table.ngt.path, pattern = '*.txt$',full.name = TRUE)

#spearman correlation
correlation = function(path){
for(i in 1:100){
   
    
    #read data
   data = read.table(path[i], sep = '\t',header =1, row.names = 1)
    data = as.matrix(data)
    # spearman correlation
    r.p = rcorr(data,type = 'spearman')
    pval = r.p$P
    rval = r.p$r
    
    #filter r > 0.3, p>0.05
    rval[which(pval> 0.05)] = 0
    rval[which(abs(rval) < 0.3 )] =0
    #rval = rval[which(rowSums(rval)!=1),]
    rval[which(rval == 1)] = 0
    #rval = rval[which(rSums(rval)!=1),]
    
    melt_data = subset(melt(rval), value != 0)
    melt_data$R = abs(melt_data$value)
    
    melt_data$R[melt_data$R == 1] = NA
    
    d = na.omit(melt_data)
    
    #write.table(d,file=paste("/home/jiating/GDM_correlation/gdm_cor_",c(seq(1,100,1))[i],".txt",sep=""),sep = '\t',col.names =NA)

    #write.table(d,file=paste("/home/jiating/NGT_correlation/ngt_cor_",c(seq(1,100,1))[i],".txt",sep=""),sep = '\t',col.names =NA)
}
}

#bootstrap correlation table link
cor1.path = '/home/jiating/GDM_correlation'
cor2.path = '/home/jiating/NGT_correlation'

cor.gdm = list.files(cor1.path,pattern = '*.txt$',full.name = TRUE)
cor.ngt = list.files(cor2.path,pattern = '*.txt$',full.name = TRUE)


# network construction analysis

#transitivity
T = function(cor.path){
transitivity = data.frame()
for(i in 1:100){
   
    
    #read data
    cor = read.table(cor.path[i], sep = '\t',header =1, row.names = 1)
    #degree
   cor.G = graph.data.frame(cor, directed = FALSE)
   G.degrees = transitivity(cor.G)

   G.degrees <- data.frame(value = G.degrees)
   #G.degrees$Species = rownames(G.degrees)
   #G.degree = G.degrees[order(G.degrees$eccentricity,decreasing = TRUE),]
    
    transitivity = rbind(transitivity ,data.frame(G.degrees))
  
}
    return(transitivity)
    }



#eigen centrality
EC = function(cor.path){

eigen_centrality = data.frame()
#eign_vector = data.frame()
    
for(i in 1:100){
   
    
    #read data
    cor = read.table(cor.path[i], sep = '\t',header =1, row.names = 1)
    #degree
   cor.G = graph.data.frame(cor, directed = FALSE)
   G.degrees = eigen_centrality(cor.G)

   G.value <- data.frame(value = G.degrees$value)
   # G.vector = data.frame(eign_vector =G.degrees$vector )
   # G.vector$Species= rownames(G.vector)
    #G.value$Species = rownames(G.value)
   #G.degrees$Species = rownames(G.degrees)
   #G.degree = G.degrees[order(G.degrees$eccentricity,decreasing = TRUE),]
    
    eigen_centrality = rbind(eigen_centrality ,data.frame(G.value))
    #eign_vector = rbind(eign_vector,data.frame(G.vector))
    
  
}

return(eigen_centrality)
}


#eign vector
EV = function(cor.path){

eigen_centrality = data.frame()
eign_vector = data.frame()
    
for(i in 1:100){
   
    
    #read data
    cor = read.table(cor.path[i], sep = '\t',header =1, row.names = 1)
    #degree
   cor.G = graph.data.frame(cor, directed = FALSE)
   G.degrees = eigen_centrality(cor.G)

   #G.value <- data.frame(value = G.degrees$value)
   G.vector = data.frame(value =G.degrees$vector )
   G.vector$Species= rownames(G.vector)
    #G.value$Species = rownames(G.value)
   #G.degrees$Species = rownames(G.degrees)
   #G.degree = G.degrees[order(G.degrees$eccentricity,decreasing = TRUE),]
    
    #eigen_centrality = rbind(eigen_centrality ,data.frame(G.value))
    eign_vector = rbind(eign_vector,data.frame(G.vector))
    
  
}

return( eign_vector)
}



#degree
D = function(cor.path){
G = data.frame()
for(i in 1:100){
   
    
    #read data
    cor = read.table(cor.path[i], sep = '\t',header =1, row.names = 1)
    #degree
   cor.G = graph.data.frame(cor)
   G.degrees = degree(cor.G)

   G.degrees <- data.frame(value = G.degrees)
   #G.degrees$Species = rownames(G.degrees)
   #G.degree = G.degrees[order(G.degrees$D,decreasing = TRUE),]
    
    G = rbind(G,data.frame(G.degree))
  
   
}
    return(G)
    }




# weighted degree
WD = function(cor.path){
WG = data.frame()
for(i in 1:100){
   
    
    #read data
    cor = read.table(cor.path[i], sep = '\t',header =1, row.names = 1)
    #degree
   cor.G = graph.data.frame(cor, directed = FALSE)
   G.degrees = strength(cor.G, vids = V(cor.G),loops = TRUE,weights = E(cor.G)$R)

   G.degrees <- data.frame(value = G.degrees)
   G.degrees$Species = rownames(G.degrees)
   G.degree = G.degrees[order(G.degrees$value,decreasing = TRUE),]
    
    WG = rbind(WG,data.frame(G.degree))
  
}
    return(WG)

}



#closeness
C = function(cor.path){

closeness = data.frame()
for(i in 1:100){
   
    
    #read data
    cor = read.table(cor.path[i], sep = '\t',header =1, row.names = 1)
    #degree
   cor.G = graph.data.frame(cor, directed = FALSE)
   G.degrees = closeness(cor.G)

   G.degrees <- data.frame(value = G.degrees)
   G.degrees$Species = rownames(G.degrees)
   G.degree = G.degrees[order(G.degrees$value,decreasing = TRUE),]
    
    closeness = rbind(closeness,data.frame(G.degree))
  
}
    return(closeness)
    }



#weighted closeness
WC = function(cor.path){

closeness = data.frame()
for(i in 1:100){
   
    
    #read data
    cor = read.table(cor.path[i], sep = '\t',header =1, row.names = 1)
    #degree
   cor.G = graph.data.frame(cor, directed = FALSE)
   G.degrees = closeness(cor.G,weights = E(cor.G)$R)

   G.degrees <- data.frame(value = G.degrees)
   G.degrees$Species = rownames(G.degrees)
   G.degree = G.degrees[order(G.degrees$value,decreasing = TRUE),]
    
    closeness = rbind(closeness,data.frame(G.degree))
  
}
    return(closeness)
    }




#betweenness
B = function(cor.path){
Betweenness = data.frame()
for(i in 1:100){
   
    
    #read data
    cor = read.table(cor.path[i], sep = '\t',header =1, row.names = 1)
    #degree
   cor.G = graph.data.frame(cor, directed = FALSE)
   G.degrees = betweenness(cor.G)

   G.degrees <- data.frame(value = G.degrees)
   G.degrees$Species = rownames(G.degrees)
   G.degree = G.degrees[order(G.degrees$value,decreasing = TRUE),]
    
    Betweenness = rbind(Betweenness,data.frame(G.degree))
  
}
    return(Betweenness)
    }




gdm.t = T(cor.gdm)
ngt.t = T(cor.ngt)

gdm.ec = EC(cor.gdm)
ngt.ec = EC(cor.ngt)

gdm.ev = EV(cor.gdm)
ngt.ev = EV(cor.ngt)

gdm.b = B(cor.gdm)
ngt.b = B(cor.ngt)

gdm.wd = WD(cor.gdm)
ngt.wd = WD(cor.ngt)

gdm.wc = WC(cor.gdm)
ngt.wc = WC(cor.ngt)


#T EC #EV #B #WD #WC
#compare using wilcoxon test
wilcoxon = function(data1, data2){
    r= wilcox.test(data1[,1],data2[,1],paired = FALSE)
    return(r)
}

wilcoxon(gdm.t,ngt.t)
wilcoxon(gdm.ec,ngt.ec)
wilcoxon(gdm.ev,ngt.ev)
wilcoxon(gdm.b,ngt.b)
wilcoxon(gdm.wd,ngt.wd)
wilcoxon(gdm.wc,ngt.wc)

#T EC #EV #B #WD #WC
# combine the gdm and ngt data
bind = function(gdm,ngt){
    gdm$Group = c('GDM')
    ngt$Group = c('NGT')
    
    data = rbind(gdm,ngt)
    
    return(data)
}

t = bind(gdm.t, ngt.t)
ec = bind(gdm.ec, ngt.ec)
ev = bind(gdm.ev, ngt.ev)
b = bind(gdm.b, ngt.b)
wd = bind(gdm.wd,ngt.wd)
wc = bind(gdm.wc, ngt.wc)

#mean rank in these networks
otu = function(data){
species = data$Species[!duplicated(data$Species)]
    return(species)
}
    
DataMean = function(data,species){
    mean = c()
for(n in species){
    dd = data[which(data$Species ==n),]
    d = data.frame(dd[1])
    mean = c(mean,colMeans(d))
}
    dataframe = data.frame(Mean = mean)
    dataframe$Species = species
    data_order = dataframe[order(dataframe$Mean,decreasing = TRUE),]
    data_order$species = factor(data_order$Species, levels = data_order$Species)
    
    return(data_order)

}

#EV #B #WD #WC
gdm.ev.s = otu(gdm.ev)
ngt.ev.s = otu(ngt.ev)

gdm.b.s = otu(gdm.b)
ngt.b.s = otu(ngt.b)

gdm.wd.s = otu(gdm.wd)
ngt.wd.s = otu(ngt.wd)

gdm.wc.s = otu(gdm.wc)
ngt.wc.s = otu(ngt.wc)

gdm.ev.mean = DataMean(gdm.ev,gdm.ev.s)
ngt.ev.mean = DataMean(ngt.ev,ngt.ev.s)

gdm.b.mean = DataMean(gdm.b,gdm.b.s)
ngt.b.mean = DataMean(ngt.b,ngt.b.s)

gdm.wd.mean = DataMean(gdm.wd,gdm.wd.s)
ngt.wd.mean = DataMean(ngt.wd,ngt.wd.s)

gdm.wc.mean = DataMean(gdm.wc,gdm.wc.s)
ngt.wc.mean = DataMean(ngt.wc,ngt.wc.s)

gdmrank = function(data){
ggplot(data=data, mapping=aes(x=species ,y=Mean))+
  geom_bar(stat="identity",width = 0.8, position = 'dodge',fill = '#F8766D')+theme(axis.text.x = element_text(size=10,angle=90,color = 'black'))+
xlab('Species')+ylab('Mean of 100 Bootstraps in GDM Network (25 samples)')
}

ngtrank = function(data){
ggplot(data, mapping=aes(x=species ,y=Mean))+
  geom_bar(stat="identity",width = 0.8, position = 'dodge',fill = '#619CFF')+theme(axis.text.x = element_text(size=10,angle=90,color = 'black'))+
xlab('Species')+ylab('Degree Mean of 100 Bootstraps in NGT Network (30 samples)')
    }

#example plot
pdf(file = 'ranking_gdm_ev.pdf' )
gdmrank(gdm.ev.mean)
dev.off()

pdf(file = 'ranking_ngt_wc.pdf' )
ngtrank(ngt.wc.mean)
dev.off()


