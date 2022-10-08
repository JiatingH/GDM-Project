# load packages
reticulate::py_module_available("shap")
library("shapper")
library("DALEX")
library("randomForest")

library(Hmisc)
library(reshape2)
library(ggsignif)

#data for analysis
args = commandArgs(T)
species = read.table(args[1],sep = ',',header = 1, row.names = 1)
papi = read.table(args[2],sep = ',',header = 1, row.names = 1)

# share sample data
co.sample = intersect(rownames(species),rownames(papi))
#print(co.sample)

sample.info = data.frame(Group = c(rep('GDM',25),rep('NGT',24)))
rownames(sample.info) = co.sample

# papi pathway with significant abundance
DA_PAPI = function(data,p_val = c()){
    for(n in rownames(data)){
        gdm = as.numeric(data[1:33,n])
        ngt = as.numeric(data[34:70,n])
        #print(data[n,47:88])
        wilcox = wilcox.test(gdm,ngt,exact = FALSE,paired = FALSE, correct = FALSE)
    
        p_val= c(p_val,wilcox$p.value)
        frame = data.frame(p = p_val)
    }
    #print(frame)
    frame$rownames = colnames(data)
    order = frame[order(frame$p),]
    #p.adjust
    order$fdr = p.adjust(order$p, method = 'BH')
    order[which(order$fdr<0.2),]
    
}

da.papi = DA_PAPI(papi)
da.papi.data = papi[,da.papi$rownames]
co.papi = da.papi[co.sample,]
sample.metadata = cbind(sample.info, co.papi)

species.profiles.input = species[co.sample,]
species.profiles.input$path1 = co.papi[,da.papi$rownames[1]]
species.profiles.input$path2 = co.papi[,da.papi$rownames[2]]
species.profiles.input$path3 = co.papi[,da.papi$rownames[3]]

# path 1
Y_train <- species.profiles.input$path1
x_train <- species.profiles.input[,1:50]

set.seed(123)
model_rf <- randomForest(x = x_train, y = Y_train)
ive_rf.predict.all=c()
ive_rf.importance.all=matrix(,ncol(x_train),1)

for (i in 1:nrow(x_train))
{

ive_rf <- individual_variable_effect(model_rf, data = x_train,
                                     new_observation = x_train[i,])

ive_rf.results=cbind(as.matrix(ive_rf$"_yhat_"),as.matrix(ive_rf$"_vname_"),as.matrix(ive_rf$"_attribution_"),as.matrix(ive_rf$"_sign_"))
colnames(ive_rf.results)=c("predict.score","variable.name","attribution","sign")

#write.table(ive_rf.results,file=paste("species.SBP.shap_results.",rownames(x_train)[i],".txt",sep=""),sep="\t",col.names=NA)


ive_rf.predict.all=c(ive_rf.predict.all,ive_rf$"_yhat_"[1])

ive_rf.importance.all=cbind(ive_rf.importance.all,as.matrix(ive_rf$"_attribution_"))
rownames(ive_rf.importance.all)=ive_rf$"_vname_"


}


ive_rf.predict.all=as.matrix(ive_rf.predict.all)
rownames(ive_rf.predict.all)=rownames(x_train)

ive_rf.importance.all=ive_rf.importance.all[,-1]

colnames(ive_rf.importance.all)=rownames(x_train)

write.table(ive_rf.importance.all,file="species_papi_path1.shap_model.importance.txt",sep="\t",col.names=NA)
#write.table(ive_rf.predict.all,file="species_papi_path1.shap_model.predict_score.txt",sep="\t",col.names=NA)

p1 = ive_rf.importance.all

DA_SHAP = function(data,p_val = c()){
    for(n in rownames(data)){
        gdm = as.numeric(data[n,1:21])
        ngt = as.numeric(data[n,22:46])
        #print(data[n,47:88])
        wilcox = t.test(gdm,ngt,paired = FALSE,exact = FALSE)
    
        p_val= c(p_val,wilcox$p.value)
        frame = data.frame(p = p_val)
    }
    #print(frame)
    frame$rownames = rownames(data)
    order = frame[order(frame$p),]
    #p.adjust
    order$fdr = p.adjust(order$p, method = 'BH')
    f = order[which(order$fdr<0.05),]
    return(f)
}

da.p1 = DA_SHAP(p1)
da.path1 = p1[da.p1$rownames,]
da.path1$species = rownames(da.path1)

melt.p1 = melt(da.path1, id.vars = 'species')
melt.p1$Group = c(rep('GDM',21*length(rownames(da.path1))), rep('NGT', 25*length(rownames(da.path1))))

pdf(file = 'SHAP.value.for.path1.pdf')
ggplot(melt.p1, aes(x=species, y=value,color=Group))+ 
  geom_boxplot(width=0.8,aes(fill = Group),alpha = 0.3)+theme_classic(base_size = 16)+scale_color_brewer(palette="Set1")+
  theme(axis.text.x = element_text(size=16,angle=90,color = 'black'),
        axis.text.y = element_text(size = 16,color = 'black'))+
     scale_y_continuous(name = 'Importance in SHAP model',limits = c(-25, 28))+
geom_signif(y_position=c(18,12,25), xmin=c(0.8,1.8,2.8), xmax=c(1.2,2.2,3.2),  
              annotation=c("**","**",'*'), tip_length=0.03, size=0.6, textsize = 6,
              vjust = 0.3,color = 'black')   
dev.off()

#repeat the analysis for more differential abundant pathway
