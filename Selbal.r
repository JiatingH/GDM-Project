library(selbal)
library(devtools)
library(Hmisc)

#setwd("set your work data path")

args=commandArgs(T)

data = read.table(args[1], sep = ',',header =1,row.names =1)
tab = read.table(args[2],sep = ',',header = 1, row.names =1)


type = tab[,1][!duplicated(tab[,1])]

if(length(type)<5){
tab_order = data.frame(data = tab[colnames(data),])
rownames(tab_order) = colnames(data)  

y = factor(tab_order[,1])}else{

tab_order = data.frame(data = tab[colnames(data),])
rownames(tab_order) = colnames(data)    
tab_fillna = impute(tab_order[,1],mean)

y = as.vector(tab_fillna)
}

        

#d = t(data)

group.dic = selbal.cv(x = data, y = y, n.fold = 5, n.iter = 10, logit.acc = 'AUC', user_numVar = 10)

pdf(file = 'selbal_accuracy.pdf' )
group.dic$accuracy.nvar
dev.off()

pdf(file = 'selbal_cross_validation.pdf' )
group.dic$var.barplot
dev.off()

pdf(file = 'selbal_global.pdf' )
grid.draw(group.dic$global.plot)
dev.off()

pdf(file = 'selbal_cross_table.pdf' )
plot.tab(group.dic$cv.tab)
dev.off()
