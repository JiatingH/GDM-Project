#load packages
library(ggplot2)
library(pheatmap)

# load data for analysis
args = commandArgs(T)
d = read.table(args[1], sep = ',',header = 1, row.names = 1)
data =t(d)
#samples as columns, species as rows
#DA function
DA = function(data,p_val = c()){
    for(n in rownames(data)){
        gdm = as.numeric(data[n,1:42]) # group1
        ngt = as.numeric(data[n,43:84]) # group2
        #print(data[n,47:88])
        wilcox = wilcox.test(gdm,ngt,exact = FALSE,paired = FALSE, correct = FALSE)
    
        p_val= c(p_val,wilcox$p.value)
        frame = data.frame(p = p_val)
    }
    #print(frame)
    frame$rownames = rownames(data)
    order
    order = frame[order(frame$p),]
    #p.adjust
    order$fdr = p.adjust(order$p, method = 'BH')
    order[which(order$fdr<0.2),]
}

da.data = DA(data)

DA.Data = data[da.data$rownames,]

pheatmap(log(DA.Data+0.05),
	cluster_row = FALSE,
	cluster_col = FALSE,
	cellheight = 12,
	cellwidth = 5,
	fontsize_col = 5,
	file = 'DA.Data.Heatmap.pdf' #save the plot as pdf
	)


