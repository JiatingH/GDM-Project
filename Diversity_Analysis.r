# load packages
library(vegan)
library(ggplot2)
library(ggsignif)

# load the relative abundance table for analysis
args = commandArgs(T)

species = read.table(args[1],sep = ',',header = 1, row.names = 1)

# please transform the data to keep the table rowname as sample ID, colname as species 

#Bray-Curtis distance
bray = vegdist(species, method = 'bray')
mds.bray = cmdscale(bray)
MDS = data.frame(MDS1 = mds.bray[,1], MDS2 = mds.bray[,2])

#group tab
group = c(rep('GDM', 42),rep('NGT',42))
MDS$group = group

# The ordinated plot of Bray-Curtis distance of the species relative abundance 

pdf(file = 'MDS.pdf') #save plot as pdf file
ggplot(MDS, aes(x = MDS1, y = MDS2, group = group, color = group))+
geom_point()+scale_color_brewer(palette = 'Set1')+
stat_ellipse(alpha = 1, linetype = 5, type = 't')+
theme_classic(base_size = 16)
dev.off()

# ANOSIM test for the difference of beta diversity
anosim.result = anosim(bray, group, permutations = 999)

# Shannon index & Richness

alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
    if (method == 'richness') result <- rowSums(x > 0)    #richness index
    else if (method == 'chao1') result <- estimateR(x)[2, ]    #Chao1 index
    else if (method == 'ace') result <- estimateR(x)[4, ]    #ACE index
    else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)    #Shannon index
    else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson index
    else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou evenness
    else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
    
    result
}

shannon.index = alpha_index(species, method = 'shannon')
rich = alpha_index(species, method = 'richness', base = exp(1))

#Wilcoxon test
shannon.test = wilcox.test(shannon.index[1:42],shannon.index[43:84],paired = FALSE)
#print(shannon.test)
rich.test = wilcox.test(rich[1:42], rich[43:84], paired = FALSE)
#print(rich.test)

Shannon.Index = data.frame(shannon.index = shannon.index, group = c(rep('GDM',42),rep('NGT',42)))
Richness = data.frame(richness = rich, group = c(rep('GDM',42),rep('NGT',42)))

pdf(file = 'shannon.index.pdf')
ggplot(Shannon.Index, aes(x=group, y=shannon.index,color=group))+ 
  geom_boxplot(width=0.8,aes(fill = group),alpha = 0.3)+theme_classic(base_size = 16)+scale_color_brewer(palette="Set1")+
  theme(axis.text.x = element_text(size=16,angle=0,color = 'black'),
        axis.text.y = element_text(size = 16,color = 'black'))+geom_jitter(aes(fill=group),width =0.2)
dev.off()


pdf(file = 'richness.pdf') # save the plot as pdf file
ggplot(Richness, aes(x=group, y=richness,color=group))+ 
  geom_boxplot(width=0.8,aes(fill = group),alpha = 0.3)+theme_classic(base_size = 16)+scale_color_brewer(palette="Set1")+
  theme(axis.text.x = element_text(size=16,angle=0,color = 'black'),
        axis.text.y = element_text(size = 16,color = 'black'))+geom_jitter(aes(fill=group),width =0.2)+
geom_signif(y_position=c(120), xmin=c(1), xmax=c(2),  
              annotation=c("*"), tip_length=0.03, size=0.6, textsize = 12,  
              vjust = 0.3,color = 'black')  
dev.off()




