# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("impute")
install.packages('spls')

# load packages
library(pheatmap)
library(impute)
library(spls)
library(mice)
library(vegan)

# data for analysis
args = commandArgs(T)

species = read.table(args[1], sep = ',',header = 1, row.names = 1)
metabolites = read.table(args[2],sep = ',',header = 1, row.names = 1)
papi = read.table(args[3],sep = ',', header = 1, row.names = 1)
phenotype = read.table(args[4], sep = ',',header = 1, row.names = 1)

# Please promise the sample id as rownames

#pick the share samples in these datasets to correlate
co.id = intersect(rownames(species), rownames(metabolites))
co.all.id = intersect(co.id, rownames(phenotype))

co.species = species[co.all.id,]
co.metabolites = metabolites[co.all.id,]
co.papi = papi[co.all.id,]
co.phenotype = phenotype[co.all.id,]

#fillna
data = mice(co.phenotype. m = 5, method = 'pmm', maxit = 100, seed =1)
co.pheno = complete(co.phenotype)

# log transfer
log.species = log(co.species+0.01)
log.metabolites = log(co.metabolites+0.01)
log.papi = log(co.papi+0.01)
log.pheno = log(co.pheno+0.01)

# standardize log data
std.species = decostand(log.species,'standardize')
std.metabolites = decostand(log.metabolites, 'standardize')
std.papi = decostand(log.papi,'standardize')
std.pheno = decostand(log.pheno, 'standardize')

# spls between species and phenotype
species.phenotype = spls(std.species, std.pheno, K = 12, eta = 0.8)
sp = ci.spls(species.phenotype)
sp.correct =correct.spls( sp )

# spls between metabolites and phenotype
metabolites.phenotype = spls(std.metabolites,std.pheno,K =6, eta = 0.8)
mp = ci.spls(metabolites.phenotype)
mp.correct =correct.spls( mp )

# spls between species and papi pathway
species.papi = spls(std.species, std.papi, K = 12, eta = 0.8)
spapi = ci.spls(species.papi)
spapi.correct = correct.spls(spapi)

# spls between species and metabolites
species.metabolites = spls(std.species, std.metabolites,K = 12, eta = 0.8)
sm = ci.spls(species.metabolites)
sm.correct = correct.spls(sm)

# filter zero data for visualization
sp.correlated = sp.correct[which(rowSums(sp.correct)!= 0),which(colSums(sp.correct)!= 0)]
sm.correlated = sm.correct[which(rowSums(sm.correct)!=0),which(colSums(sm.correct)!=0)]
spapi.correlated = spapi.correct[which(rowSums(spapi.correct)!=0),which(colSums(spapi.correct)!=0)]
mp.correlated = mp.correct[which(rowSums(mp.correct)!= 0),which(colSums(mp.correct)!= 0)]

# set the range of the figure legend
bk_1 = c(seq(-0.6,0,by= 0.001),seq(0,0.6,by = 0.001))
bk = bk_1[!duplicated(bk_1)]

pheatmap(mp.correlated,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
          cellwidth = 30,
         cellheight =30,
         color = c(colorRampPalette(colors = c("dodgerblue","white"))(length(bk_1)/2),
                   colorRampPalette(colors = c("white","orangered"))(length(bk_1)/2)),
         legend_breaks = seq(-0.6,0.1,0.6),
         breaks = bk,
        display_numbers = TRUE,
        filename = 'metabolites.phenotype.correlation.pdf'
)

pheatmap(sp.correlated,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
          cellwidth = 30,
         cellheight =30,
         color = c(colorRampPalette(colors = c("dodgerblue","white"))(length(bk_1)/2),
                   colorRampPalette(colors = c("white","orangered"))(length(bk_1)/2)),
         legend_breaks = seq(-0.6,0.1,0.6),
         breaks = bk,
        display_numbers = TRUE,
        filename = 'species.phenotype.correlation.pdf'
)

pheatmap(sm.correlated,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
          cellwidth = 30,
         cellheight =30,
         color = c(colorRampPalette(colors = c("dodgerblue","white"))(length(bk_1)/2),
                   colorRampPalette(colors = c("white","orangered"))(length(bk_1)/2)),
         legend_breaks = seq(-0.6,0.1,0.6),
         breaks = bk,
        display_numbers = TRUE,
        filename = 'species.metabolites.correlation.pdf'
)

pdf(file = 'species.papi.correlation.pdf')
corrplot(spapi.correlated,
	tl.col = 'black'
	)
dev.off()
