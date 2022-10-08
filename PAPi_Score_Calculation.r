# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("PAPi")

#load package
library(PAPi)

###PAPi Activity Score Calculation###
# data
args = commandArgs(T)
metabolites.with.keggid = read.table(args[1],sep = ',',header = 1, row.names = 1)

# log transformation
log.data = log(metabolites.with.keggid)

# set the group tab for the data
Replicates = data.frame(replicates = c(rep('GDM', 36), rep('NGT',46))) #gdm 36 vs. ngt 46

bind_data = cbind(Replicates, log.data)
data.for.papi = t(bind_data) # keggid as rownames, sampleid as colnames
name_data = data.frame(name = rownames(data.for.papi))
papi.df = cbind(name_data,data.for.papi)
papi.data = data.frame(papi.df,stringsAsFactors = FALSE)
print(papi.data)
papi(papi.data, save = TRUE, folder = '/home/jiating/PAPi') #folder = 'the path you want to save your data'

