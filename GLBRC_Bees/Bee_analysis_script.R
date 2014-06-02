bee.specimens<-read.csv(file="C:/Rdata/GLBRC/GLBRC_bees.csv", header=TRUE, na.strings="NA")
library (reshape2)
colnames(bee.specimens)


# melt data
specimen.list<-melt(bee.specimens, id=1:17, na.rm=TRUE)

#see how many species we have
species.list<-dcast(specimen.list, Taxon~Year, length)

#cast a cross-tab  by species
bee.matrix<-dcast(specimen.list, Year+Site.ID~Taxon, length)
#remove 'unknown' column- these are the specimens we were 
#unable to ID to species so they're being dropped as to not
#overinflate diversity metrics

#cast raw abundance by year and site
#creates a base data frame where we can add summary stats
#from our bee diversity matrix later
metrics.matrix<-dcast(specimen.list, Year+Site.ID+Treatment~., length)
#rename variables
names(metrics.matrx)[names(metrix.matrix)=="NA"] <- "Bee.abundance"

#cast a cross-tab with counts by pollen deposition group
group.matrix<-dcast(specimen.list, Year+Site.ID+Treatment~Group, length)

#cast a cross-tab with counts by nest biology
nest.matrix<-dcast(specimen.list, Year+Site.ID+Treatment~Nest.biology, length)

#cast a cross-tab with counts by sociality
nest.matrix<-dcast(specimen.list, Year+Site.ID+Treatment~Sociality, length)

#rename variables
#names(masterlist)[names(masterlist)=="variable"] <- "Aphid.species"
#names(masterlist)[names(masterlist)=="value"] <- "Captures"
#head(masterlist)
#row.names(masterlist)<-NULL
#export data
#write.csv(masterlist, "C:/Rdata/suction/lists/Mastersuction.csv")