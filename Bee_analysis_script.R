bee.specimens<-read.csv(file="C:/Rdata/GLBRC/GLBRC_bees.csv", header=TRUE, na.strings="NA")
landscape.var<-read.csv(file="C:/Rdata/GLBRC/GLBRC_bees_landscape.csv", header=TRUE, na.strings="NA")
site.coords<-read.csv(file="C:/Rdata/GLBRC/GLBRC_bees_site_coords.csv", header=TRUE, na.strings="NA")
library (reshape2)
colnames(bee.specimens)


# melt data
specimen.list<-melt(bee.specimens, id=1:17, na.rm=TRUE)

#see how many species we have
species.list<-dcast(specimen.list, Taxon~Year, length)
species.list.bystate<-dcast(specimen.list, Taxon~State, length)

#cast raw abundance by year and site
#creates a base data frame where we can add summary stats
#from our bee diversity matrix later
metrics.matrix<-dcast(specimen.list, Year+State+Site.ID+Treatment~., length)
#rename variables
names(metrics.matrix)[id=5] <- "Bee.abundance"

#cast a cross-tab with counts by pollen deposition group
group.matrix<-dcast(specimen.list, Year+State+Site.ID+Treatment~Group, length)
#develop total pollen deposition estimate based on winfree 2007 simulations
#and put it in the oerall metrics matrix
metrics.matrix$pollen.dep<-4335*group.matrix$Honey+2814*group.matrix$Bumble+1515*group.matrix$Small+624*group.matrix$Green+94*group.matrix$Large

#cast a cross-tab with counts by nest biology
nest.matrix<-dcast(specimen.list, Year+State+Site.ID+Treatment~Nest.biology, length)
# Quick PCA to reduce dimensionality of nesting biology data
nest.pca<-prcomp(nest.matrix[5:10], scale.=T)
biplot(nest.pca)
summary(nest.pca)
#nest biology dimensionality does not reduce well. Guess each guild needs to be its own response variable
#Do a pearson correlation analysis to see if we can at least pull out individual correlations
library(Hmisc)
rcorr(as.matrix(nest.matrix[5:10]), type="pearson")

#cast a cross-tab with counts by sociality
sociality.matrix<-dcast(specimen.list, Year+State+Site.ID+Treatment~Sociality, length)
sociality.pca<-prcomp(sociality.matrix[5:8], scale.=T)
biplot(sociality.pca)
summary(sociality.pca)

#socialbiology dimensionality does not reduce well either. As above.
#Do a pearson correlation analysis as above
rcorr(as.matrix(sociality.matrix[5:8]), type="pearson")

#do a pca on landscape variables- maybe we can reduce that?
landscape.pca<-prcomp(landscape.var[3:8], scale.=T)
biplot(landscape.pca)
summary(landscape.pca)

#aggregate all these data into the metrix.matrix table
metrics.matrix<-merge(metrics.matrix, nest.matrix, by=c("Year","State", "Site.ID", "Treatment"))
metrics.matrix<-merge(metrics.matrix, sociality.matrix, by=c("Year","State", "Site.ID", "Treatment"))
metrics.matrix<-merge(metrics.matrix, landscape.var, by=c("Year","Site.ID"))

#cast a cross-tab  by species
bee.matrix<-dcast(specimen.list, Year+State+Site.ID~Taxon, length)
#remove 'unknown' column- these are the specimens we were 
#unable to ID to species so they're being dropped as to not
#overinflate diversity metrics
bee.matrix$Unknown<-NULL
#also  remove site, year columns from data (required for diversity computations)
bee.matrix$Year<-NULL
bee.matrix$State<-NULL
bee.matrix$Site.ID<-NULL
head(bee.matrix)

#compute diversity metrics, append them to overall metrics matrix
#raw diversity and richness by site
library(vegan)
metrics.matrix$H<-diversity(bee.matrix)
metrics.matrix$simp<-diversity(bee.matrix, "simpson")
metrics.matrix$S<-specnumber(bee.matrix)
#rarified richness
raremax <- min(rowSums(bee.matrix))
metrics.matrix$Srare <- rarefy(bee.matrix, raremax)

#now that all the data is together and in the correct form, let's do a scatterplot matrix 
#to see if there's any obvious relationships
pairs(metrics.matrix[5:26])

#create a distance matrix so we can check for spatial autocorrelation in resiuduals in models
#Mantel test for spatial autocorrelation for each data set
#ade4 manual http://cran.r-project.org/web/packages/ade4/ade4.pdf
library(ade4)
site.dist<-dist(cbind(site.coords$LON, site.coords$LAT))




#model selection time! Same general procedure for all variables. Specify the complete model, 
#using best error structure possible for the response variable 
# (3 stepwise procedures per response variable)
# reference manual for MuMIn http://cran.r-project.org/web/packages/MuMIn/MuMIn.pdf
#note forest and grassland percentages dropped from analysis because 
# both wer at very low percentages in the landscape
library(MuMIn)
library(MASS)

#model total bee abundance first, use negative binomial error structure. 
abundance.model<-glm.nb(Bee.abundance~as.factor(Year)+State+Treatment+Annual+Perennial+Forest+Urban+Wetland+OpenWater, data=metrics.matrix)
dredge(abundance.model, extra="R^2")
summary(model.avg(dredge(abundance.model)))
#check for spatial autocorrelation in residuals of best model
pred.dist<-dist(residuals(abundance.model))
mantel.rtest(site.dist, pred.dist, nrepet =9999)
anova(abundance.model)
