######################################################
#
#     Analysis of GLBRC 2008-9 bee data
#
#
#model selection experiment to determine what factors are 
#most closely associated with bee community 
#responses in bioenergy fields and landscapes and to determine
#association between bee community data and pollen deposition
#
#bee community data will be used to generate response variables, 
#site attributes will be used to create 
#predictor variables 
#
# multivariate analysis will be used to determine associations
#between nesting guilds, sociality and envuironmental/lanscape variables
###################################################




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

#compute diversity metrics, append them to overall metrics matrix
#raw diversity and richness by site
library(vegan)
metrics.matrix$H<-diversity(bee.matrix)
metrics.matrix$simp<-diversity(bee.matrix, "simpson")
metrics.matrix$S<-specnumber(bee.matrix)


#now that all the data is together and in the correct form, let's do a scatterplot matrix 
#to see if there's any obvious relationships
pairs(metrics.matrix[5:25])

#create a distance matrix so we can check for spatial autocorrelation in resiuduals in models
#Mantel test for spatial autocorrelation for each data set
#ade4 manual http://cran.r-project.org/web/packages/ade4/ade4.pdf
library(ade4)
site.dist<-dist(cbind(site.coords$LON, site.coords$LAT))




#model selection time! Same general procedure for all variables. Specify the complete model, 
#using best error structure possible for the response variable 
# (3 stepwise procedures per response variable)
# reference manual for MuMIn http://cran.r-project.org/web/packages/MuMIn/MuMIn.pdf
library(MuMIn)
library(MASS)

#model total bee abundance first, use negative binomial error structure. 
abundance.model<-glm.nb(Bee.abundance~as.factor(Year)+Treatment+Annual+Perennial+Forest+Urban+Wetland+OpenWater, na.action=na.fail, data=metrics.matrix)
summary(model.avg(dredge(abundance.model)))
#check for spatial autocorrelation in residuals of best model
pred.dist<-dist(residuals(abundance.model))
mantel.rtest(site.dist, pred.dist, nrepet =9999)
summary(abundance.model)
anova(abundance.model)
TukeyHSD(aov(abundance.model), "Treatment")

#model pollen deposition, use negative binomial error structure. 
pollen.model<-glm.nb(pollen.dep~as.factor(Year)+Treatment+Annual+Perennial+Forest+Urban+Wetland+OpenWater, na.action=na.fail, data=metrics.matrix)
summary(model.avg(dredge(pollen.model)))
#check for spatial autocorrelation in residuals of best model
pred.dist<-dist(residuals(pollen.model))
mantel.rtest(site.dist, pred.dist, nrepet =9999)
summary(pollen.model)
anova(pollen.model)
TukeyHSD(aov(pollen.model), "Treatment")

#model Apis mellifera abundance, use negative binomial error structure. 
apis.model<-glm.nb(bee.matrix$Apis.mellifera~as.factor(Year)+Treatment+Annual+Perennial+Forest+Urban+Wetland+OpenWater,  na.action=na.fail, data=metrics.matrix)
summary(model.avg(dredge(apis.model)))
#check for spatial autocorrelation in residuals of best model
pred.dist<-dist(residuals(apis.model))
mantel.rtest(site.dist, pred.dist, nrepet =9999)
summary(apis.model)
anova(apis.model)
TukeyHSD(aov(apis.model), "Treatment")

#subsequent analyses are sample size dependant so
#exclude location years with fewer than  10 bees captured
metrics.matrix.new<-metrics.matrix[which(metrics.matrix$Bee.abundance>9),] 
bee.matrix.new<-bee.matrix[which(rowSums(bee.matrix)>9),]
#compute rarefied ricness
raremax <- min(rowSums(bee.matrix.new))
metrics.matrix.new$Srare <- rarefy(bee.matrix.new, raremax)

#model Shannon's H, use normal error structure. 
H.model<-glm(H~as.factor(Year)+Treatment+Annual+Perennial+Forest+Urban+Wetland+OpenWater,  na.action=na.fail, data=metrics.matrix.new)
summary(model.avg(dredge(H.model)))
summary(H.model)
anova(H.model)
TukeyHSD(aov(H.model), "Treatment")

#model Simpson's D, use normal error structure. 
D.model<-glm(simp~as.factor(Year)+Treatment+Annual+Perennial+Forest+Urban+Wetland+OpenWater,  na.action=na.fail, data=metrics.matrix.new)
summary(model.avg(dredge(D.model)))
summary(D.model)
anova(D.model)
TukeyHSD(aov(D.model), "Treatment")

#model richness, use normal error structure
S.model<-glm(Srare~as.factor(Year)+Treatment+Annual+Perennial+Forest+Urban+Wetland+OpenWater,  na.action=na.fail, data=metrics.matrix.new)
summary(model.avg(dredge(S.model)))
summary(S.model)
anova(S.model)
TukeyHSD(aov(S.model), "Treatment")



#################################################################
#
#Multivariate analysis plotting communities, nesting guilds, sociality
#
#################################################################

# bee communities
library(vegan)

bees.no.singletons <- bee.matrix.new[, which(colSums(bee.matrix.new)>1) ]

ord.bees<-metaMDS(bees.no.singletons)
ord.bees
most_abund<-colSums(bees.no.singletons)>100
plot(ord.bees, disp='sites', type="n")
points(ord.bees, display="species", select=which(most_abund==FALSE), pch=21, cex=1, col="red")
text(ord.bees, display="species", select=which(most_abund==TRUE), cex=0.75, col="red")

#looks like the community doesn't conform to 2-D NMDS very well, but maybe functional groups of bees will work better

#set some style parameters
colvec <- c("brown", "darkgoldenrod1", "chartreuse4")
shapevec<-c(15, 16, 18)

#first by bee groups, after Winfree
ord.groups<-metaMDS(group.matrix[5:9])
ord.groups
ordiplot(ord.groups, disp='sites', type="n", xlim=c(-1.5, 1.5))
with(metrics.matrix, points(ord.nest, display = "sites", col = colvec[Treatment], pch = shapevec[Treatment], 
                            bg = colvec[Treatment]))
with(metrics.matrix, legend("topright", legend = levels(Treatment), bty = "n",
                            col = colvec, pch = shapevec[Treatment], pt.bg = colvec))
text(ord.groups, display="species", cex=1, col="red")

#CCA to overlay environmental variables on bee group ordination

ordfit.group<-envfit(ord.groups~Annual+Perennial+Forest+Urban+Wetland+OpenWater, data=metrics.matrix, perm=1000)
plot(ordfit.group, cex=1, col="blue4")
summary(ordfit.group)
ordfit.group

#now by nesting guilds


ord.nest<-metaMDS(metrics.matrix[7:12])
ord.nest
ordiplot(ord.nest, disp='sites', type="n")
with(metrics.matrix, points(ord.nest, display = "sites", col = colvec[Treatment], pch = shapevec[Treatment], 
                            bg = colvec[Treatment]))
with(metrics.matrix, legend("topright", legend = levels(Treatment), bty = "n",
                      col = colvec, pch = shapevec[Treatment], pt.bg = colvec))
text(ord.nest, display="species", cex=1, col="red")

#CCA to overlay environmental variables on bee nest ordination

ordfit.nest<-envfit(ord.nest~Annual+Perennial+Forest+Urban+Wetland+OpenWater, data=metrics.matrix, perm=1000)
plot(ordfit.nest, cex=1, col="blue4")
summary(ordfit.nest)
ordfit.nest


#now by sociality


ord.soc<-metaMDS(metrics.matrix[13:16])
ord.soc
ordiplot(ord.soc, disp='sites', type="n")
with(metrics.matrix, points(ord.nest, display = "sites", col = colvec[Treatment], pch = shapevec[Treatment], 
                            bg = colvec[Treatment]))
with(metrics.matrix, legend("topright", legend = levels(Treatment), bty = "n",
                            col = colvec, pch = shapevec[Treatment], pt.bg = colvec))
text(ord.soc, display="species", cex=1, col="red")

#CCA to overlay environmental variables on bee nest ordination

ordfit.soc<-envfit(ord.soc~Annual+Perennial+Forest+Urban+Wetland+OpenWater, data=metrics.matrix, perm=1000)
plot(ordfit.soc, cex=1, col="blue4")
summary(ordfit.soc)
ordfit.soc

