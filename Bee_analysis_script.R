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
#between nesting guilds, sociality and environmental/lanscape variables
###################################################




bee.specimens<-read.csv(file="C:/Rdata/GLBRC/GLBRC_bees.csv", header=TRUE, na.strings="NA")
landscape.var<-read.csv(file="C:/Rdata/GLBRC/GLBRC_bees_landscape.csv", header=TRUE, na.strings="NA")
site.coords<-read.csv(file="C:/Rdata/GLBRC/GLBRC_bees_site_coords.csv", header=TRUE, na.strings="NA")
library (reshape2)
colnames(bee.specimens)


# melt data
specimen.list<-melt(bee.specimens, id=1:16, na.rm=TRUE)

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
metrics.matrix<-merge(metrics.matrix, site.coords, by=c("Year","Site.ID", "State"))
#re-order treatment variables so they appear in order of increasing diversity
metrics.matrix$Treatment<-factor(metrics.matrix$Treatment, levels=c("corn", "switchgrass","prairie"))

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



#now that all the data is together and in the correct form, let's do a scatterplot matrix 
#to see if there's any obvious relationships
pairs(metrics.matrix[5:24])

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
site.dist.new<-dist(cbind(metrics.matrix.new$LON, metrics.matrix.new$LAT))
#compute rarefied ricness
library(vegan)
raremax <- min(rowSums(bee.matrix.new))
raremax
metrics.matrix.new$Srare <- rarefy(bee.matrix.new, raremax)
metrics.matrix.new$H<-diversity(bee.matrix.new)
metrics.matrix.new$simp<-diversity(bee.matrix.new, "simpson")

#model Shannon's H, use normal error structure. 
H.model<-glm(H~as.factor(Year)+Treatment+Annual+Perennial+Forest+Urban+Wetland+OpenWater,  na.action=na.fail, data=metrics.matrix.new)
summary(model.avg(dredge(H.model)))
#check for spatial autocorrelation in residuals of best model
pred.dist<-dist(residuals(H.model))
mantel.rtest(site.dist.new, pred.dist, nrepet =9999)
summary(H.model)
anova(H.model)
TukeyHSD(aov(H.model), "Treatment")

#model Simpson's D, use normal error structure. 
D.model<-glm(simp~as.factor(Year)+Treatment+Annual+Perennial+Forest+Urban+Wetland+OpenWater,  na.action=na.fail, data=metrics.matrix.new)
summary(model.avg(dredge(D.model)))
#check for spatial autocorrelation in residuals of best model
pred.dist<-dist(residuals(D.model))
mantel.rtest(site.dist.new, pred.dist, nrepet =9999)
summary(D.model)
anova(D.model)
TukeyHSD(aov(D.model), "Treatment")

#model richness, use normal error structure
S.model<-glm(Srare~as.factor(Year)+Treatment+Annual+Perennial+Forest+Urban+Wetland+OpenWater,  na.action=na.fail, data=metrics.matrix.new)
summary(model.avg(dredge(S.model)))
#check for spatial autocorrelation in residuals of best model
pred.dist<-dist(residuals(S.model))
mantel.rtest(site.dist.new, pred.dist, nrepet =9999)
summary(S.model)
anova(S.model)
TukeyHSD(aov(S.model), "Treatment")

library(gplots)
par(mfrow=c(2, 3), cex.axis=1.1, cex.lab=1.2, font.lab=2, cex.main=2)
plotmeans(Bee.abundance~Treatment, data=metrics.matrix, cex=2,  xlab="", ylab="Total bee abundance", n.label=FALSE, barcol="black", pch=16)
mtext("A", adj=0, font=2, cex=1.5)
plotmeans(pollen.dep~Treatment, data=metrics.matrix, cex=2,  xlab="", ylab="Pollen depostion, grains/day", n.label=FALSE, barcol="black", pch=16)
mtext("B", adj=0, font=2, cex=1.5)
plotmeans(bee.matrix$Apis.mellifera~Treatment, cex=2, data=metrics.matrix, xlab="", ylab="Honey bee abundance", n.label=FALSE, barcol="black", pch=16)
mtext("C", adj=0, font=2, cex=1.5)
plotmeans(H~Treatment, data=metrics.matrix.new, cex=2, xlab="", ylab="Shannon's H", n.label=FALSE, barcol="black", pch=16)
mtext("D", adj=0, font=2, cex=1.5)
plotmeans(simp~Treatment, data=metrics.matrix.new, cex=2,  xlab="Crop", ylab="       Simpson's D", n.label=FALSE, barcol="black", pch=16)
mtext("E", adj=0, font=2, cex=1.5)
plotmeans(Srare~Treatment, data=metrics.matrix.new, cex=2,   xlab="", ylab="Rarefied richness", n.label=FALSE, barcol="black", pch=16)
mtext("F", adj=0, font=2, cex=1.5)
#################################################################
#
#Multivariate analysis plotting communities, nesting guilds, sociality
#
#################################################################

# bee communities
library(vegan)
par(mfrow=c(2, 2), cex.axis=1.1, cex.lab=1.2, font.lab=2, cex.main=2, adj=0.5, oma = c(4, 1, 0, 1),mar = c(4, 4, 2, 4))
#set some style parameters
colvec <- c("brown", "darkgoldenrod1", "chartreuse4")
shapevec<-c(15, 17, 18)

bees.no.singletons <- bee.matrix.new[, which(colSums(bee.matrix.new)>1) ]

ord.bees<-metaMDS(bees.no.singletons)
ord.bees
most_abund<-colSums(bees.no.singletons)>250
ordiplot(ord.bees, disp='sites', type="n")
mtext("A", adj=0, font=2, cex=1.5)
with(metrics.matrix, points(ord.bees, display = "sites", col = colvec[Treatment], pch = shapevec[Treatment], 
                            bg = colvec[Treatment]))
text(ord.bees, display="species", select=which(most_abund==TRUE), cex=1, col="red")

#looks like the community doesn't conform to 2-D NMDS very well, but maybe functional groups of bees will work better



#first by bee groups, after Winfree
ord.groups<-metaMDS(group.matrix[5:9])
ord.groups
ordiplot(ord.groups, disp='sites', type="n", xlim=c(-1.5, 1.5))
mtext("B", adj=0, font=2, cex=1.5)
with(metrics.matrix, points(ord.groups, display = "sites", col = colvec[Treatment], pch = shapevec[Treatment], 
                            bg = colvec[Treatment]))
text(ord.groups, display="species", cex=1, col="red")

#CCA to overlay environmental variables on bee group ordination

ordfit.group<-envfit(ord.groups~Annual+Perennial+Forest+Urban+Wetland+OpenWater, data=metrics.matrix, perm=1000)
plot(ordfit.group, cex=1, col="blue4")
summary(ordfit.group)
ordfit.group

#now by nesting guilds


ord.nest<-metaMDS(metrics.matrix[7:12])
ord.nest
ordiplot(ord.nest, disp='sites', type="n", ylab="NMDS2", xlab="NMDS1")
mtext("C", adj=0, font=2, cex=1.5)
with(metrics.matrix, points(ord.nest, display = "sites", col = colvec[Treatment], pch = shapevec[Treatment], 
                            bg = colvec[Treatment]))
text(ord.nest, display="species", cex=1, col="red")

#CCA to overlay environmental variables on bee nest ordination

ordfit.nest<-envfit(ord.nest~Annual+Perennial+Forest+Urban+Wetland+OpenWater, data=metrics.matrix, perm=1000)
plot(ordfit.nest, cex=1, col="blue4")
summary(ordfit.nest)
ordfit.nest


#now by sociality


ord.soc<-metaMDS(metrics.matrix[13:16])
ord.soc
ordiplot(ord.soc, disp='sites', type="n", xlab="NMDS1", xlim=c(-2, 2))
mtext("D", adj=0, font=2, cex=1.5)
with(metrics.matrix, points(ord.nest, display = "sites", col = colvec[Treatment], pch = shapevec[Treatment], 
                            bg = colvec[Treatment]))
text(ord.soc, display="species", cex=1, col="red")

#CCA to overlay environmental variables on bee nest ordination

ordfit.soc<-envfit(ord.soc~Annual+Perennial+Forest+Urban+Wetland+OpenWater, data=metrics.matrix, perm=1000)
plot(ordfit.soc, cex=1, col="blue4")
summary(ordfit.soc)
ordfit.soc
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("corn", "switchgrass", "prairie"), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n", pch =shapevec, col = colvec, cex = 1.5)

