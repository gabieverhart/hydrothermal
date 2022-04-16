

# read community table
com <- read.table('./data/SelectMGC_ASV_TaxAbs.txt', sep='\t', header = TRUE)
names(com) <- c(names(com)[1:9], sub('.', '', names(com)[10:24], fixed=TRUE))
# transponse community matrix
tcom <- t(com)
# subset community matrix to just site by otu with abundance
comm <- tcom[10:24, ]
comm <- apply(comm, 2, as.numeric)
colnames(comm) <- tcom[2, ]
rownames(comm) <- rownames(tcom)[10:24]
# read environmental data
dat <- read.csv('./data/MGC_Chars.csv')
# drop sites from dat we don't have data for
dat_sub <- dat[c(1, 2, 4, 5, 6, 8, 3, 7, 9:15), ]
names(dat_sub) <- c('mgc' ,'seamount', 'marker', 'duration','temp','lat', 'long', 'depth')

# example of how to reorder samples based upon a key
key <- letters
unordered_sample <- sample(key)
unordered_sample[match(key, unordered_sample)]

# multivariate analysis ------------
library(vegan)
library(mobr)

#direct ordination
cca_mod <- cca(comm ~ seamount + temp + depth, data = dat_sub)
RsquareAdj(cca_mod) 
anova(cca_mod) # full model test
anova(cca_mod, by = 'terms')
plot(cca_mod, display = c('sp', 'bp'), type='n')
text(cca_mod, display = 'cn')
text(cca_mod, display='si', col='red', cex = .5)
#points(cca_mod, display='sp', col='blue', cex =.5, pch = 19)


# indirect ordination
mds_mod <- ca(comm, trymax = 1e3)
plot(mds_mod, display = c('sp', 'si'), type = 'n')
text(mds_mod, display = 'si', col='red', cex = .5, labels = rownames(comm))
points(mds_mod, display = 'sp', col='blue', cex = .5, )

# diversity analysis ----------------
# make mob in object
mob_in <- make_mob_in(comm, dat_sub, coord_names = c('long', 'lat'), latlong = TRUE)

plot_rarefaction(mob_in, 'seamount', ref_level = 'Lōʻihi', 'SBR', lwd = 4, log = TRUE)

rarefaction(comm, method = 'IBR', effort = )

#general info about seamount

#depth,temp average Loihi, EPR
mean(dat[1:8, 8])
mean(dat[9:15, 8])
mean(dat_sub[1:8, 5])
mean(dat_sub[9:15, 8])

#graph initial trial
plot(Temperature..C. ~ Latitude, data = dat)
boxplot(1:39680,1:24)
boxplot(Temperature..C. ~ Seamount, data = dat)

#subset zeta by itself
zetas <- subset(com, subset = class=="Zetaproteobacteria")



#graph second trial
plot(Temperature..C. ~ MGC, data = dat)
boxplot(Temperature..C. ~ MGC, data = dat) 
barplot(Temperature..C. ~ MGC, data = dat, col="orange")
boxplot(temp ~ seamount, data = dat_sub, col="orange")

#trying something new to merge data for experimental reasons
plot(dat_sub[1:15, 2], type="o", col="orange")
par(new=TRUE)
plot(dat_sub[1:15, 5], type ="o", col="orange")

#rearrange to get columns and rows i DID THIS ALL MYSELF
seamount.col <- dat_sub[1:15, 2]
temp.col <- dat_sub[1:15,5]
barchart(seamount.col ~ temp.col, data= dat_sub, col="light blue", border="orange", background="gray")

#rearrange matrices
summary(comm)
summary(com)
summary(A)
summary(B)


#make my boxplot pretty
x<-dat(100) boxplot(x)
boxplot('dat',col="red")
seamount<-rnorm(100) boxplot(seamount,col="red")

#let try again
boxplot(Temperature..C. ~ Seamount, data = dat, col='red')
boxplot(Temperature..C. ~ Seamount, data = dat, col='sky blue')
rbind(comm, dat_sub)



