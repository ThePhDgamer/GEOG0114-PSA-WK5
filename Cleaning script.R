library("mice")
library("olsrr")
library("spdep")
library("spatialreg")
library("sf")
library("tmap")
library("sp")

rm(list = ls()); gc()
setwd("/Users/anwarmusah/Documents/GITHUB/GEOG0114-PSA-WK5")

temp = list.files(pattern="*.csv")
csvfiles = lapply(temp, read.csv)

head(csvfiles[[1]])
head(csvfiles[[2]])
head(csvfiles[[3]])
head(csvfiles[[4]])
head(csvfiles[[5]])
head(csvfiles[[6]])
head(csvfiles[[7]])

for (i in c(4, 5, 6)) {
  csvfiles[[3]] <- merge(csvfiles[[3]], csvfiles[[i]], by.x = "LSOACODE", by.y = "LSOACODE", all.x= TRUE)
  }

head(csvfiles[[3]])

datafile <- csvfiles[[3]]

datafile$PTACATEGORY[datafile$PTACATEGORY == "1a"] <- "1"
datafile$PTACATEGORY[datafile$PTACATEGORY == "1b"] <- "1"
datafile$PTACATEGORY[datafile$PTACATEGORY == "6a"] <- "6"
datafile$PTACATEGORY[datafile$PTACATEGORY == "6b"] <- "6"

datafile$LSOAPRICE <- as.numeric(datafile$LSOAPRICE)

datafile <- datafile[,-c(2,7)]

# merge to spatial data to get all LSOACODES
LSOAshp <- read_sf("London_LSOA_areas.shp")
names(LSOAshp)[1] <- "LSOACODE"
names(LSOAshp)[2] <- "LSOANAME"

datafile <- merge(LSOAshp, datafile, by.x = "LSOACODE", by.y = "LSOACODE", all.x = TRUE)
st_geometry(datafile) <- NULL

datafile <- datafile[,-2]

md.pattern(datafile)

datafile <- mice(datafile, m=5, maxit=50, method="pmm", seed=500)
datafile <- complete(datafile, 2)

datafile$PTACAT[datafile$PTAINDEX >= 0 & datafile$PTAINDEX < 5] <- "1. Lowest"
datafile$PTACAT[datafile$PTAINDEX >= 5 & datafile$PTAINDEX < 10] <- "2. Low"
datafile$PTACAT[datafile$PTAINDEX >= 10 & datafile$PTAINDEX < 15] <- "3. Medium"
datafile$PTACAT[datafile$PTAINDEX >= 15 & datafile$PTAINDEX < 20] <- "3. Medium"
datafile$PTACAT[datafile$PTAINDEX >= 20 & datafile$PTAINDEX < 25] <- "4. High"
datafile$PTACAT[datafile$PTAINDEX >= 25 & datafile$PTAINDEX < 101] <- "5. Highest"

names(datafile)[2] <- "AVEPRICE"
names(datafile)[3] <- "AVEINCOME"

st_write(LSOAshp, "London LSOA Areas.shp")


BOROUGHshp <- read_sf("London_Boroughs_shapefile.shp")
BOROUGHshp <- BOROUGHshp[,c(2,3)]
names(BOROUGHshp)[1] <- "BOROUGHCODE"
names(BOROUGHshp)[2] <- "BOROUGHNAME"
st_write(BOROUGHshp, "London Borough Areas.shp")

write.csv(datafile, file = "London LSOA 2015 data.csv", row.names = FALSE)
##########################################################################

#tab2
LSOAshp <- read_sf("London LSOA Areas.shp")
BOROUGHshp <- read_sf("London Borough Areas.shp")

spatialdatafile <- read.csv(file = "London LSOA 2015 data.csv", header = TRUE, sep = ",")

plot1 <- tm_shape(spatialdatafile) + tm_fill("AVEPRICE", style = "quantile", n = 7, palette = "Greens") +
	tm_shape(BOROUGHshp) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
	tm_text("BOROUGHN", size = "AREA") +
	tm_compass(position = c("right", "top")) +
	tm_scale_bar(position = c("left", "bottom")) +
	tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

plot2 <- tm_shape(spatialdatafile) + tm_fill("AVEINCOME", style = "quantile", n = 7, palette = "Greens") +
	tm_shape(BOROUGHshp) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
	tm_text("BOROUGHN", size = "AREA") +
	tm_compass(position = c("right", "top")) +
	tm_scale_bar(position = c("left", "bottom")) +
	tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

tmap_arrange(plot1, plot2, ncol = 2)
####################################

# Regression

modelMLR <- lm(log10(AVEPRICE) ~ log10(IMDSCORE) + log10(AVEINCOME) + log10(PTAINDEX), data = spatialdatafile)
options(scipen = 7)
summary(modelMLR)

spatialdatafile$RESIDUALS <- modelMLR$residuals

# ADD THESE TO TUTORIALS
ols_plot_resid_qq(modelMLR)
ols_plot_resid_fit(modelMLR)
ols_plot_resid_hist(modelMLR)

a <- as.data.frame(residuals(modelMLR))
hist(a$`residuals(modelMLR)`)


spatialdatafile <- merge(LSOAshp, spatialdatafile, by.x = "LSOACODE", by.y = "LSOACODE", all.x = TRUE)


ResidualBreaks <- c(-4,-3,-2,-1,1,2,3,4,7)

tm_shape(spatialdatafile) + tm_fill("RESIDUALS", style = "cont", midpoint = 0, palette = "-RdBu") +
	tm_shape(BOROUGHshp) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
	tm_text("BOROUGHN", size = "AREA") +
	tm_compass(position = c("right", "top")) +
	tm_scale_bar(position = c("left", "bottom")) +
	tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

#generate unique number for each row
spatialdatafile$ROWNUM <- 1:nrow(spatialdatafile)
spatialdatafile_2.0 <- as(spatialdatafile, "Spatial")
Weights <- poly2nb(spatialdatafile_2.0, row.names = spatialdatafile_2.0$ROWNUM)
WeightsMatrix <- nb2mat(Weights, style='B')
Residual_WeightMatrix <- mat2listw(WeightsMatrix , style='W')
lm.morantest(modelMLR, Residual_WeightMatrix, alternative="two.sided")


# reuse spatial weight matrix created in the object as "Residual_WeighMatrix" 
modelSLY <- lagsarlm(log10(AVEPRICE) ~ log10(IMDSCORE) + log10(AVEINCOME) + log10(PTAINDEX), data = spatialdatafile, Residual_WeightMatrix)
summary(modelSLY)

# Interprete Rho, LR and it P-value.

# extract the residuals for modelSLY object and dump back to original sf spatialdatafile object
spatialdatafile$RESID_SLY <- modelSLY$residuals
# use Moran's I test using moran.mc() function
moran.mc(spatialdatafile$RESID_SLY, Residual_WeightMatrix, 1000, zero.policy = T)
# generate the map
tm_shape(spatialdatafile) + tm_fill("RESID_SLY", style = "cont", midpoint = 0, palette = "-RdBu") +
	tm_shape(BOROUGHshp) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
	tm_text("BOROUGHN", size = "AREA") +
	tm_compass(position = c("right", "top")) +
	tm_scale_bar(position = c("left", "bottom")) +
	tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

# impacts
Weights_2.0 <- as(Residual_WeightMatrix, "CsparseMatrix")
trMC <- trW(Weights_2.0, type="MC")
impacts(modelSLY, tr = trMC, R=100)

# reuse spatial weight matrix created in the object as "Residual_WeighMatrix" 
modelSLX <- lmSLX(log10(AVEPRICE) ~ log10(IMDSCORE) + log10(AVEINCOME) + log10(PTAINDEX), data = spatialdatafile_2.0, Residual_WeightMatrix)
summary(modelSLX)

IMPACTS_SLX <- impacts(modelSLX, tr = trMC, R=100)
IMPACTS_SLX


modelSER <- errorsarlm(log10(AVEPRICE) ~ log10(IMDSCORE) + log10(AVEINCOME) + log10(PTAINDEX), data = spatialdatafile_2.0, Residual_WeightMatrix)

summary(modelSER)

# extract the residuals for modelSLY object and dump back to original sf spatialdatafile object
spatialdatafile$RESID_SER <- modelSER$residuals
# use Moran's I test using moran.mc() function
moran.mc(spatialdatafile$RESID_SER, Residual_WeightMatrix, 1000, zero.policy = T)
# generate the map
tm_shape(spatialdatafile) + tm_fill("RESID_SER", style = "cont", midpoint = 0, palette = "-RdBu") +
	tm_shape(BOROUGHshp) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
	tm_text("BOROUGHN", size = "AREA") +
	tm_compass(position = c("right", "top")) +
	tm_scale_bar(position = c("left", "bottom")) +
	tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)


lm.LMtests(modelMLR, Residual_WeightMatrix, test = c("LMerr","LMlag"))


# Generate an empty map to visualise the spatial configuration and hierarchy of LSOA and Boroughs
# First add LSOA layer 
tm_shape(LSOAshp) + tm_polygons() +
	# Add Borough layer on top of LSOA layer and make it transparent with alpha = 0
	tm_shape(BOROUGHshp) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
	# Apply cosmetics by adding compass and scale
	tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("left", "bottom"))


tm_shape(spatialdatafile) + tm_fill("RESIDUALS", style = "cont", midpoint = 0, palette = "-RdBu") +
	tm_shape(BOROUGHshp) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
	tm_text("BOROUGHN", size = "AREA") +
	tm_compass(position = c("right", "top")) +
	tm_scale_bar(position = c("left", "bottom")) +
	tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

tm_shape(spatialdatafile) + tm_fill("RESID_SLY", style = "cont", midpoint = 0, palette = "-RdBu") +
	tm_shape(BOROUGHshp) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
	tm_text("BOROUGHN", size = "AREA") +
	tm_compass(position = c("right", "top")) +
	tm_scale_bar(position = c("left", "bottom")) +
	tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)


tm_shape(spatialdatafile) + tm_fill("RESID_SER", style = "cont", midpoint = 0, palette = "-RdBu") +
	tm_shape(BOROUGHshp) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black") +
	tm_text("BOROUGHN", size = "AREA") +
	tm_compass(position = c("right", "top")) +
	tm_scale_bar(position = c("left", "bottom")) +
	tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)