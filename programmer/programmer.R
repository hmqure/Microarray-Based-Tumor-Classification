library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)
library(dplyr)

#STEP 2
#Reading the CEL files

# Set working directory to CEL file locations
setwd("/projectnb/bf528/users/vangogh2022/project_1/samples")

# Save list of files as character vector
files <- list.files()

# Read CEL files
data <- ReadAffy(filenames = files)

#STEP 3
#Normalizing all of the CEL files together
eset_rma <- affy::rma(data)

#STEP 4
#Fitting the Probe Level Model 
data_fit <- fitPLM(data, background=TRUE, normalize=TRUE)

#Plotting histograms to visualize the distribution of the medians
hist(RLE(data_fit,type="stats")[1,], main="RLE - Distribution of medians",
     xlab="Median", ylab="Frequency", col="indianred2")
hist(NUSE(data_fit,type="stats")[1,], main="NUSE - Distribution of medians",
     xlab="Median", ylab="Frequency", col="lightblue")

#STEP 5
exprs = exprs(eset_rma)
#Reading the metadata
proj_metadata <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")

#Getting the batch column
batch <- proj_metadata$normalizationcombatbatch

mod <- model.matrix(~normalizationcombatmod, data=proj_metadata)

#Using combat to correct for batch effects
batch_corrected_data <- ComBat(dat = exprs, 
                               batch = batch, 
                               mod = mod)

write.csv(batch_corrected_data, file="/projectnb/bf528/users/vangogh2022/project_1/batch_corrected_data.csv")

#STEP 6
transposed_data <- t(batch_corrected_data)
scaled_data <- scale(transposed_data, center = TRUE, scale = TRUE)
result <- t(scaled_data)
prcomp_object <- prcomp(result, center = FALSE, scale = FALSE)
pcs <- data.frame(prcomp_object$rotation)
#Getting the variation
variation <- data.frame(summary(prcomp_object)$importance)
#Getting the variability explained by PC1 and PC2 and converting to percentage
variation_pc1 <- toString(variation$PC1[2]*100)
variation_pc2 <- toString(variation$PC2[2]*100)
#Adding the cancer subtype column
pcs_with_subtype <- cbind(pcs, subtype=proj_metadata$SixSubtypesClassification)

#STEP 7 
#Plotting PC1 vs PC2
pca_plot <- ggplot(pcs_with_subtype,aes(x=PC1, y=PC2,color=subtype)) + geom_point()
pca_plot + ggtitle("PCA plot") + 
  ylab(paste("PC2",variation_pc2,"%")) + xlab(paste("PC1",variation_pc1,"%"))