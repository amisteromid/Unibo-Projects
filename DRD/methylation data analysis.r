# removing the objects and setting the working directory
rm(list=ls())
setwd("C://Users/OM3D/Desktop/DRD/REPORT")


### 1st Step; Load raw data with minfi and create an object called RGset storing the RGChannelSet object 
library(minfi)
baseDir <- ("C://Users/OM3D/Desktop/DRD/REPORT/Input_data_report/")
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets) #RGChannelSet object
save(RGset,file="RGset.RData")

### 2nd Step; Create a dataframe Red and Green to store the red and green fluorescences respectively
# GetRed and GetGreen functions to store fluorescent intensitiy values.
Red <- data.frame(getRed(RGset))
Green <- data.frame(getGreen(RGset))

### 3rd Step; what are the Red and Green fluorescences for the address assigned to you? Optional: check in the manifest file if the address corresponds to a Type I or a Type II probe and, in case of Type I probe, report its color.
load("C://Users/OM3D/Desktop/DRD/REPORT/Illumina450Manifest_clean.RData") #load Illumina_manifest450_clean file

Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="45652402",]  #checking the address probe inside the manifest file
# Type I, next Base A, Color channel Red - AddressB_ID: 64689504
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="45652402",] # No addressB_ID with this address
Illumina450Manifest_clean[Illumina450Manifest_clean$IlmnID=="cg01707559",]

Red[rownames(Red)=="45652402",]
Red[rownames(Red)=="64689504",]
# Below are noises and not main channel (green)
Green[rownames(Green)=="45652402",]
Green[rownames(Green)=="64689504",]


### 3rd Step; CREATE THE OBJECT MSet.raw
MSet.raw <- preprocessRaw(RGset)
MSet.raw # 485512 probes and 8 samples
save(MSet.raw,file="MSet_raw.RData")


### 4th Step; Perform the following quality checks and provide a brief comment to each step
qc <- getQC(MSet.raw) # consider the median of methylation and unmethylation channels for each sample
#1
plotQC(qc)
#2
controlStripPlot(RGset, controls="NEGATIVE") # the intensity values of each type of controls probes in our samples
#3
detP <- detectionP(RGset) # calculating detecion p-value
failed <- detP>0.01
summary(failed)


### 6th Step; Getting raw M value from the MSet.raw object using getM(MSet.raw) with values ranging from -infinite to + infinite 
M <- getM(MSet.raw)
Beta <- getBeta(MSet.raw)
targets$Group # to check the groups

# grouping both M and beta by their subsets
MUT_M <- M[,c(2,5,6,8)]
MUT_Beta <- Beta[,c(2,5,6,8)]
WT_M <- M[,c(1,3,4,7)]
WT_Beta <- Beta[,c(1,3,4,7)]

# Calculating the means for both M and Beta values
Beta_WT_mean <- apply(WT_Beta,1,mean,na.rm=T)
M_WT_mean <- apply(WT_M,1,mean,na.rm=T)
beta_MUT_mean <- apply(MUT_Beta,1,mean,na.rm=T)
M_MUT_mean <-apply(MUT_M,1,mean,na.rm=T)

# Removing NA values if they remained in the output vector
Beta_WT_mean <- Beta_WT_mean[!is.na(Beta_WT_mean)]
M_WT_mean <- M_WT_mean[!is.na(M_WT_mean)]
beta_MUT_mean <-beta_MUT_mean [!is.na(beta_MUT_mean)]
M_MUT_mean <- M_MUT_mean[!is.na(M_MUT_mean)]

# Calculate the density distribution for each created vector using density function
dens_mean_beta_WT <-density(Beta_WT_mean)
dens_mean_M_WT <- density(M_WT_mean)
dens_mean_beta_MUT <- density(beta_MUT_mean)
dens_mean_M_MUT <- density(M_MUT_mean )

# Plotting the density distribution for M and Beta values
par(mfrow=c(1,2))  #pairing the 2 graphs in the same image

plot(dens_mean_beta_WT,main="Density Beta values", col='Green')
lines(dens_mean_beta_MUT,col="Red")

plot(dens_mean_M_WT, main="Density M values", col='Green')
lines(dens_mean_M_MUT,col="Red")

### 7th Step;
# subsetting based on probe type
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)

# Now subsetting the beta matrix to obtain only rows whose name is in the first column of inf_I and inf_II
beta_I <- Beta[rownames(Beta) %in% dfI$IlmnID,]
beta_II <- Beta[rownames(Beta) %in% dfII$IlmnID,]

# Calculating the mean and SD
mean_of_beta_I <- apply(beta_I,1,mean)
mean_of_beta_II <- apply(beta_II,1,mean)
sd_of_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_of_beta_II <- apply(beta_II,1,sd,na.rm=T)

# preprocessQuantile Normalization
preprocessQuantile_data <- preprocessQuantile(RGset)

# Getting the beta_processQuantile matrix
beta_preprocessQuantile <- getBeta(preprocessQuantile_data)

# Dividing the previous matrix according to typeI and typeII probes
beta_preprocessQuantile_I <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfI$IlmnID,]
beta_preprocessQuantile_II <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfII$IlmnID,]
# Computing the mean and SD of the previous matrices
beta_preprocessQuantile_I_mean <- apply(beta_preprocessQuantile_I, 1, mean)
beta_preprocessQuantile_II_mean <- apply(beta_preprocessQuantile_II, 1, mean)
beta_preprocessQuantile_I_sd <- apply(beta_preprocessQuantile_I, 1, sd)
beta_preprocessQuantile_II_sd <- apply(beta_preprocessQuantile_II, 1, sd)

# Calculate the density distribution of the 2 vectors of mean and SD values:
d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)
d_sd_of_beta_I <- density(sd_of_beta_I,)
d_sd_of_beta_II <- density(sd_of_beta_II,na.rm=T)

# Plotting
par(mfrow = c(2,3)) #partition
{plot(d_mean_of_beta_I, col = "black", main = "Mean of raw Beta values", xlim = c(0,1), ylim = c(0,6), xlab = "Beta mean", ylab = "Density")
  lines(d_mean_of_beta_II, col = "#756bb1")
  legend(legend = c("TypeI", "TypeII"), col = c("black", "#756bb1"), "topright", lty = 1, bty = "n")}
{plot(d_sd_of_beta_I, col = "black", main = "SD of raw Beta values", xlim = c(0,0.5), ylim = c(0,45), xlab = "Beta SD", ylab = "Density")
  lines(d_sd_of_beta_II, col = "#756bb1")
  legend(legend = c("TypeI", "TypeII"), col = c("black", "#756bb1"), "topright", lty = 1, bty = "n")}
colorbar <- c("blue","red","green","blue","red","green","blue","red")
boxplot(Beta,ylim=c(0,1),main="Raw Beta values",las = 2, cex.axis = 0.4,col=colorbar)
{plot(density(beta_preprocessQuantile_I_mean, na.rm = T), col = "black", main = "Mean of normalized Beta values", xlim = c(0,1), ylim = c(0,4), xlab = "Beta mean", ylab = "Density")
  lines(density(beta_preprocessQuantile_II_mean, na.rm = T), col = "#756bb1")
  legend(legend = c("TypeI", "TypeII"), col = c("black", "#756bb1"), "topright", lty = 1:1, bty = "n")}
{plot(density(beta_preprocessQuantile_I_sd, na.rm = T), col = "black", main = "SD of normalize Beta values", xlim = c(0,0.5), ylim = c(0,30), xlab = "Beta SD", ylab = "Density")
  lines(density(beta_preprocessQuantile_II_sd, na.rm = T), col = "#756bb1")
  legend(legend = c("TypeI", "TypeII"), col = c("black", "#756bb1"), "topright", lty = 1:1, bty = "n")}
boxplot(beta_preprocessQuantile, ylim = c(0,1), main = "Normalized Beta values", las = 2, cex.axis = 0.4, col = colorbar)

# boxplot according to the group
colorbar <- c("blue","red","blue","blue","red","red","blue","red")
boxplot(Beta,ylim=c(0,1),main="Raw Beta values",las = 2, cex.axis = 0.4,col=colorbar)
boxplot(beta_preprocessQuantile, ylim = c(0,1), main = "Normalized Beta values", las = 2, cex.axis = 0.4, col = colorbar)

### 8th Step; PCA analysis
pca <- prcomp(beta_preprocessQuantile, scale = T)
par(mfrow = c(1,1))
{plot(pca$x[, 1], pca$x[, 2], cex = 2, pch = 18, xlab = "PC1", ylab = "PC2", xlim = c(-500, 800), main = "PCA")
  text(pca$x[, 1], pca$x[,2], labels = rownames(pca$x), pos = 4, cex = 0.5)}

library(ggplot2)
pcr <- data.frame(pca$rotation[,1:3], Group=targets$Group)
ggplot(pcr,aes(PC1,PC2, color=Group))+geom_point(size=4)
pcr <- data.frame(pca$rotation[,1:3], Group=targets$Sex)
ggplot(pcr,aes(PC1,PC2, color=Group))+geom_point(size=4)

### Step 9 Differentially methylated positions
# we create the t-test function and then apply it to whole data
t_test_func <- function(x){
  t_test <- t.test(x ~ targets$Group)
  return (t_test$p.value)
}
p_value_t_test <- apply(beta_preprocessQuantile[1:10000,], 1, t_test_func)
head(p_value_t_test, n = 10) # The first 10 entries


### 10th Step; MULTIPLE TEST CORRECTION
t_test_df <- data.frame(beta_preprocessQuantile[1:10000,], p_value_t_test)
t_test_df <- t_test_df[order(t_test_df$p_value_t_test),] # contains sorted p-value of all probes

# Bonferroni correction 
p_value_Bon <- p.adjust(t_test_df$p_value_t_test, "bonferroni")
# BH correction 
p_value_BH <- p.adjust(t_test_df$p_value_t_test, "BH")
t_test_df_corrected <- data.frame(t_test_df, p_value_Bon, p_value_BH) # We can see the corrections in this dataframe

#Number of probes with p-value less than 0.05
dim(t_test_df_corrected[t_test_df_corrected$p_value_t_test <= 0.05,])

#Number of probes that has passed the Bonferroni correction 
dim(t_test_df_corrected[t_test_df_corrected$p_value_Bon <= 0.05,])

#Number of probes that has passed the BH correction 
dim(t_test_df_corrected[t_test_df_corrected$p_value_BH <= 0.05,])


### VOLCANO PLOT and MANHATTEN PLOT
beta_corrected <- t_test_df_corrected[, 1:8]
beta_corrected_wt <- beta_corrected[, factor(targets$Group)=="WT"] #Beta matrix with p-value corrected of WTs
beta_corrected_mut <- beta_corrected[, factor(targets$Group)=="MUT"] #Beta matrix with p-value corrected of WTs
mean_beta_corrected_wt <- apply(beta_corrected_wt, 1, mean) #Mean of WTs
mean_beta_corrected_mut <- apply(beta_corrected_mut, 1, mean) #Mean of MUTs
delta <- mean_beta_corrected_wt - mean_beta_corrected_mut #Difference between the avg of both groups

# Creating a df containing the delta values and the -log10 of the p-values
volc_plot <- data.frame(delta, -log10(t_test_df_corrected$p_value_t_test))
# Volcano plot
{plot(volc_plot[, 1], volc_plot[, 2],pch = 18, cex = 0.5, xlab = "Delta", ylab = "-log10(p-value)")
  abline(h = -log10(0.05), b = 0, col = "#756bb1")
  hyper_met <- volc_plot[abs(volc_plot[, 1] > 0.2) & volc_plot[, 2] > (-log10(0.05)),]
  points(hyper_met[, 1], hyper_met[, 2], pch = 18, cex = 0.5, col = "#e34a33")
  hypo_met <- volc_plot[abs(volc_plot[, 1] < (-0.2)) & volc_plot[, 2] > (-log10(0.05)),]
  points(hypo_met[, 1], hypo_met[, 2], pch = 18, cex = 0.5, col = "#3182bd")
  legend(legend = c("Hypomethylation", "Hypermethylation"), col = c("#3182bd", "#e34a33"), "topright", bty = "n", pch = 16)}

# Converting rows with columns
t_test_df_corrected <- data.frame(rownames(t_test_df_corrected), t_test_df_corrected)
colnames(t_test_df_corrected)[1] <- "IlmnID"
t_test_df_corrected_annotated <- merge(t_test_df_corrected, Illumina450Manifest_clean, by = "IlmnID")
# Manhattan plot
library(qqman)
library(RColorBrewer)
Man_data <- t_test_df_corrected_annotated[colnames(t_test_df_corrected_annotated) %in% c("IlmnID","CHR","MAPINFO","p_value_t_test")]
order_chr <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
Man_data$CHR <- as.numeric(factor(Man_data$CHR, levels = order_chr))
manhattan(Man_data, snp = "IlmnID", chr = "CHR", bp = "MAPINFO", p = "p_value_t_test", col = brewer.pal(n = 6, "Purples"), genomewideline = F, cex = 0.5, cex.axis = 0.4)


### HEATMAP
#Creating a matrix containing our input data for the heatmaps
colorbar <- c("#377eb8", "#d95f02", "#377eb8", "#377eb8", "#d95f02", "#d95f02", "#377eb8", "#d95f02")
my_col_2 = colorRampPalette(c("#4daf4a", "black", "#e41a1c"))(100)
library(gplots)
t_test_df_corrected <- data.frame(t_test_df, p_value_Bon, p_value_BH)
heatmap_data = as.matrix(t_test_df_corrected[1:100, 1:8]) 

#Average linkage method
{heatmap.2(heatmap_data, col = my_col_2, Rowv = T, Colv = T, hclustfun = function(x) hclust(x, method = 'average'), dendrogram = "both", key = T, ColSideColors = colorbar, density.info = "none", trace = "none", scale = "none", symm = F, main = "Average linkage")
  legend(legend = c("MUT", "WT"), col = c("#d95f02", "#377eb8"), "topright", pch = 15, bty = "n", cex = 0.7)}

#Single linkage method
{heatmap.2(heatmap_data, col = my_col_2, Rowv = T, Colv = T, hclustfun = function(x) hclust(x, method = 'single'), dendrogram = "both", key = T, ColSideColors = colorbar, density.info = "none", trace = "none", scale = "none", symm = F, main = "Single linkage")
  legend(legend = c("MUT", "WT"), col = c("#d95f02", "#377eb8"), "topright", pch = 15, bty = "n", cex = 0.7)}

#Complete linkage method
{heatmap.2(heatmap_data, col = my_col_2, Rowv = T, Colv = T, dendrogram = "both", key = T, ColSideColors = colorbar, density.info = "none", trace = "none", scale = "none", symm = F, main = "Complete linkage")
  legend(legend = c("MUT", "WT"), col = c("#d95f02", "#377eb8"), "topright", pch = 15, bty = "n", cex = 0.7)}
