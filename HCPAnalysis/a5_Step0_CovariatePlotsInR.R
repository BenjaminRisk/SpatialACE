# Create plots for web supplement from 
# cortical thickness analysis
library(ggplot2)
library(gridExtra)
library(R.matlab)

# create histograms of pvalues used in variable selection
mlevarselect = readMat('~/Dropbox/SpatialACDE/Results/Results_1200SubjectRelease/mle_cortThick676_SelectCovariates.mat')
# takes a few minutes to read

allpvalues = data.frame(t(mlevarselect$mleresultsR.HIGH[[9]]))
names(allpvalues) = c('intercept','gender','age','handedness','height','weight','BMI','ICV')

# overlay histogram and normal density
b11 = ggplot(allpvalues, aes(intercept)) + geom_histogram(aes(y = ..density..),binwidth=0.05) + stat_function(fun = dunif, args = list(min = 0, max=1), lwd = 1, col = 'red')+xlab('Intercept')

b12 = ggplot(allpvalues, aes(gender)) + geom_histogram(aes(y = ..density..),binwidth=0.05) +
  stat_function(fun = dunif, args = list(min=0,max=1), lwd = 1, col = 'red')+xlab('Gender')

b13 = ggplot(allpvalues, aes(age)) + geom_histogram(aes(y = ..density..),binwidth=0.05) +
  stat_function(fun = dunif, args = list(min=0,max=1), lwd = 1, col = 'red')+xlab('Age')

b14 = ggplot(allpvalues, aes(handedness)) + geom_histogram(aes(y = ..density..),binwidth=0.05) +
  stat_function(fun = dunif, args = list(min=0,max=1), lwd = 1, col = 'red')+xlab('Handedness')


b21 = ggplot(allpvalues, aes(height)) + geom_histogram(aes(y = ..density..),binwidth=0.05) +
  stat_function(fun = dunif, args = list(min=0,max=1), lwd = 1, col = 'red')+xlab('Height')

b22 = ggplot(allpvalues, aes(weight)) + geom_histogram(aes(y = ..density..),binwidth=0.05) +
  stat_function(fun = dunif, args = list(min=0,max=1), lwd = 1, col = 'red')+xlab('Weight')

b23 = ggplot(allpvalues, aes(BMI)) + geom_histogram(aes(y = ..density..),binwidth=0.05) +
  stat_function(fun = dunif, args = list(min=0,max=1), lwd = 1, col = 'red')+xlab('BMI')

b24 = ggplot(allpvalues, aes(ICV)) + geom_histogram(aes(y = ..density..),binwidth=0.05) +
  stat_function(fun = dunif, args = list(min=0,max=1), lwd = 1, col = 'red')+xlab('ICV')

pdf(file='~/Dropbox/SpatialACDE/Documents/Figures/Figures_1200SubjectRelease/Examine_pvalues_allcovariates_cortThick.pdf',width = 10,height=5)
grid.arrange(b11,b12,b13,b14,b21,b22,b23,b24,nrow=2,ncol=4)
dev.off()




#------------------------------------
# OLD CODE
# look at fixed effects
library(R.matlab)
cortthick = readMat('~/Dropbox/SpatialACDE/Data/subjectData_cortThick.mat')
cortdata = cortthick$dataMatR
#covnames = {'intercept','gender','age','handedness','height','weight','BMI','ICV'};

outvertex = lm(cortthick$dataMatR[199,]~subjData$dGender+subjData$Age+subjData$Handedness+subjData$Height+subjData$Weight+subjData$BMI+subjData$FS_IntraCranial_Vol)
summary(outvertex)
plot(outvertex)
# there is a bimodel distribution in the predicted values due to gender in some subset of the locations.
# this bimodal distribution is not present in the original data
hist(predict(outvertex))
hist(cortthick$dataMatR[199,])
