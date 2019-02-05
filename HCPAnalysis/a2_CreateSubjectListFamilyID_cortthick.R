#-------------------------------------
# Create list of MZ, DZ, and singletons from HCP1200
# release with Coritcal Thickness data. 
# Also calculates basic demographic info.
#-------------------------------------
## Create subject datasets
source('~/Dropbox/SpatialACDE/Programs/Functions/createSubjData_GT.R')
load('~/Dropbox/SpatialACDE/Data/SubjectData.RData')
    
set.seed(123)
subjData1 = createSubjData_GT(subjCovData = hcpDemog2, subjImageDataFolder = '~/Dropbox/MyHCP/Data/cortThick32kLR_thickness_MSMAll/')
subjData = subjData1$subjData

# Examine genotyped twin status:
#restrict to twin pairs in the cortical thickness data:
with(subjData1$subjDataAllCovariates[subjData1$subjDataAllCovariates$nZygosity%in%c(1,2),],table(ZygositySR,ZygosityGT,useNA = 'always'))
# 31 SR DZs actually MZ

table(subjData$zygosity)
# Checks for mle dataset:
  sum(subjData$keeploglik==1)
  # 676 subjects included in likelihood calculations
  # audit  for mle: Check that there are 2 individuals in each family with twinpair and 1 otherwise
  temp = as.data.frame.matrix(with(subjData[subjData$keeploglik==1,],table(familyID,zygosity)))
  table(apply(temp,1,sum))
  # :)


# create relatedness dataset
# We create two datasets: 
#   1. subjDataForMLE: excludes siblings, for use in vertex-wise MLE and MWLE
#   2. subjDataForCovfun: includes all subjects for use in PSE-ACEM

write.csv(x = subjData, file = '~/Dropbox/SpatialACDE/Data/subjCovariatesForCovfun.csv',row.names=FALSE)
write.csv(x = subjData[subjData$keeploglik==1,], file = '~/Dropbox/SpatialACDE/Data/subjCovariatesForMLE.csv',row.names=FALSE)


# create relatedness dataset for matlab friendly objects
# We create two datasets: 
#   1. subjRelationsForMLE: excludes siblings 
#   2. subjRelationsForCovfun: includes all subjects
subjRelations = with(subjData,data.frame(subjID, familyID, zygosity,jFamilyIndex))
subjRelations$MZtp1 = with(subjRelations,ifelse(zygosity==1 & jFamilyIndex==1,1,0))
subjRelations$MZtp2 = with(subjRelations,ifelse(zygosity==1 & jFamilyIndex==2,1,0))
subjRelations$DZtp1 = with(subjRelations,ifelse(zygosity==2 & jFamilyIndex==1,1,0))
subjRelations$DZtp2 = with(subjRelations,ifelse(zygosity==2 & jFamilyIndex==2,1,0))
subjRelations$MDti = with(subjRelations,ifelse(zygosity%in%c(3,4,5),1,0))

subjRelations = subjRelations[,c('subjID','familyID','MZtp1','MZtp2','DZtp1','DZtp2','MDti')]

# make sure the ordering is the same:
all(subjRelations$subjID == subjData$subjID)

write.csv(x = subjRelations, file = '~/Dropbox/SpatialACDE/Data/subjRelatednessForCovfun.csv',row.names=FALSE)


subjRelationsMLE = subjRelations[subjData$keeploglik==1,]

write.csv(x = subjRelationsMLE, file = '~/Dropbox/SpatialACDE/Data/subjRelatednessForMLE.csv',row.names=FALSE)




#---------------------------------
# Create tables summarizing demographics for sample data:
longcov = subjData1$subjDataAllCovariates
table(longcov$Race)
table(longcov$Ethnicity)
prop.table(table(longcov$Race))

temp = longcov[!longcov$Race%in%c("Am. Indian/Alaskan Nat.","More than one","Unknown or Not Reported","Asian/Nat. Hawaiian/Othr Pacific Is."),]
temp$Race = factor(as.character(temp$Race))
(a = table(temp$Race, temp$zygosityCorrected))
prop.table(a)

prop.table(table(longcov$Ethnicity))

table(longcov$Age_in_Yrs)
table(longcov$Age)
mean(longcov$Age_in_Yrs)
sd(longcov$Age_in_Yrs)



table(longcov$Gender) #more females

table(longcov$zygosityCorrected)
a = table(longcov$Gender,longcov$zygosityCorrected)
prop.table(a)
