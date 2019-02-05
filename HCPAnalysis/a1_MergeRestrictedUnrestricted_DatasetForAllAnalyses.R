#-------------------------------------
# 
# Merge restricted and unrestricted data
# Correct some mistakes in the data
#
# NOTE: The dataset created here is used in multiple
# projects

hcpDemog = read.csv('~/Dropbox/MyHCP/Data/Untouched/RESTRICTED_XX BLINDED XX.csv')

# merge in StudyIDs to maintain confidentiality as per
# section 4c in http://store.humanconnectome.org/data/data-use-terms/restricted-access.php/
# parts of the code below use subject specific IDs to fix data errors
load("~/Dropbox/SpatialACDE/Data/StudyIDs.RData")

hcpDemog = merge(hcpDemog,myStudyID)
# Examine genotyped twin status:
table(hcpDemog$ZygosityGT,useNA = 'always')
# genetic genotyping are only available for 489 subjects
# took away three dizygotics in update on 4/12/2017


table(hcpDemog$ZygositySR,useNA = 'always')
table(hcpDemog$ZygosityGT,hcpDemog$ZygositySR)
# In all cases in which an individual reported MZ, they were MZ (226 / 226)
# However, in 70 instances (35 families), twins reported NotMZ but were MZ


##########################################
# Create variable for twin status:
hcpDemog$myZygosity = as.character(hcpDemog$ZygositySR)
# overwrite all DZs with unknown:
hcpDemog$myZygosity[hcpDemog$myZygosity%in%c('NotMZ',' ')] = 'Unknown'
# overwrite unknowns with DZs with genotype
hcpDemog$myZygosity[hcpDemog$ZygosityGT=='DZ'] = 'DZ'
# overwrite incorrect SR DZs with MZ, also overwrites two MZs without SR
hcpDemog$myZygosity[hcpDemog$ZygosityGT=='MZ'] = 'MZ'

table(hcpDemog$myZygosity,hcpDemog$ZygosityGT)

# mother <our study ID: mySubject 333> has self-reported "triplets" that are age 22, 26, and 26. The age-22
# appears to be a mistake
View(hcpDemog[hcpDemog$myMother_ID==333,])
hcpDemog$myZygosity[hcpDemog$Subject==532]='NotTwin'



# THIS CODE HAS BEEN MODIFIED THE MEET THE TERMS OF THE RESTRICTED DATA:
# HERE, DATA ARE IDENTIFIED USING CONFIDENTIAL STUDY GENERATED IDs:

# mother <our study id: myMother_ID 361> has self-reported "triplets" but are ages 28, 28, and 32
View(hcpDemog[hcpDemog$myMother_ID==361,])
# the 32-year-old does not appear to be a twin:
hcpDemog$myZygosity[hcpDemog$mySubject==1204]='NotTwin'

# ZygositySR not reported for <our study id: mySubject 534>.
# The person is a non-twin:
View(hcpDemog[hcpDemog$myMother_ID==383,])
# <our study ID: mySubject 534> is not a twin
hcpDemog$myZygosity[hcpDemog$mySubject==534]='NotTwin'

# the NotTwin subjects <my study id: mySubject  421 and 1075> look like twins -- same age. We will treat as nontwins, such that their data are used in S0 but not in S1 or S2.


##
## Inspect how many twins have their twin in the behavior dataset
## Create a summary of mother id and children (ignoring father id and 
## possiblity of half sibs)
# create a variable that counts the number of family members in each family with image data
a = data.frame(table(hcpDemog$Mother_ID))
names(a)[1] = 'Mother_ID'
names(a)[2] = 'nFamilyMembers'
# Count number of families:
nfamilies=nrow(a)
# merge number of family members:
temp = merge(hcpDemog,a)

# Create a table with columns counting the number MZ, number DZ, and number Singletons (i.e., "NotTwin") for each family:
a = as.data.frame.matrix(table(temp$Mother_ID,temp$myZygosity))
a$Mother_ID = rownames(a)
# No families have MZ and DZ:
sum(a$MZ!=0 & a$DZ!=0)
temp = merge(temp,a,by="Mother_ID")
# Create variable that indicates whether a family has a twin pair in the
# dataset:
temp$FamilyWithTwinPair = ifelse(temp$MZ==2 | temp$DZ==2,1,0)
temp$zygosityCorrected = as.character(temp$myZygosity)
temp$zygosityCorrected[temp$myZygosity%in%c('DZ','MZ') & temp$FamilyWithTwinPair==0] = 'NotTwin'

# Diagnostic table:
table(temp$zygosityCorrected,temp$myZygosity)
# This table shows that 188 / 190 of the GT DZ individuals have their twin in the dataset
# Only 336 / 364 MZs have their twin in the dataset



# Merge with the unrestricted data:
MRIheader = read.csv('~/Dropbox/MyHCP/Data/Untouched/unrestricted_ XX BLINDED XX.csv')
hcpDemog2 = merge(x = hcpDemog, y = MRIheader, by="Subject")

#--------------------------------------
# Correct some data mistakes:
# As indicated in HCP-Users, the unknown gender is a male:
# 3/31/2017: This was corrected in the 1200-subject release
# hcpDemog2$Subject[hcpDemog2$Gender=='U']
# hcpDemog2$Gender[hcpDemog2$Gender=='U'] = 'M'
# hcpDemog2$Gender = factor(as.character(hcpDemog2$Gender))
# table(hcpDemog2$Gender)


# Calculate whether any of the DZs have different genders
a = as.data.frame.matrix(table(hcpDemog2$Mother_ID[hcpDemog2$myZygosity%in%c('DZ','MZ')],hcpDemog2$Gender[hcpDemog2$myZygosity%in%c('DZ','MZ')]))
a$Mother_ID = rownames(a)
a$Flag = 0
a$Flag[a$F!=0 & a$M!=0]=1
a[a$Flag==1,]
temp = hcpDemog2[hcpDemog2$myMother_ID==361,]
temp[,c(1:5,204)]
# the 32-year-old was previously corrected to be "Not Twin", and this
# is the only family in which the gender of the twins differed

save(hcpDemog2,file='~/Dropbox/SpatialACDE/Data/SubjectData.RData')
