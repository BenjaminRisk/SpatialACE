#################
# This file creates subject IDs to be used in place
# of the HCP assigned subject IDs to maintain confidentiality as per
# section 4c in http://store.humanconnectome.org/data/data-use-terms/restricted-access.php/
# To obtain the seed, one needs restricted data privileges


hcpDemog = read.csv('~/Dropbox/MyHCP/Data/Untouched/RESTRICTED_ XX BLINDED XX.csv')

# extract mother IDs:
Mother_ID = hcpDemog$Mother_ID
Mother_ID = Mother_ID[!duplicated(Mother_ID)]

set.seed(AVAILABLE ONLY TO RESEARCHERS WITH ACCESS TO THE RESTRICTED DATA)
myMother_ID = data.frame(Mother_ID,'myMother_ID' = sample(1:length(Mother_ID)))

Subject = hcpDemog[,c('Subject','Mother_ID')]
sum(duplicated(Subject)) # all subject IDs are unique

myStudyID = data.frame(Subject, 'mySubject' = sample(1:nrow(Subject)))

# merge the two datasets:
myStudyID = merge(myStudyID,myMother_ID)

save(myStudyID,file='~/Dropbox/SpatialACDE/Data/StudyIDs.RData')
