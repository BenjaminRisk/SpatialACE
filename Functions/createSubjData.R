
createSubjData = function(subjCovData,subjImageDataFolder) {
  
  # list subjects with image data:
  subjImgList = dir(subjImageDataFolder)
  
  subjImgList = as.numeric(substr(subjImgList,1,6))
  
  
  subjImgData = data.frame('Subject'=subjImgList,'ImgData'=1)
  
  # merge two datasets; keep only subjects with CortThick data:
  hcp = merge(subjCovData,subjImgData,by = 'Subject')

  
  # create a variable that counts the number of family members in each family with image data
  a = data.frame(table(hcp$Mother_ID))
  names(a)[1] = 'Mother_ID'
  names(a)[2] = 'nFamilyMembers'
  # Count number of families:
  nfamilies=nrow(a)
  
  hcpMerge2 = merge(hcp,a)
  
  # Create a table with columns counting the number MZ, number DZ (i.e., "NotMZ") and number Singletons (i.e., "NotTwin") for each family:
  a = as.data.frame.matrix(table(hcpMerge2$Mother_ID,hcpMerge2$Zygosity))
  a$Mother_ID = rownames(a)

# checks when debugging:  
      max(a$MZ)
      max(a$NotMZ)
      
      # No families have MZ and DZ:
      sum(a$MZ!=0 & a$NotMZ!=0)
  
  hcpMerge3 = merge(hcpMerge2,a,by="Mother_ID")
  
# checks when debugging: the nFamilyMembers is equivalent to MZ
      all((hcpMerge3$MZ+hcpMerge3$NotMZ+hcpMerge3$NotTwin) == hcpMerge3$nFamilyMembers)
  
  # Create variable that indicates whether a family has a twin pair in the
  # dataset:
  hcpMerge3$FamilyWithTwinPair = ifelse(hcpMerge3$MZ==2 | hcpMerge3$NotMZ==2,1,0)
  
  hcpMerge3$keep=0
  
  # Keep MZ twin pairs:
  hcpMerge3$keep[hcpMerge3$MZ==2 & hcpMerge3$Zygosity=='MZ'] = 1
  # sets keep to 1 if both children in an MZ pair appear in the dataset
  
  # Keep DZ twin pairs:
  hcpMerge3$keep[hcpMerge3$NotMZ==2 & hcpMerge3$Zygosity=='NotMZ'] = 1
  # sets keep to 1 if both children in a DZ pair appear in the dataset
  
  # keep individuals in families with only one member:
  hcpMerge3$keep[hcpMerge3$nFamilyMembers==1] = 1
  
  
  # Randomly assign the j index within each family.
  familyIDs = hcpMerge3$Mother_ID
  familyIDs = familyIDs[!duplicated(familyIDs)]

  hcpMerge3$jFamilyIndex=0
  for (k in familyIDs) {
    nFamilyMembers = hcpMerge3$nFamilyMembers[hcpMerge3$Mother_ID==k][1]
    hcpMerge3$jFamilyIndex[hcpMerge3$Mother_ID==k] = sample(1:nFamilyMembers)
  }  
  
  
  
  # set keep=1 for all jFamilyIndex equal to 1,
  # then set keep=0 if jFamilyIndex==1, Twin_Stat==NotTwin, and FamilyWithTwinPair==1
  hcpMerge3$keep[hcpMerge3$jFamilyIndex==1]=1
  hcpMerge3$keep[hcpMerge3$Twin_Stat=='NotTwin' & hcpMerge3$FamilyWithTwinPair==1] = 0
  hcpMerge4 = hcpMerge3[hcpMerge3$keep==1,]
  
# check when debugging: number of families should equal previous total:
      length(table(hcpMerge4$Mother_ID))==nfamilies
  
  #----------------------------
  # Re-define zygosity to be specific to the sample with cort thick data:
  # in the variable "zygosityCorrected", the value is "NotTwin" if a twin
  # singleton appears in the dataset (i.e., no data for his/her twin).
  hcpMerge5 = hcpMerge4
  hcpMerge5$zygosityCorrected = as.character(hcpMerge5$Zygosity)
  hcpMerge5$zygosityCorrected[hcpMerge5$Twin_Stat=='Twin' & hcpMerge5$FamilyWithTwinPair==0] = 'NotTwin'
  hcpMerge5$zygosityCorrected[hcpMerge5$zygosityCorrected=='NotMZ'] = 'DZ'
  

  # create a numeric vector in which 1 = MZ, 2 = DZ, and 3 = Singleton
  hcpMerge5$nZygosity = 0
  hcpMerge5$nZygosity[hcpMerge5$zygosityCorrected=='MZ'] = 1 
  hcpMerge5$nZygosity[hcpMerge5$zygosityCorrected=='DZ'] = 2 
  hcpMerge5$nZygosity[hcpMerge5$zygosityCorrected=='NotTwin'] = 3 
  hcpMerge5$dGender = ifelse(hcpMerge5$Gender=='F',1,0)
  
  # dMRI have all been updated with  r227
  #dMRI = rep(NA,nrow(hcpMerge5))
  #dMRI[hcpMerge5$dMRI_3T_ReconVrs=='r177']=0
  #dMRI[hcpMerge5$dMRI_3T_ReconVrs=='r227']=1
  
  # there are 6 subjects from Q 3 that are labeled "r177 r227", so presumably 
  # some of their runs were 177 and some 227. Make these 177 so we can
  # use their data in matlab:
  fMRI = numeric(nrow(hcpMerge5))
  fMRI[hcpMerge5$fMRI_3T_ReconVrs=='r227']=1
  
  
  subjData = with(hcpMerge5,data.frame('subjID' = Subject, 'familyID' = Mother_ID, 'zygosity'= nZygosity, 'jFamilyIndex'=jFamilyIndex, dGender, 'Age' = scale(Age_in_Yrs), 'Handedness' = scale(Handedness), 'Height' = scale(Height), 'Weight' = scale(Weight), 'BMI' = scale(BMI),'FS_IntraCranial_Vol' = scale(FS_IntraCranial_Vol),'fMRI_3T_ReconVrs'=fMRI))
  
  
  # when there is missing data, matlab programs don't work.
     sum(is.na(subjData$Height))
  # only 1 subject is missing Height, etc, data
  
  # Replace "NA" with 0, since all variables are centered and
  # this is equivalent to imputing the mean
  subjData$Height[is.na(subjData$Height)] = 0
  subjData$Weight[is.na(subjData$Weight)] = 0
  subjData$BMI[is.na(subjData$BMI)] = 0 
  
  
  # over-write the jFamilyIndex:
  familyIDs = subjData$familyID[!duplicated(subjData$familyID)]
  for (k in familyIDs) {
    nmembers = sum(subjData$familyID==k)
    subjData$jFamilyIndex[subjData$familyID==k] = 1:nmembers
  }
  
  orderUp = with(subjData,order(zygosity,familyID,jFamilyIndex))
  subjData = subjData[orderUp,]
  hcpMerge5 = hcpMerge5[orderUp,]
  
  # check subjData and hcpMerge5 have same sort order:
  all(subjData$subjID == hcpMerge5$Subject)
  
  return(list('subjData'=subjData,'subjDataAllCovariates'=  hcpMerge5[,names(hcpMerge5)!='jFamilyIndex']))
}