# Benjamin Risk
# Create boxplots of MSE from fsem simulations from matlab

library(R.matlab)
library(ggplot2)
library(boot)
library(scales)
simresults = readMat('<MODIFY WITH YOUR DATA.mat')

mseCovfunResults = simresults$simresultsCovfunMSE
mseCovfunResults = aperm(mseCovfunResults,c(3,1,2))
niter = dim(mseCovfunResults)[1]

MISE = mseCovfunResults[,1,]
Method = c('S-SW','PSD-SW','S-FSEM','PSD-FSEM','PSD-ACE','PSD-ACE-O')
MISE = as.vector(t(MISE))
Method = rep(Method,niter)
SigmaA = data.frame(MISE,Method)
SigmaA$Method = factor(SigmaA$Method,levels=c('S-FSEM','PSD-FSEM','S-SW','PSD-SW','PSD-ACE-O', 'PSD-ACE'))

pdf(file='Boxplot_SigmaA_psd_const_sharm.pdf')#,width=10,height=7)
p <- ggplot(SigmaA, aes(x=Method,y=MISE))
p+ geom_boxplot(fill='salmon',outlier.colour='salmon')+ggtitle(expression(paste('d) ',hat(Sigma)['a'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))+ylab('ISE')
dev.off()
#theme(axis.text.x = element_text(size=rel(1.75),colour='black'), axis.title.y = element_text(size=rel(2)),axis.title.x=element_blank(), plot.title = element_text(size=rel(5)),axis.text.y = element_text(size=rel(1.75)))

MISEc = mseCovfunResults[,2,]
MISEc = as.vector(t(MISEc))
SigmaC = data.frame(MISEc,Method)
SigmaC$Method = factor(SigmaC$Method,levels=c('S-FSEM','PSD-FSEM','S-SW','PSD-SW','PSD-ACE-O', 'PSD-ACE'))
pdf(file='Boxplot_SigmaC_psd_const_sharm.pdf')
p <- ggplot(SigmaC, aes(x=Method,y=MISEc))
p+ geom_boxplot(fill='salmon',outlier.colour='salmon')+ggtitle(expression(paste('e) ',hat(Sigma)['c'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))+ylab('ISE')
dev.off()





MISEeg = mseCovfunResults[,3,]
MISEeg = as.vector(t(MISEeg))
SigmaEg = data.frame(MISEeg,Method)
SigmaEg$Method = factor(SigmaEg$Method,levels=c('S-FSEM','PSD-FSEM','S-SW','PSD-SW','PSD-ACE-O', 'PSD-ACE'))
pdf(file='Boxplot_SigmaEg_psd_const_sharm.pdf')
p <- ggplot(SigmaEg, aes(x=Method,y=MISEeg))
p+ geom_boxplot(fill='salmon',outlier.colour='salmon')+ggtitle(expression(paste('f) ',hat(Sigma)['e,G'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))+ylab('ISE')
dev.off()

#########################
# Subset to S-FSEM, PSD-FSEM, and PSD-ACE
# Scaled by 1/V^2
# These figures appear in the manuscript
SigmaA_subset = SigmaA[SigmaA$Method%in%c("S-FSEM","PSD-FSEM","PSD-ACE"),]
SigmaA_subset$MISE = SigmaA_subset$MISE/1002^2
SigmaA_subset$Method = factor(SigmaA_subset$Method,levels=c('S-FSEM','PSD-FSEM','PSD-ACE'))

pdf(file='Boxplot_SigmaA_psd_const_sharm_3methods.pdf')#,width=10,height=7)
p <- ggplot(SigmaA_subset, aes(x=Method,y=MISE))
p+ geom_boxplot(fill='salmon',outlier.colour='salmon')+ggtitle(expression(paste('d) ',hat(Sigma)['a'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))+ylab('ISE')
dev.off()



SigmaC_subset = SigmaC[SigmaC$Method%in%c("S-FSEM","PSD-FSEM","PSD-ACE"),]
SigmaC_subset$MISEc = SigmaC_subset$MISEc/1002^2
SigmaC_subset$Method = factor(SigmaC_subset$Method,levels=c('S-FSEM','PSD-FSEM','PSD-ACE'))

pdf(file='Boxplot_SigmaC_psd_const_sharm_3methods.pdf')
p <- ggplot(SigmaC_subset, aes(x=Method,y=MISEc))
p+ geom_boxplot(fill='salmon',outlier.colour='salmon')+ggtitle(expression(paste('e) ',hat(Sigma)['c'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))+ylab('ISE')
dev.off()

SigmaEg_subset = SigmaEg[SigmaEg$Method%in%c("S-FSEM","PSD-FSEM","PSD-ACE"),]
SigmaEg_subset$MISEeg = SigmaEg_subset$MISEeg/1002^2
SigmaEg_subset$Method = factor(SigmaEg_subset$Method,levels=c('S-FSEM','PSD-FSEM','PSD-ACE'))

pdf(file='Boxplot_SigmaEg_psd_const_sharm_3methods.pdf')
p <- ggplot(SigmaEg_subset, aes(x=Method,y=MISEeg))
p+ geom_boxplot(fill='salmon',outlier.colour='salmon')+ggtitle(expression(paste('f) ',hat(Sigma)['e,G'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))+ylab('ISE')
dev.off()






##############################
# Create stacked barplots()
attach(simresults)

reSortIndexSix = c(3,4,1,2,6,5)
dlabels = unlist(simresults[['depthLabels']])[reSortIndexSix]

stackbar = data.frame('Source'= c(rep('Bias-sq',6),rep('Variance',6)),'Method' = c(dlabels,dlabels),'MISE' = c(simresults[['simresultsCovMatBiassq.sumV']][1,reSortIndexSix],simresults[['simresultsCovMatVariance.sumV']][1,reSortIndexSix]))
stackbar$Method = factor(stackbar$Method,levels=c('S-FSEM','PSD-FSEM','S-SW','PSD-SW','PSD-ACE-O', 'PSD-ACE'))
pdf(file = 'Stackedbarplot_SigmaA.pdf')
ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('a) ',hat(Sigma)['a'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))
dev.off()

stackbar = data.frame('Source'= c(rep('Bias-sq',6),rep('Variance',6)),'Method' = c(dlabels,dlabels),'MISE' = c(simresults[['simresultsCovMatBiassq.sumV']][2,reSortIndexSix],simresults[['simresultsCovMatVariance.sumV']][2,reSortIndexSix]))
stackbar$Method = factor(stackbar$Method,levels=c('S-FSEM','PSD-FSEM','S-SW','PSD-SW','PSD-ACE-O', 'PSD-ACE'))
pdf(file = 'Stackedbarplot_SigmaC.pdf')
ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('b) ',hat(Sigma)['c'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))
dev.off()

stackbareg = data.frame('Source'= c(rep('Bias-sq',6),rep('Variance',6)),'Method' = c(dlabels,dlabels),'MISE' = c(simresults[['simresultsCovMatBiassq.sumV']][3,reSortIndexSix],simresults[['simresultsCovMatVariance.sumV']][3,reSortIndexSix]))
stackbareg$Method = factor(stackbareg$Method,levels=c('S-FSEM','PSD-FSEM','S-SW','PSD-SW','PSD-ACE-O', 'PSD-ACE'))
pdf(file = 'Stackedbarplot_SigmaEg.pdf')
ggplot(data=stackbareg,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('c) ',hat(Sigma)['e,G'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))
dev.off()
#########################################


###########
# re-do for three methods:
stackbar = data.frame('Source'= c(rep('Bias-sq',6),rep('Variance',6)),'Method' = c(dlabels,dlabels),'MISE' = c(simresults[['simresultsCovMatBiassq.sumV']][1,reSortIndexSix],simresults[['simresultsCovMatVariance.sumV']][1,reSortIndexSix]))
stackbar = stackbar[stackbar$Method%in%c("S-FSEM","PSD-FSEM","PSD-ACE"),]
stackbar$Method = factor(stackbar$Method,levels=c('S-FSEM','PSD-FSEM','PSD-ACE'))
stackbar$MISE = stackbar$MISE/1002^2

pdf(file = 'Stackedbarplot_SigmaA_3methods.pdf')
ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('a) ',hat(Sigma)['a'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))
dev.off()
###

stackbar = data.frame('Source'= c(rep('Bias-sq',6),rep('Variance',6)),'Method' = c(dlabels,dlabels),'MISE' = c(simresults[['simresultsCovMatBiassq.sumV']][2,reSortIndexSix],simresults[['simresultsCovMatVariance.sumV']][2,reSortIndexSix]))
stackbar = stackbar[stackbar$Method%in%c("S-FSEM","PSD-FSEM","PSD-ACE"),]
stackbar$Method = factor(stackbar$Method,levels=c('S-FSEM','PSD-FSEM','PSD-ACE'))
stackbar$MISE = stackbar$MISE/1002^2

pdf(file = 'Stackedbarplot_SigmaC_3methods.pdf')
ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('b) ',hat(Sigma)['c'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))
dev.off()


stackbar = data.frame('Source'= c(rep('Bias-sq',6),rep('Variance',6)),'Method' = c(dlabels,dlabels),'MISE' = c(simresults[['simresultsCovMatBiassq.sumV']][3,reSortIndexSix],simresults[['simresultsCovMatVariance.sumV']][3,reSortIndexSix]))
stackbar = stackbar[stackbar$Method%in%c("S-FSEM","PSD-FSEM","PSD-ACE"),]
stackbar$Method = factor(stackbar$Method,levels=c('S-FSEM','PSD-FSEM','PSD-ACE'))
stackbar$MISE = stackbar$MISE/1002^2

pdf(file = 'Stackedbarplot_SigmaEg_3methods.pdf')
ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('c) ',hat(Sigma)['e,G'])))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))
dev.off()
#########################################







######## VARIANCES AND h2


##########################################3
# Create plots for MISE of the h2
mseh2Results = simresults$simresultsh2MSE
Methodnine = unlist(simresults$depthLabels)
#Methodnine = c('S-SW','PSD-SW','S-FSEM','PSD-FSEM','PSD-ACE','PSD-ACE-O','MLE','SMLE','MWLE')
MISEh2 = as.vector(mseh2Results)
Methodnine = rep(Methodnine,niter)
h2 = data.frame(MISEh2,Methodnine)
h2$Methodnine = factor(h2$Methodnine,levels=c('MLE','MWLE','SMLE','S-FSEM','PSD-FSEM','S-SW','PSD-SW','PSD-ACE-O', 'PSD-ACE'))
pdf(file='Boxplot_h2.pdf')
p <- ggplot(h2, aes(x=Methodnine,y=MISEh2))
p+ geom_boxplot(fill='salmon',outlier.colour='salmon')+ggtitle('e) heritability')+theme(axis.text.x = element_text(size=rel(1),colour='black'), axis.title.y = element_text(size=rel(1)),axis.title.x=element_blank(), plot.title = element_text(size=rel(2)),axis.text.y = element_text(size=rel(1.75)))+ylab('ISE')+scale_y_continuous(trans=log10_trans())
dev.off()

##########
# MLE, MWLE, S-FSEM, PSD-FSEM, and PSD-ACE
h2_method5 = h2[h2$Methodnine%in%c("MLE","MWLE","S-FSEM","PSD-FSEM","PSD-ACE"),]
h2_method5$Method = factor(h2_method5$Method,levels=c('MLE','MWLE','S-FSEM','PSD-FSEM','PSD-ACE'))
h2_method5$MISEh2 = h2_method5$MISEh2/1002
pdf(file='Boxplot_h2_5methods.pdf')
p <- ggplot(h2_method5, aes(x=Method,y=MISEh2))
p+ geom_boxplot(fill='salmon',outlier.colour='salmon')+ggtitle('e) heritability')+theme(axis.text.x = element_text(size=rel(1.75),colour='black'), axis.title.y = element_text(size=rel(1)),axis.title.x=element_blank(), plot.title = element_text(size=rel(2)),axis.text.y = element_text(size=rel(1.75)))+ylab('ISE')+scale_y_continuous(trans=log10_trans())
dev.off()





#reSortIndexNine = c(7,9,8,3,4,1,2,6,5)
reSortIndexSeven = c(7,9,8,3,1,6,5)
dlabels = unlist(simresults$depthLabels)[reSortIndexSeven]
stackbar = data.frame('Source'= c(rep('Bias-sq',7),rep('Variance',7)),'Method' = c(dlabels,dlabels),'MISE' = c(simresults$simresultsVarVecBiassq.sumV[1,reSortIndexSeven],simresults$simresultsVarVecVariance.sumV[1,reSortIndexSeven]))
stackbar$Method = factor(stackbar$Method,levels=c('MLE','MWLE','SMLE','S-FSEM','S-SW','PSD-ACE-O', 'PSD-ACE'))

pdf(file = 'Stackedbarplot_sigmasqA.pdf')
ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('a) ',{hat(sigma)['a']}^2)))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))
dev.off()


#-----------------
stackbar = data.frame('Source'= c(rep('Bias-sq',7),rep('Variance',7)),'Method' = c(dlabels,dlabels),'MISE' = c(simresults$simresultsVarVecBiassq.sumV[2,reSortIndexSeven],simresults$simresultsVarVecVariance.sumV[2,reSortIndexSeven]))
stackbar$Method = factor(stackbar$Method,levels=c('MLE','MWLE','SMLE','S-FSEM','S-SW','PSD-ACE-O', 'PSD-ACE'))

pdf(file = 'Stackedbarplot_sigmasqC.pdf')
ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('b) ',{hat(sigma)^2}['c']))) +theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))+ylim(0,1.5)
dev.off()

#---------------------
# subset to 5 methods

# SigmaA
adlabels = unlist(simresults$depthLabels)
stackbar = data.frame('Source'= c(rep('Bias-sq',9),rep('Variance',9)),'Method' = c(adlabels,adlabels),'MISE' = c(simresults$simresultsVarVecBiassq.sumV[1,],simresults$simresultsVarVecVariance.sumV[1,]))

stackbar = stackbar[stackbar$Method%in%c("MLE","MWLE","S-FSEM","PSD-FSEM","PSD-ACE"),]
stackbar$Method = factor(stackbar$Method,levels=c('MLE','MWLE','S-FSEM','PSD-FSEM','PSD-ACE'))
stackbar$MISE = stackbar$MISE/1002
pdf(file = 'Stackedbarplot_sigmasqA_5methods.pdf')
ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('a) ',{hat(sigma)['a']}^2)))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1.5)),plot.title=element_text(size=rel(2)))
dev.off()

# SigmaC
stackbar = data.frame('Source'= c(rep('Bias-sq',9),rep('Variance',9)),'Method' = c(adlabels,adlabels),'MISE' = c(simresults$simresultsVarVecBiassq.sumV[2,],simresults$simresultsVarVecVariance.sumV[2,]))
stackbar = stackbar[stackbar$Method%in%c("MLE","MWLE","S-FSEM","PSD-FSEM","PSD-ACE"),]
stackbar$Method = factor(stackbar$Method,levels=c('MLE','MWLE','S-FSEM','PSD-FSEM','PSD-ACE'))
stackbar$MISE = stackbar$MISE/1002
pdf(file = 'Stackedbarplot_sigmasqC_5methods.pdf')
ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('b) ',{hat(sigma)['c']}^2)))+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1.5)),plot.title=element_text(size=rel(2)))
dev.off()
#------


adlabels = unlist(simresults$depthLabels)
 stackbar = data.frame('Source'= c(rep('Bias-sq',9),rep('Variance',9)),'Method' = c(adlabels,adlabels),'MISE' = c(simresults$simresultsVarVecBiassq.sumV[3,],simresults$simresultsVarVecVariance.sumV[3,]))
 stackbar$Method = factor(stackbar$Method,levels=c('MLE','MWLE','SMLE','S-FSEM','PSD-FSEM','S-SW','PSD-SW','PSD-ACE-O', 'PSD-ACE'))

 pdf(file = 'Stackedbarplot_sigmasqEg.pdf')
 ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('c) ',{hat(sigma)^2}['e,G']))) +theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))
 dev.off()

 # Subset to 5 methods
 stackbar = stackbar[stackbar$Method%in%c("MLE","MWLE","S-FSEM","PSD-FSEM","PSD-ACE"),]
 stackbar$Method = factor(stackbar$Method,levels=c('MLE','MWLE','S-FSEM','PSD-FSEM','PSD-ACE'))
 stackbar$MISE = stackbar$MISE/1002
 pdf(file = 'Stackedbarplot_sigmasqEg_5methods.pdf')
 ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title=expression(paste('c) ',{hat(sigma)^2}['e,G']))) +theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1.5)),plot.title=element_text(size=rel(2)))
 dev.off()
 
 
#-------------------------
stackbar = data.frame('Source'= c(rep('Bias-sq',9),rep('Variance',9)),'Method' = c(adlabels,adlabels),'MISE' = c(simresults$simresultsh2Biassq.sumV,simresults$simresultsh2Variance.sumV))
 stackbar$Method = factor(stackbar$Method,levels=c('MLE','MWLE','SMLE','S-FSEM','PSD-FSEM','S-SW','PSD-SW','PSD-ACE-O', 'PSD-ACE'))
pdf(file = 'Stackedbarplot_h2.pdf')
ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title='d) heritability')+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1)),plot.title=element_text(size=rel(2)))
dev.off()

#-------
# Subset to 5 methods
stackbar = stackbar[stackbar$Method%in%c("MLE","MWLE","S-FSEM","PSD-FSEM","PSD-ACE"),]
stackbar$Method = factor(stackbar$Method,levels=c('MLE','MWLE','S-FSEM','PSD-FSEM','PSD-ACE'))
stackbar$MISE = stackbar$MISE / 1002
pdf(file = 'Stackedbarplot_h2_5methods.pdf')
ggplot(data=stackbar,aes(x=Method,y=MISE,fill=Source))+geom_bar(stat='identity')+labs(title='d) heritability')+theme(axis.text.y = element_text(size=rel(1.75)),axis.text.x = element_text(size=rel(1.4)),plot.title=element_text(size=rel(2)))
dev.off()

