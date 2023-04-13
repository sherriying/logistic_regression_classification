
library(dplyr)
library(openxlsx) 
library(ggplot2) 
library(ggsci) 
library(cowplot)
library(pROC) 
library(plyr) 
library(tibble)
library(caret)

############################################################################
############### function for predictive models: logistic regression ########
# exprssData: expression table with samples in rows and variables(e.g. proteins) in columns
# n is sample numbers in the dataset
# response variable need be changed based on your data
##############################################################################
predict.crossvalidation<-function(exprssData,n,x){
  training.samples.index<-NULL
  pred.accracy.all<-NULL
  sensitivity.all<-NULL
  specificity.all<-NULL
  glm.probes.all<-NULL
  testsample.PID.all<-NULL
  testsample.response.all<-NULL
  #predictor.x.expr.all<-NULL
  outputfile<-paste("output",".txt",sep="")
  set.seed(2)
  #set.seed(125)
  # 10-fold cross validations
  for (i in c(1:10)){
    print(sprintf("The %i time cross-validation",i))
    # taking the predictors and dependent variable (outcome:"denovo5")
    dat.sub.new<-data.frame(denovo5=exprssData[,"denovo5"],exprssData[,x])
    dat.sub.new$denovo5<-factor(dat.sub.new$denovo5,levels=c("lesion_0","lesion_5"))
    #print(dim(dat.sub.new))
    # randomly split samples into training (80%) and test samples 
    sample.index<-sample(c(1:n),size=n*0.8)
    # recording all index from 10 cross validations(bootstrap)
    training.samples.index[[i]]<-sample.index
    # taking training and test samples
    training<-dat.sub.new[sample.index,]
    testsample<-dat.sub.new[-sample.index,]
    #print(testsample[,"denovo5"])
    predictor.x.expr.all<-testsample[,-1]
    # for ROC curve
    testsample.PID.all<-append(testsample.PID.all,row.names(testsample))
    testsample.response.all<-append(testsample.response.all,testsample[,"denovo5"])
    # logistic regression model 
    glm.fit<-glm(denovo5~.,data=training,family=binomial(link="logit"))
    #print(summary(glm.fit))
    # predict on test samples and give precicted probabilities 
    glm.probes<-predict(glm.fit,testsample,type="response")
    # for ROC 
    glm.probes.all<-append(glm.probes.all,glm.probes)
    # assign phenotype "lesion_0" or "lesion_5", if glm.probes>0.5, assign to "lesion_5"
    glm.pred<-rep("lesion_0",(n-n*0.8))
    glm.pred[glm.probes>0.5]="lesion_5"
    # convert to factor
    status.code<-c("lesion_0","lesion_5")
    glm.pred<-factor(glm.pred,levels=rev(status.code))
    # real class of observation/phenotype in test samples
    actual.outcome<-factor(testsample$denovo5,levels=rev(status.code))
    # confusion table 
    conf.table<-table(glm.pred,actual.outcome)
    print(conf.table)
    sensitivity.ind<-sensitivity(conf.table)
    sensitivity.all[[i]]<-sensitivity.ind
    specificity.ind<-specificity(conf.table)
    specificity.all[[i]]<-specificity.ind
    # mean()function can be used to compute the fration of patients for which the prediction was correct
    pred.accuracy<-mean(glm.pred==testsample$denovo5)
    print(pred.accuracy)
    pred.accracy.all[[i]]<-pred.accuracy
    write(sprintf("the %i time cross validation",i),file=outputfile,append=TRUE)
    write(sprintf(
      "Confusion matrix: lesion_5  lesion_0
     lesion_5    %i        %i
     lesion_0    %i        %i
    ",conf.table[1],
      conf.table[3],
      conf.table[2],
      conf.table[4]),
      file=outputfile,append=TRUE)
    write(sprintf("predict accuracy:  %f",pred.accuracy),file=outputfile,append=TRUE)
    write(sprintf("Sensitivity:  %f",sensitivity.ind),file=outputfile,append=TRUE)
    write(sprintf("Specificity:  %f",specificity.ind),file=outputfile,append=TRUE)
    write(print("\n"),file=outputfile,append=TRUE)
  }
  # overall prediction accuracy for 10 cross-validation rounds
  pred.accuracy.mean<-mean(unlist(pred.accracy.all))
  sensitivity.mean<-mean(unlist(sensitivity.all))
  specificity.mean<-mean(unlist(specificity.all))
  write(sprintf("overall predict accuracy:  %f",pred.accuracy.mean),file=outputfile,append=TRUE)
  write(sprintf("overall sensitivity:  %f",sensitivity.mean),file=outputfile,append=TRUE)
  write(sprintf("overall specificity:  %f",specificity.mean),file=outputfile,append=TRUE)
  # ROC analysis
  # response is the real class of observation, pred is the predict probablity
  response.pred.ROC<-data.frame(PID=testsample.PID.all,response=testsample.response.all,pred=glm.probes.all)
  roc.obj<-roc(response=response.pred.ROC$response,predictor=response.pred.ROC$pred,percent=T)
  write(sprintf("AUC %s",format(roc.obj$auc)),file=outputfile,append=TRUE)
  pdf(paste("AUC",".pdf",sep=""))
  plot(roc.obj)
  dev.off()
  # output training.samples.index, convert to dataframe
  training.samples.index<-ldply(training.samples.index,rbind)
  return(list(training.samples.index=training.samples.index,response_pred=response.pred.ROC))
}

traniningSamples.index.panelGenes<-predict.crossvalidation(exprssData,65,c("FGF-21","Flt3L","CXCL9","CDCP1"))
