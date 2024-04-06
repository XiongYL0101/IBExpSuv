#' IBSurAD
#'
#' IBSurAD is used for predicting the survival rates of lung adenocarcinoma patients of stage IB
#' @param x a data.frame containing gene expression matrix of lung adenocarcinoma patients of stage IB,while each column represents a patient, and each row represents a gene
#'
#' @import rms
#' @import pec
#' @import survival

IBSurAD<-function(x){
  x<-x[as.character(coef$Gene),]
  x<-as.data.frame(apply(x,2,as.numeric))
  x<-as.data.frame(sapply(x,function(y)sum(y*(coef$coef))))
  x$Patients<-rownames(x)
  colnames(x)<-c("LAS","Patients")
  a1<-rms::cph(survival::Surv(OS.time,OS.status)~LAS,data=mm1,
               x=T, y=T, surv=T)
  t <- c(1*12,3*12,5*12)
  newd<-as.data.frame(x[,-2])
  colnames(newd)<-"LAS"
  rownames(newd)<-x$Patients
  survprob <- pec::predictSurvProb(a1,newd=newd,times=t)
  survprob1<-as.data.frame(survprob)
  colnames(survprob1)<-c("1-Year OS","3-Year OS","5-Year OS")
  survprob1
}
