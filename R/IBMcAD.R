#' IBMcAD
#'
#' IBMcAD is used to classify stage IB lung adenocarcinoma patients into high-risk (IB H), low-risk (IB L) and medium-risk (IB M) subtypes
#' @param x a data.frame containing gene expression matrix of lung adenocarcinoma patients of stage IB,while each column represents a patient, and each row represents a gene
#'
#' @import dplyr
#' @import caret


IBMcAD<-function(x){
  x1<-x[as.character(rownames(trainxtran)),]
  x2<-as.data.frame(t(x1))
  testx<-x2
  testy<-rep("1B",ncol(x1))
  pred1 <- predict(models1[[1]], newdata=cbind(testy,testx))
  pred2 <- predict(models1[[2]], newdata=cbind(testy,testx))
  pred3 <- predict(models1[[3]], newdata=cbind(testy,testx))
  pred4 <- predict(models2[[1]], newdata=cbind(testy,testx))
  pred5 <- predict(models2[[2]], newdata=cbind(testy,testx))
  pred6 <- predict(models2[[3]], newdata=cbind(testy,testx))

  tey1<-data.frame(pred1)
  tey1$Patients<-rownames(testx)
  tey1<-dplyr::rename(tey1,Type=pred1)
  tey2<-data.frame(pred2)
  tey2$Patients<-rownames(testx)
  tey2<-dplyr::rename(tey2,Type=pred2)
  tey3<-data.frame(pred3)
  tey3$Patients<-rownames(testx)
  tey3<-dplyr::rename(tey3,Type=pred3)
  tey4<-data.frame(pred4)
  tey4$Patients<-rownames(testx)
  tey4<-dplyr::rename(tey4,Type=pred4)
  tey5<-data.frame(pred5)
  tey5$Patients<-rownames(testx)
  tey5<-dplyr::rename(tey5,Type=pred5)
  tey6<-data.frame(pred6)
  tey6$Patients<-rownames(testx)
  tey6<-dplyr::rename(tey6,Type=pred6)

  H1<-subset(tey1,tey1$Type=="2H",select=c(1:2))
  H2<-subset(tey2,tey2$Type=="2H",select=c(1:2))
  H3<-subset(tey3,tey3$Type=="2H",select=c(1:2))
  H4<-subset(tey4,tey4$Type=="2H",select=c(1:2))
  H5<-subset(tey5,tey5$Type=="2H",select=c(1:2))
  H6<-subset(tey6,tey6$Type=="2H",select=c(1:2))


  inteH<-Reduce(intersect,list(H1$Patients,H2$Patients,H3$Patients,
                               H4$Patients,H5$Patients,H6$Patients))
  inteH
  H<-data.frame(inteH,rep("IB H",length(inteH)))
  colnames(H)<-c("Patients","Type")

  x3<-x[as.character(coef$Gene),]
  x3<-as.data.frame(apply(x3,2,as.numeric))
  x3<-as.data.frame(sapply(x3,function(z)sum(z*(coef$coef))))
  x3$Patients<-rownames(x3)
  colnames(x3)<-c("LAS","Patients")


  te2<-merge(x3,tey1[,c("Patients","Type")],by=c("Patients"))
  te2$Type<-1
  te2$Type[te2$Patients %in% inteH]<-"2H"
  te3<-subset(te2,te2$Type !="2H",select=c(1:ncol(te2)))
  te3$LAS1[te3$LAS>= quantile(te3$LAS, 0.319) ]<-"IB M"
  te3$LAS1[te3$LAS< quantile(te3$LAS, 0.319) ]<-"IB L"
  te4<-subset(te3,te3$LAS1=="IB M",select=c(1:ncol(te3)))
  te5<-subset(te3,te3$LAS1=="IB L",select=c(1:ncol(te3)))

  M<-dplyr::select(te4,Patients,LAS1)
  M<-dplyr::rename(M,Type=LAS1)

  L<-dplyr::select(te5,Patients,LAS1)
  L<-dplyr::rename(L,Type=LAS1)

  P<-rbind(H,M,L)

  P

}
