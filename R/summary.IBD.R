#' Summary of IBD
#'
#' @description Summary function for IBDInfer
#' @keywords internal
#' @return No return value, called for summary.
#' @export
summary.SS<- function(object,...){
  IBDInfer <- object
  cat(rep("_", 30), "\n")
  result<-matrix(NA, ncol=5, nrow=2)
  result <- data.frame(result)
  colnames(result)<-c("Estimator","Estimate","Std.Error.wb","Std.Error.bb",
                      paste("CI.wb(",round(object$alpha/2*100, digits=2), "%)", sep=""),
                      paste("CI.wb(",round((1-object$alpha/2)*100, digits=2), "%)", sep=""),
                      paste("CI.bb(",round(object$alpha/2*100, digits=2), "%)", sep=""),
                      paste("CI.bb(",round((1-object$alpha/2)*100, digits=2), "%)", sep=""))
  rownames(result)<-""
  result[,1] <- c("HT","Haj")
  result[,2] <- c(object$tau.ht,object$tau.haj)
  result[,3] <- c(sqrt(object$var_tau_ht_wb),sqrt(object$var_tau_haj_wb))
  result[,4] <- c(sqrt(object$var_tau_ht_bb),sqrt(object$var_tau_haj_bb))
  result[,5] <- c(object$CI_ht_wb[1],object$CI_haj_wb[1])
  result[,6] <- c(object$CI_ht_wb[2],object$CI_haj_wb[2])
  result[,7] <- c(object$CI_ht_bb[1],object$CI_haj_bb[1])
  result[,8] <- c(object$CI_ht_bb[2],object$CI_haj_bb[2])
  print(round(result,3),right=F)
}
