#' Summary of IBD
#'
#' @description Summary function for IBDInfer
#' @keywords internal
#' @return No return value, called for summary.
#' @export
summary.IBD<- function(object,...){
  IBDInfer <- object
  cat(rep("_", 50), "\n")
  result<-matrix(NA, ncol=7, nrow=2)
  result <- data.frame(result)
  colnames(result)<-c("Estimate","Std.Error.wb","Std.Error.bb",
                      paste("CI.wb(",round(IBDInfer$alpha/2*100, digits=2), "%)", sep=""),
                      paste("CI.wb(",round((1-IBDInfer$alpha/2)*100, digits=2), "%)", sep=""),
                      paste("CI.bb(",round(IBDInfer$alpha/2*100, digits=2), "%)", sep=""),
                      paste("CI.bb(",round((1-IBDInfer$alpha/2)*100, digits=2), "%)", sep=""))
  rownames(result)<-c("HT","Haj")
  result[,1] <- c(IBDInfer$tau.ht,IBDInfer$tau.haj)
  result[,2] <- c(ifelse(is.na(IBDInfer$var_tau_ht_wb),NA,sqrt(IBDInfer$var_tau_ht_wb)),
                       ifelse(is.na(IBDInfer$var_tau_haj_wb),NA,sqrt(IBDInfer$var_tau_haj_wb)))
  result[,3] <- c(ifelse(is.na(IBDInfer$var_tau_ht_bb),NA,sqrt(IBDInfer$var_tau_ht_bb)),
                  ifelse(is.na(IBDInfer$var_tau_haj_bb),NA,sqrt(IBDInfer$var_tau_haj_bb)))
  result[,4] <- c(ifelse(is.na(IBDInfer$CI_ht_wb[1]),NA,IBDInfer$CI_ht_wb[1]),
                  ifelse(is.na(IBDInfer$CI_haj_wb[1]),NA,IBDInfer$CI_haj_wb[1]))
  result[,5] <- c(ifelse(is.na(IBDInfer$CI_ht_wb[2]),NA,IBDInfer$CI_ht_wb[2]),
                  ifelse(is.na(IBDInfer$CI_haj_wb[2]),NA,IBDInfer$CI_haj_wb[2]))
  result[,6] <- c(ifelse(is.na(IBDInfer$CI_ht_bb[1]),NA,IBDInfer$CI_ht_bb[1]),
                  ifelse(is.na(IBDInfer$CI_haj_bb[1]),NA,IBDInfer$CI_haj_bb[1]))
  result[,7] <- c(ifelse(is.na(IBDInfer$CI_ht_bb[2]),NA,IBDInfer$CI_ht_bb[2]),
                  ifelse(is.na(IBDInfer$CI_haj_bb[2]),NA,IBDInfer$CI_haj_bb[2]))
  print(round(result,3),right=F)
}
