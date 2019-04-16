#' Robust Censored Distance Correlation Screening (RCDCS)
#'
#' @author Thomas Hielscher
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @details .
#' @return a named vector with squared distance correlation for each feature
#' @references Chen, X., Chen, X., and Wang, H. (2018). Robust feature screening for ultra-high dimensional right censored data via distance correlation. Computational Statistics & Data Analysis, 119,118-138.
#' @examples
#' \dontrun{
#' require("survival")
#' require("Rfast")
#' require("penalized")
#' data(nki70)
#' RCDCS(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @importFrom Rfast dcor
#' @export

# Chen, Chen & Wang, CSDA 2018
RCDCS <- function(df, time, status, features){
  df     <- df[order(df[,time]),]
  # Fn(t)
  kme    <- survfit(formula(paste("Surv(",time,",",status,")~1")), df)
  df$FnT <- 1 - summary(kme, censored=T, times=df[,time])$surv
  # replace xj by Fjn(x)
  df[,features] <- apply(df[,features], 2, function(d) ecdf(d)(d))
  # distance correlation w/o asymptotic results
  out <- apply(df[,features], 2, function(d) Rfast::dcor(as.matrix(d), df$FnT)$dcor^2) # ^2(?)
  return(out)
}

#' Composite Robust Censored Distance Correlation Screening (CRCDCS)
#'
#' @author Thomas Hielscher
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @param q  vector of quantiles of survival function to be evaluated, default is \code{seq(0.1,0.9,0.05)}
#' @details .
#' @return a named vector with squared distance correlation for each feature
#' @references Chen, X., Chen, X., and Wang, H. (2018). Robust feature screening for ultra-high dimensional right censored data via distance correlation. Computational Statistics & Data Analysis, 119, 118-138.
#' @examples
#' \dontrun{
#' require("survival")
#' require("Rfast")
#' require("penalized")
#' data(nki70)
#' CRCDCS(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @importFrom Rfast dcor
#' @export

# Chen, Chen & Wang, CSDA 2018
CRCDCS <- function(df, time, status, features, q=seq(0.1,0.9,0.05)){
  df     <- df[order(df[,time]),]
  # Fn(t)
  kme    <- survfit(formula(paste("Surv(",time,",",status,")~1")), df)
  df$FnT <- 1 - summary(kme, censored=T, times=df[,time])$surv
  # Qt(t)
  qt     <- quantile(kme,q)$quantile # probs apply to F(t)=1-S(t)
  # nonmissing sample distribution estimates
  q      <- q[!is.na(qt)]
  qt     <- qt[!is.na(qt)]
  # wih
  wmat <- matrix(NA, ncol=length(q),nrow=nrow(df))
  for (i in 1:length(q)) {
    wh       <- ifelse(df[,status]==1 | df$FnT > q[i], 1, (q[i] - df$FnT)/(1 - df$FnT))
    wmat[,i] <- q[i] - wh * as.numeric(df[,time] <= qt[i])
  }
  # replace xj by Fjn(x)
  df[,features] <- apply(df[,features], 2, function(d) ecdf(d)(d))
  # distance correlation
  out <- apply(df[,features], 2, function(d) Rfast::dcor(as.matrix(d), wmat)$dcor^2) # ^2(?)
  return(out)
}

#' SIS-Cox
#'
#' @author .
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @details wrapper for \code{coxph} function
#' @return a named vector with partial log-likelihood for each feature
#' @references Fan, J., Feng, Y., and Wu, Y. (2010). High-dimensional variable selection for Cox’s proportional hazards model. In Borrowing Strength: Theory Powering Applications–A Festschrift for Lawrence D. Brown, pages 70-86. Institute of Mathematical Statistics.
#' @examples
#' \dontrun{
#' require("survival")
#' require("penalized")
#' data(nki70)
#' COX(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @export

COX <- function(df, time, status, features){
  y <- Surv(df[,time], df[,status])
  out <- apply(df[,features], 2, function(d) coxph.fit(as.matrix(d), y, strata=NULL, offset=NULL, init=NULL, control=coxph.control(), weights=NULL, method="efron",
                                                       rownames=NULL)$loglik[2])
  return(out)
}

#' SIS-Cox using absolute Pearson correlation of feature values and survival time
#'
#' @author .
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @details Ignores survival status
#' @return a named vector with absolute correlation coefficient for each feature
#' @references Saldana, D. F. and Feng, Y. (2018). SIS: An R package for sure independence screening in ultrahigh-dimensional statistical models. Journal of Statistical Software, 83, 1-25.
#' @examples
#' \dontrun{
#' require("survival")
#' require("penalized")
#' data(nki70)
#' SIS(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @export

SIS <- function(df, time, status, features){
  out <- apply(df[,features], 2, function(d) abs(cor(d, df[,time])))
  return(out)
}

#' MFP based SIS-Cox
#'
#' @author .
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @param fp degrees of freedom of the FP model. Default is 4.
#' @details wrapper for \code{mfp} function
#' @return a named vector with partial log-likelihood for each feature
#' @references Ambler, G. and Benner, A. (2015). mfp: Multivariable Fractional Polynomials. R package version 1.5.2.
#' @examples
#' \dontrun{
#' require("survival")
#' require("penalized")
#' require("mfp")
#' data(nki70)
#' FP(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @export

FP <- function(df, time, status, features, fp=4){
  #  FP
  y <- Surv(df[,time], df[,status])
  out <- apply(df[,features], 2, function(d) {
    x <- try(mfp:::mfp.fit(as.matrix(d), y, cox=TRUE, gauss=FALSE, df=fp, scaling=TRUE, alpha=0.05, select=1, verbose=FALSE, xnames="x", maxits = 20, strata=NULL, offset=NULL, init=NULL, control=coxph.control(), weights=NULL, method="efron", rownames=NULL), silent=TRUE)
    if(inherits(x, "try-error")) {NA}
    else {x$loglik[2]}
  })
  return(out)
}

#' FAST (feature aberration at survival times)
#'
#' @author .
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @details only a wrapper for \code{ahaz} function in \code{ahaz} package
#' @return a named vector with absolute coefficient for each feature
#' @references Gorst-Rasmussen, A. and Scheike, T. (2013). Independent screening for single-index hazard rate models with ultrahigh dimensional features. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75, 217-245.
#' @examples
#' \dontrun{
#' require("survival")
#' require("penalized")
#' require("ahaz")
#' data(nki70)
#' FAST(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @export

# Gorst-Rasmussen & Scheike, JRSSB 2013
FAST <- function(df, time, status, features){
  uniahaz <- ahaz::ahaz(surv = Surv(df[,time],df[,status]), X = df[,features], univariate = TRUE)
  out     <- abs(uniahaz$d)
  return(out)
}

#' C-Index for each feature
#'
#' @author .
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @details only a wrapper for \code{coxph} function
#' @return a named vector with c-index for each feature
#' @examples
#' \dontrun{
#' require("survival")
#' require("penalized")
#' data(nki70)
#' CoxConc(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @export

CoxConc <- function(df, time, status, features){
  out <- apply(df[,features], 2, function(d) summary(coxph(Surv(df[,time], df[,status]) ~ d))$concordance[1])
  return(out)
}

#' LASSO via glmnet, selecting a pre-specified number of features with non-zero coefficient
#'
#' @author .
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @param dfmax dfmax parameter in call to \code{glmnet}
#' @details only a wrapper for \code{glmnet} function
#' @return a vector with absolute coefficient for each feature
#' @references Simon, N., Friedman, J., Hastie, T., and Tibshirani, R. (2011). Regularization paths for Cox's proportional hazards model via coordinate descent. Journal of Statistical Software, 39, 1-13.
#' @examples
#' \dontrun{
#' require("survival")
#' require("penalized")
#' require("glmnet")
#' data(nki70)
#' LASSO(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @export

# LASSO via glmnet
LASSO <- function(df, time, status, features, dfmax=100){
  y = Surv(df[,time], df[,status])
  fit = glmnet::glmnet(as.matrix(df[,features]), y, family="cox", standardize=FALSE, dfmax=dfmax, nlambda=200, lambda.min.ratio=0.001)
  beta <- as.vector(coef(fit, s=min(fit$lambda)))
  return(out=abs(beta))
}


#' Ball Correlation Sure Independence Screening (BCORSIS)
#'
#' @author .
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @param d see \code{bcorsis}
#' @details only a wrapper for \code{bcorsis} function in the \code{Ball} R package
#' @return a vector with Ball correlation for each feature
#' @references Pan, W., Wang, X., Xiao, W., and Zhu, H. (2018). A generic sure independence screening procedure. Journal of the American Statistical Association, to appear.
#' @examples
#' \dontrun{
#' require("survival")
#' require("penalized")
#' require("Ball")
#' data(nki70)
#' BCORSIS(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @export


BCORSIS <- function(df, time, status, features, d="large")
{
  X <- df[,features]
  out <- Ball::bcorsis(x=X, y=df[,c(time,status)], method = "survival", d=d)$complete.info[[1]]
  return(out)
}

#' Integrated powered density (IPOD)
#'
#' @author .
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @param gamma power of interest, see \code{IPOD.cont}
#' @details only a wrapper for \code{IPOD.cont} function, which is downloaded from \url{https://github.com/younghhk/software/blob/master/R/IPOD.R}
#' @return a named vector
#' @references Hong, H. G., Chen, X., Christiani, D. C., and Li, Y. (2018a). Integrated powered density: Screening ultrahigh dimensional covariates with survival outcomes. Biometrics, 74, 421-429.
#' @examples
#' \dontrun{
#' require("survival")
#' require("penalized")
#' require("survPresmooth")
#' data(nki70)
#' IPOD(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @export

IPOD <- function(df, time, status, features, gamma=1)
{
  #  IPOD(x, delta, time, gamma=1)
  X <- df[,features]
  p = ncol(X)
  cep=numeric(p)
  #source("https://github.com/younghhk/software/blob/master/R/IPOD.R", encoding="UTF-8")
  one_model = function(j) IPOD.cont(j,X,df[,status],df[,time],gamma=gamma)
  cep = sapply(1:p, one_model)
  names(cep) <- names(X)
  return(out = cep)
}

pairwise.difference <- function(m){
  npairs <- choose( ncol(m), 2 )
  results <- matrix( NA, nc=npairs, nr=nrow(m) )
  cnames <- rep(NA, npairs)
  if(is.null(colnames(m))) colnames(m) <- paste("col", 1:ncol(m), sep="")
  k <- 1
  for(i in 1:ncol(m)){
    for(j in 1:ncol(m)){
      if(j <= i) next;
      results[ ,k] <- m[ ,i] - m[ ,j]
      cnames[k] <- paste(colnames(m)[ c(i, j) ], collapse=".vs.")
      k <- k + 1
    }
  }
  colnames(results) <- cnames
  rownames(results) <- rownames(m)
  return(results)
}

IPOD.cont= function(j,x,delta,time,gamma){
  #tau=quantile(time, prob=.9)

  N=nrow(x); Lambda=(3:round(log(N),0))
  out=array(0,dim=c(length(gamma), length(Lambda)))

  for( i in 1: length(Lambda)){##i
    R=Lambda[i] #of slicings, R=3,4,5,6,...
    q=c(quantile(x[,j],probs=c((1:(R-1))/R),na.rm = TRUE))

    ##create index for subgroup
    index=rep(1,nrow(x))
    for (r in 1:length(q)){
      ind=which(x[,j]>=q[r])
      index[ind]=r+1
    }

    time.pool=numeric(R)
    for(r in 1:R){##r
      time.pool[r]=max(time[index==r])
    }
    t=seq(min(time),min(time.pool),.1)
    h=c(0,diff(range(t))/10)

    temp.a=array(0,dim=c(length(t)-1,length(gamma),R))
    R.result=numeric(length(gamma))
    for(r in 1:R){##r
      sub.time=time[index==r]; sub.delta=delta[index==r]
      f=survPresmooth::presmooth(sub.time, sub.delta, x.est=t, estimand = "f", bw.selec="fixed",fixed.bw =h)

      for(a in 1:length(gamma)){ ##a
        g=f$estimate^gamma[a]
        rec=(t[-1]-t[-length(t)])*g[-length(t)]
        temp.a[,a,r]<- cumsum(rec[1:length(rec)])
      } ##a
    }##r

    for( a in 1:length(gamma)){#a
      R.result[a]=max(apply(abs(pairwise.difference (temp.a[,a,])),1,max))
    }##a
    out[,i]=R.result
  }##i
  out=apply(out,1,sum)
  return(out)
}


#' Censored rank independence screening (CRIS)
#'
#' @author .
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @details Code downloaded from \url{https://rsong.wordpress.ncsu.edu/files/2018/09/cris.zip}
#' @return a named vector
#' @references Song, R., Lu, W., Ma, S., and Jessie Jeng, X. (2014). Censored rank independence screening for high-dimensional survival data. Biometrika, 101, 799-814.
#' @examples
#' \dontrun{
#' require("survival")
#' require("penalized")
#' data(nki70)
#' CRIS(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @export

# Song et al. Biometrika 2014
CRIS <- function(df, time, status, features){
  # Censored rank independence screening: CRIS
  n = nrow(df)
  ix = order(df[,time])
  ord.delta = df[,status][ix]
  fitc = survfit(Surv(df[,time], 1-df[,status]) ~ 1)
  Sc = fitc$surv ###estimated survival function for censoring time
  out <- apply(df[,features], 2, function(d) {
    ord.x = as.numeric(d[ix])
    rk1 = 0
    for(k in 1:(n-1)){
      if(ord.delta[k]==1){
        for(h in (k+1):n){
          rk1 = rk1 + as.numeric(ord.x[h] > ord.x[k])/(Sc[k])^2
          #          if(ord.delta[h]==1) rk2 = rk2 + as.numeric(ord.x[h] > ord.x[k])/(Sc[k]*Sc[h])
        }
      }
    }
    abs(rk1/(n*(n-1))-0.25)
    #    rk1/(n*(n-1))-0.25
  })
  #  return(sort(out, decreasing=TRUE))
  return(out)
}

#' Distance correlation based on martingale residuals (RESI)
#'
#' @author Dominic Edelmann
#' @param df data frame including columns with features and time-to-event endpoint information
#' @param time character/column name containing survival time
#' @param status  character/column name containing survival status
#' @param features  vector of column names or column indices containing features
#' @details just a wrapper for the \code{distcorr} function
#' @return a named vector of distance correlations
#' @references .
#' @examples
#' \dontrun{
#' require("survival")
#' require("penalized")
#' data(nki70)
#' RESI2(nki70,"time","event", colnames(nki70)[-c(1:7)])[1:10]
#' }
#' @export


RESI2 <- function(df, time, status, features)
{
  X <- df[,features]
  fitcox <- coxph(Surv(df[,time], df[,status]) ~ 1)
  resid <- residuals(fitcox,type="martingale")
  n=nrow(df)

  residist <- distmat(X = resid, metr.X = "euclidean", n = n, p = 1)

  dcorlist<- sapply(1:ncol(X), function(i) distcorr(X[,i], residist, type.Y = "distance"))
  names(dcorlist) <- names(X)

  return(out = dcorlist)
}
