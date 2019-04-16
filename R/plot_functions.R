#' performance plot, visualization of performance vs model size
#'
#' @author Manuela Hummel
#' @param scores scores matrix; output of screening methods for simulated data, larger value means better ranking, larger value means better ranking
#' @param trueVars vector of indices of true variables
#' @param what performance criterion to be plotted (TPR, TDR, FDR, propallselected)
#' @param minm minimal model size
#' @param maxm maximal model size
#' @param col colour
#' @param lwd lwd
#' @param ylab ylab
#' @param ylimit limit of y-axis
#' @param add add to existing plot
#' @param ... further parameters pased to plot function
#' @details y = average TPR/TDR/FDR/propallselected, x = respective model size. vertical line indicates MMMS (median of minimum model size needed to detect all true variables). can be shown for all variables (m = 1,...,p) or up to 'maxm' variables (m = 1,...,maxm). show e.g. TPR and FDR for one method or compare e.g. TPR across methods
#' @export

perfplot <- function(scores, trueVars, what=c("TPR", "TDR", "FDR", "propallselected"), minm=25, maxm, col=1, lwd=2, ylab, ylimit=c(0,1), add=FALSE, ...){
  # add: if line should be added to existing plot
  # ...: arguments passed to plot() or lines()

  what <- match.arg(what)
  p <- ncol(scores)
  if(missing(maxm))
    maxm <- p

  # median of minimum model size needed to detect all true variables
  MMMS <- mmms(scores, trueVars)

  # performance
  perf <- switch(what,
                 TPR = TPR(scores, minm:maxm, trueVars),
                 TDR = 1 - FDR(scores, minm:maxm, trueVars),
                 FDR = FDR(scores, minm:maxm, trueVars),
                 propallselected = propallselected(scores, minm:maxm, trueVars))

  if(missing(ylab))
    ylab <- what

  if(!add)
    plot(minm:maxm, perf, type="l", xlab="Model size", ylab=ylab, ylim=ylimit, xaxt="n", col=col, lwd=lwd, ...)
  else
    lines(minm:maxm, perf, col=col, lwd=lwd, ...)
  axis(side=1, at=seq(minm, maxm, length.out=4))
  abline(v=MMMS, col=col)
}

#' performance plot, visualization of performance vs model size for several methods
#'
#' @author Manuela Hummel
#' @param files paths + names of .RData files containing simulation results
#' @param trueVars vector of indices of true variables
#' @param what performance criterion to be plotted (TPR, TDR, FDR, propallselected)
#' @param whichMethods methods to be plotted
#' @param maxm maximal model size
#' @param ylimit limit of y-axis
#' @details .
#' @export

plotall.perf <- function(files, trueVars, whichMethods, what=c("TPR", "TDR", "FDR", "propallselected"), ylimit=c(0,1), maxm=100){
  # files: paths + names of .RData files containing simulation results
  what <- match.arg(what)
  for(l in files){
    sim <- sapply(strsplit(l, "/"), function(x) x[length(x)])
    sim <- sapply(strsplit(sim, "\\."), function(x) x[1])
    print(l)
    load(l)
    if(!missing(whichMethods))
      output <- output[whichMethods,,]
    n.methods <- dim(output)[1]
    N <- dimnames(output)[[1]]
    lty <- rep(1:3, length.out=n.methods)
    par(mar=c(4,4,4,8))
    perfplot(output[1,,], trueVars, what=what, maxm=maxm, ylimit=ylimit, main=sim)
    for(i in 2:n.methods)
      perfplot(output[i,,], trueVars, what=what, maxm=maxm, ylimit=ylimit, col=i, lty=lty[i], add=TRUE)
    legend(par("usr")[2], par("usr")[4], N, col=1:n.methods, lwd=2, lty=lty, xpd=TRUE)
  }
}


#' performance plot, visualization of performance vs model size
#'
#' @author Manuela Hummel
#' @param scores scores matrix; output of screening methods for simulated data, larger value means better ranking, larger value means better ranking
#' @param trueVars vector of indices of true variables
#' @param what performance criterion to be plotted (TPR, TDR, FDR, propallselected)
#' @param minm minimal model size
#' @param maxm maximal model size
#' @param col colour
#' @param lwd lwd
#' @param ylab ylab
#' @param ylimit limit of y-axis
#' @param add add to existing plot
#' @param ... further parameters pased to plot function
#' @details y = average TPR/TDR/FDR/propallselected, x = respective model size. vertical line indicates MMMS (median of minimum model size needed to detect all true variables). can be shown for all variables (m = 1,...,p) or up to 'maxm' variables (m = 1,...,maxm). show e.g. TPR and FDR for one method or compare e.g. TPR across methods
#' @export

my.perfplot <- function(scores, trueVars, what=c("TPR", "TDR", "FDR", "propallselected"), minm=25, maxm, col=1, lwd=2, ylab, ylimit=c(0,1), add=FALSE, ...){

  what <- match.arg(what)
  p <- ncol(scores)
  if(missing(maxm))
    maxm <- p

  # median of minimum model size needed to detect all true variables
  MMMS <- mmms(scores, trueVars)

  # performance
  perf <- switch(what,
                 TPR = TPR(scores, minm:maxm, trueVars),
                 TDR = 1 - FDR(scores, minm:maxm, trueVars),
                 FDR = FDR(scores, minm:maxm, trueVars),
                 propallselected = propallselected(scores, minm:maxm, trueVars))

  if(missing(ylab))
    ylab <- what

  if(!add)
    plot(minm:maxm, perf, type="l", xlab="Model size", ylab=ylab, ylim=ylimit, col=col, lwd=lwd, xaxt="n", ...)
  else
    lines(minm:maxm, perf, col=col, lwd=lwd, ...)
  axis(side=1, at=seq(minm, maxm, length.out=4))
  abline(v=MMMS, col=col)
}

#' performance plot, visualization of performance vs model size for several methods
#'
#' @author Manuela Hummel
#' @param files paths + names of .RData files containing simulation results
#' @param trueVars vector of indices of true variables
#' @param what performance criterion to be plotted (TPR, TDR, FDR, propallselected)
#' @param whichMethods methods to be plotted
#' @param maxm maximal model size
#' @param ylimit limit of y-axis
#' @details .
#' @export

my.plotall.perf <- function(files, trueVars, whichMethods, what=c("TPR", "TDR", "FDR", "propallselected"), ylimit=c(0,1), maxm=100){
  # files: paths + names of .RData files containing simulation results
  color=c("black","red","darkgreen","blue","magenta","darkgrey")
  lty=c(rep(1,6),rep(2,6),rep(3,6))
  what <- match.arg(what)
  for(l in files){
    sim <- sapply(strsplit(l, "/"), function(x) x[length(x)])
    sim <- sapply(strsplit(sim, "\\."), function(x) x[1])
    print(l)
    load(l)
    if(!missing(whichMethods))
      output <- output[whichMethods,,]
    n.methods <- dim(output)[1]
    N <- dimnames(output)[[1]]
    par(mar=c(4,4,1,1))
    my.perfplot(output[1,,], trueVars, what=what, maxm=maxm, ylimit=ylimit, col=color[1], lty=lty[1], xlim=c(25,maxm))
    for(i in 2:n.methods)
      my.perfplot(output[i,,], trueVars, what=what, maxm=maxm, ylimit=ylimit, col=color[i], lty=lty[i], add=TRUE)
#    legend(par("usr")[2], par("usr")[4], N, col=1:n.methods, lwd=2, lty=lty, xpd=TRUE)
  }
}

#' similarities btw. all methods, corrplot
#'
#' @author Manuela Hummel
#' @param simoutput array w. output of simulation study from Axel
#' @param trueVars vector of indices of true variables
#' @param sim.method measure of similarity btw. methods; either Jaccard index for detected features or correlation of feature ranks
#' @param whichMethods names of methods that shall be included in the plot; if missing, all methods in simoutput are shown
#' @param order if/how methods shall be ordered in corrplot, see ?corrplot, here default is "hclust"
#' @param hclust.method agglomeration method in case methods are ordered by 'hclust', see ?corrplot
#' @param m scalar or vector of model sizes
#' @param ... further parameters passed to corrplot()
#' @details wrapper for \code{corrplot} function R package corrplot
#' @export

# similarities btw. all methods, corrplot
methods.simplot <- function(simoutput, trueVars, whichMethods, sim.method=c("Jaccard","corr"), m=50,
                            order=c("hclust", "original", "AOE", "FPC", "alphabet"),
                            hclust.method=c("complete", "ward", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"), ...){
# simoutput: array w. output of simulation study from Axel
# trueVars: vector of indices of true variables
# whichMethods: names of methods that shall be included in the plot; if missing, all methods in simoutput are shown
# sim.method: measure of similarity btw. methods; either Jaccard index for detected features or correlation of feature ranks
# m: model size for which to determine selected features; ignored if sim.method="corr"
# order: if/how methods shall be ordered in corrplot, see ?corrplot, here default is "hclust"
# hclust.method: agglomeration method in case methods are ordered by 'hclust', see ?corrplot
# ...: further parameters passed to corrplot()



  sim.method <- match.arg(sim.method)
  order <- match.arg(order)
  hclust.method <- match.arg(hclust.method)

  if(!missing(whichMethods))
    simoutput <- simoutput[whichMethods,,]

  n.methods <- dim(simoutput)[1]
  N <- dimnames(simoutput)[[1]]
  D <- matrix(NA, nrow=n.methods, ncol=n.methods)
  dimnames(D) <- list(N, N)
  diag(D) <- 1

  # similarity matrix
  for(i in 1:n.methods){
    for(j in 1:n.methods)
      if(j < i)
        D[i,j] <- method.similarity(simoutput[i,,], simoutput[j,,], trueVars, sim.method, m)
  }
  D[upper.tri(D)] <- t(D)[upper.tri(D)]

  # hierarchical clustering
  #hc <- hclust(as.dist(1-D), method=hclust.method)
  #plot(hc, xlab="", ylab="", sub="", ...)

  # corrplot
  lim <- switch(sim.method == "Jaccard", c(0, 1), NULL)  # ifelse instead of switch would not return NULL and give an error !!
  corrplot::corrplot(D, order=order, hclust.method=hclust.method, cl.lim=lim, ...)
}


#' similarities btw. all methods, corrplot
#'
#' @author Manuela Hummel
#' @param files paths + names of .RData files containing simulation results
#' @param trueVars vector of indices of true variables
#' @param sim.method measure of similarity btw. methods; either Jaccard index for detected features or correlation of feature ranks
#' @param whichMethods names of methods that shall be included in the plot; if missing, all methods in simoutput are shown
#' @param order if/how methods shall be ordered in corrplot, see ?corrplot, here default is "hclust"
#' @param hclust.method agglomeration method in case methods are ordered by 'hclust', see ?corrplot
#' @param m scalar or vector of model sizes
#' @param ... further parameters passed to corrplot()
#' @details wrapper for \code{corrplot} function R package corrplot
#' @export

# plot for all simulation settings
plotall.methodsim <- function(files, trueVars, whichMethods, sim.method=c("Jaccard","corr"), m=50,
                              order=c("hclust", "original", "AOE", "FPC", "alphabet"),
                              hclust.method=c("complete", "ward", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"), ...){
# files: paths + names of .RData files containing simulation results
# trueVars: vector of indices of true variables
# whichMethods: names of methods that shall be included in the plot; if missing, all methods in simoutput are shown
# sim.method: measure of similarity btw. methods; either Jaccard index for detected features or correlation of feature ranks
# m: model size for which to determine selected features; ignored if sim.method="corr"
# order: if/how methods shall be ordered in corrplot, see ?corrplot, here default is "hclust"
# hclust.method: agglomeration method in case methods are ordered by 'hclust', see ?corrplot
# ...: further parameters passed to corrplot()

  sim.method <- match.arg(sim.method)
  order <- match.arg(order)
  hclust.method <- match.arg(hclust.method)

  for(l in files){
    sim <- sapply(strsplit(l, "/"), function(x) x[length(x)])
    sim <- sapply(strsplit(sim, "\\."), function(x) x[1])
    print(l)
    load(l)
    methods.simplot(output, trueVars, whichMethods, sim.method, m, order, hclust.method, main=sim, mar=c(0,0,2,0))
  }
}


