#' Median of Minimum Model Size needed to detect all true features (MMMS)
#'
#' @author Manuela Hummel
#' @param scores scores matrix; output of screening methods for simulated data, larger value means better ranking, larger value means better ranking
#' @param trueVars vector of indices of true variables
#' @param fractrue  mmms for selecting \code{fractrue} * \code{trueVars} (e.g. 90\%) of the true variables
#' @details ties are handled using the 'random' approach in \code{rank}
#' @return scalar
#' @export

mmms <- function(scores, trueVars, fractrue=1){
  p <- ncol(scores)

  # variable ranks (larger value = better rank)
  ranks <- apply(scores, 1, function(x) p + 1 - rank(x, ties.method="random", na.last=FALSE))

  # ranks of true variables
  rankstrue <- apply(ranks, 2, function(x) x[trueVars])

  # keep only fraction of fractrue true variables
  keep <- round(fractrue * length(trueVars))
  rankstrue <- apply(rankstrue, 2, function(x) sort(x)[1:keep])

  # mmms = maximal rank within those
  median(apply(rankstrue, 2, max))
}

#' Proportion of simulation runs where single relevant variables are selected for given model size(s)
#'
#' @author Manuela Hummel
#' @param scores scores matrix; output of screening methods for simulated data, larger value means better ranking, larger value means better ranking
#' @param trueVars vector of indices of true variables
#' @param m scalar or vector of model sizes
#' @details ties are handled using the 'random' approach in \code{rank}
#' @return vector with selection proportions for true variables
#' @export

propselected <- function(scores, m, trueVars){

  p <- ncol(scores)
  B <- nrow(scores)

  # variable ranks (larger value = better rank)
  ranks <- apply(scores, 1, function(x) p + 1 - rank(x, ties.method="random", na.last=FALSE))

  # ranks of true variables
  rankstrue <- apply(ranks, 2, function(x) x[trueVars])

  # selection frequency of relevant variables amongst the top m
  tmp <- matrix(rankstrue %in% 1:m, nrow=nrow(rankstrue))
  out <- rowSums(tmp) / B
  names(out) <- trueVars

  return(out)
}

#' Proportion of simulation runs where all relevant variables are selected for given model size(s)
#'
#' @author Manuela Hummel
#' @param scores scores matrix; output of screening methods for simulated data, larger value means better ranking, larger value means better ranking
#' @param trueVars vector of indices of true variables
#' @param m scalar or vector of model sizes
#' @details ties are handled using the 'random' approach in \code{rank}
#' @return scalar or vector with proportion of runs with all true variables being selected
#' @export

propallselected <- function(scores, m, trueVars){
  # m: scalar or vector of model sizes

  p <- ncol(scores)

  # variable ranks (larger value = better rank)
  ranks <- apply(scores, 1, function(x) p + 1 - rank(x, ties.method="random", na.last=FALSE))

  # ranks of true variables
  rankstrue <- apply(ranks, 2, function(x) x[trueVars])

  out <- NULL
  for(i in m){
    # overlap w. top m variables
    tmp <- matrix(rankstrue %in% 1:i, nrow=nrow(rankstrue))

    # fraction where all true variables are within the top m
    out <- c(out, mean(apply(tmp, 2, all)))
  }
  names(out) <- m

  return(out)
}


#' True Positive Rate (TPR) for given model size(s)
#'
#' @author Manuela Hummel
#' @param scores scores matrix; output of screening methods for simulated data, larger value means better ranking, larger value means better ranking
#' @param m scalar or vector of model sizes
#' @param trueVars vector of indices of true variables
#' @details ties are handled using the 'random' approach in \code{rank}. the average number of true variables contained in the top m selected we get by multiplying TPR by the number of true variables
#' @return scalar/vector with TPR for each model size
#' @export

TPR <- function(scores, m, trueVars){
  # m: scalar or vector of model sizes

  p <- ncol(scores)

  # variable ranks (larger value = better rank)
  ranks <- apply(scores, 1, function(x) p + 1 - rank(x, ties.method="random", na.last=FALSE))

  # ranks of true variables
  rankstrue <- apply(ranks, 2, function(x) x[trueVars])

  out <- NULL
  for(i in m){
    # overlap w. top m variables
    tmp <- matrix(rankstrue %in% 1:i, nrow=nrow(rankstrue))

    # average TPR
    out <- c(out, mean(colMeans(tmp)))
  }
  names(out) <- m

  return(out)
}


#' False Discovery Rate (FDR) for given model size(s)
#'
#' @author Manuela Hummel
#' @param scores scores matrix; output of screening methods for simulated data, larger value means better ranking, larger value means better ranking
#' @param m scalar or vector of model sizes
#' @param trueVars vector of indices of true variables
#' @details ties are handled using the 'random' approach in \code{rank}.
#' @return scalar/vector with FDR for each model size
#' @export
#'

FDR <- function(scores, m, trueVars){
  # m: scalar or vector of model sizes

  p <- ncol(scores)

  # variable ranks (larger value = better rank)
  ranks <- apply(scores, 1, function(x) p + 1 - rank(x, ties.method="random", na.last=FALSE))

  # ranks of true variables
  #rankstrue <- apply(ranks, 2, function(x) x[trueVars])

  out <- NULL
  for(i in m){
    # ranks of selected variables
    selected <- apply(ranks, 2, function(x) order(x)[1:i])

    # false discoveries
    tmp <- matrix(!(selected %in% trueVars), nrow=nrow(selected))

    # average FDR
    out <- c(out, mean(colMeans(tmp)))
    #tmp <- apply(rankstrue, 2, function(x) length(setdiff(1:i, x)))
    #out <- c(out, mean(tmp / i))
  }
  names(out) <- m

  return(out)
}


#' method similarity between two methods
#'
#' @author Manuela Hummel
#' @param scores1 scores matrix; output of screening methods for simulated data, larger value means better ranking, rows = simulation runs, columns = variables
#' @param scores2 scores matrix; output of screening methods for simulated data, larger value means better ranking, rows = simulation runs, columns = variables
#' @param trueVars vector of indices of true variables
#' @param sim.method measure of similarity btw. methods; either Jaccard index for detected features or correlation of feature ranks
#' @param m model size for which to determine selected features; ignored if sim.method="corr"
#' @details ties are handled using the 'random' approach in \code{rank}. average Pearson correlation btw. rankings of 2 methods for b'th simulation run. average Jaccard index for detected true variables at given model size
#' @return scalar/vector with FDR for each model size
#' @export


method.similarity <- function(scores1, scores2, trueVars, sim.method=c("Jaccard","corr"), m=50){

  sim.method <- match.arg(sim.method)

  # ranks of true variables
  p <- ncol(scores1)
  ranks1 <- apply(scores1, 1, function(x) p + 1 - rank(x, ties.method="random", na.last=FALSE))
  ranks1 <- apply(ranks1, 2, function(x) x[trueVars])
  ranks2 <- apply(scores2, 1, function(x) p + 1 - rank(x, ties.method="random", na.last=FALSE))
  ranks2 <- apply(ranks2, 2, function(x) x[trueVars])

  B <- nrow(scores1)
  sim <- numeric(B)
  for(b in 1:B){
    if(sim.method == "Jaccard"){
      # selected features at model size m
      sel1 <- names(ranks1[,b])[ranks1[,b] %in% 1:m]
      sel2 <- names(ranks2[,b])[ranks2[,b] %in% 1:m]

      # Jaccard index
      sim[b] <- length(intersect(sel1, sel2)) / length(union(sel1, sel2))
    }

    else if(sim.method == "corr")
      # Pearson correlation of ranks
      sim[b] <- cor(ranks1[,b], ranks2[,b], use="complete.obs")
  }

  return(mean(sim, na.rm=TRUE))
}

