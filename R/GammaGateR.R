#' Fits a gamma mixture model to aid clustering of mIF imaging data
#'
#' Takes a cfGMM model and forces the right-most component
#' probabilities to be monotonic increasing with respect
#' to the input values.
#'
#' @param expressionMarkers An ncell by nmarker matrix of expression values. Assumed to be coarsely normalized and transformed.
#' @param subBatch If there are multiple subBatch on a slide, subBatch can be used to return probability estimates independently for each region.
#' @param boundaryMarkers A nmarker list of 4x4 matrices giving the boundaries for the modes of the unexpressed and expressed cell distributions.
#' @param qboundaryMarkers A nmarker list of 4x4 matrices giving the qauntile boundaries for the modes for the unexpressed and expressed cell distributions.
#' @param ... Arguments passed to cfGMM function
#'
#' @return
#' \item{expressionZ}{The monotone proportional density of second component}
#' \item{expressionX}{The input expression value}
#' \item{expressionW}{The raw proportional density of second component}
#' \item{params}{The monotone proportional density of second component}
#' \item{fit}{cfGMM::cfGMM fitting output}
#' \item{subBatch}{same as the subBatch in the argument}
#'
#' @importFrom cfGMM cfGMM
#' @importFrom stats quantile
#' @export
#' @details Fits cfGMM models to each marker channel in a matrix of marker channels for one slide
GammaGateR <- function(expressionMarkers, boundaryMarkers=NULL, qboundaryMarkers=NULL, subBatch=NULL, ...){
  if(is.null(subBatch)) subBatch = rep(1, nrow(expressionMarkers))
  expressionMarkers = as.data.frame(expressionMarkers)
  # There could be better checks here
  if(!is.null(qboundaryMarkers)){
    newBoundaryMarkers =  mapply(function(expressionMarker, quantileMat, boundaryMat){
      # computes quantiles from qboundaryMarkers
      if(!is.null(quantileMat)){
        matrix(quantile(expressionMarker, quantileMat, na.rm=TRUE), nrow=nrow(quantileMat), ncol=ncol(quantileMat))
      } else {
        matrix(c(-Inf, -Inf, Inf, Inf), nrow=2, ncol=2)
      }
    },
    expressionMarker=expressionMarkers, quantileMat=qboundaryMarkers, SIMPLIFY = FALSE)
    if(!is.null(boundaryMarkers)){
      boundaryMarkers = mapply(function(newBoundaryMat, boundaryMat){
        # chooses most conservative combination of two boundary options
        if(is.null(boundaryMat)){
          boundaryMat = matrix(c(-Inf, -Inf, Inf, Inf), nrow=2, ncol=2)
        }
        cbind(pmax(newBoundaryMat[,1], boundaryMat[,1]), pmin(newBoundaryMat[,2], boundaryMat[,2]) )
      }, newBoundaryMat=newBoundaryMarkers, boundaryMat=boundaryMarkers, SIMPLIFY = FALSE)
    } else {
      boundaryMarkers = newBoundaryMarkers
    }
  }
  # run models
  if(is.null(boundaryMarkers)){ boundaryMarkers = rep(list(boundaryMarkers), ncol(expressionMarkers)) }
  result = mapply(GammaGateRX, x=expressionMarkers, constraints=boundaryMarkers, MoreArgs=list(subBatch=subBatch, ...=...))
  result = list(expressionZ = as.data.frame(do.call(cbind, result[3,])),
                expressionX=as.data.frame(do.call(cbind, result[4,])),
                expressionW=as.data.frame(do.call(cbind, result[5,])),
                params = result[2,], fit=result[1,], subBatch=subBatch)
  class(result) = c('GammaGateR', class(result))
  return(result)
}

#' Fits a gamma mixture model to a vector of expression data
#'
#'
#' @param x A vector of expression values. Assumed to be coarsely normalized and transformed.
#' @param constraints A nmarker list of 4x4 matrices giving the boundaries for the modes of the unexpressed and expressed cell distributions.
#' @param subBatch If there are multiple subBatch on a slide, the vector subBatch can be used to return probability estimates independently for each region.
#' @param nn0 If there are fewer than nn0, then it will not fit the cfGMM. Default is 200.
#' @param ... Arguments passed to cfGMM function
#' @importFrom cfGMM cfGMM
#' @importFrom stats pgamma
#' @details Takes a cfGMM model and forces the right-most component
#' probabilities to be monotonic increasing with respect
#' to the input values.
GammaGateRX = function(x, constraints=NULL, subBatch=NULL, nn0=200, ...){
  n = length(x)
  naInds = is.na(x)
  zeroInds = (x==0)
  if(is.null(subBatch)){
    subBatch = rep(1, length(x))
  }
  y = x
  x = x[!naInds & !zeroInds ]
  u.region <- unique(subBatch)
  # fit without constraints if none are passed
  if((sum(!zeroInds, na.rm=TRUE)>=nn0 & mean(naInds)!=1)|mode(x)!="character"){
    if(is.null(constraints)){
      fit <- cfGMM::cfGMM(x, k=2, ...)
    } else {
      if(all(dim(constraints)!=c(2,2))){
        stop('constraints must be NULL or a 2x2 matrix.')
      } else {
        fit <- cfGMM::cfGMM(x, k=2, constraint = constraints, ...)
      }
    }
    fit <- monoZ(fit)
  } else {
    # cfGMM not fit
    fit = list(posterior=matrix(NA, nrow=sum(!zeroInds, na.rm=TRUE), ncol=2),
               lambda=rep(NA,2),
               gamma.pars=matrix(NA, nrow=2, ncol=2, dimnames=list(c('alpha', 'beta'))),
               convergence=NA)
  }

  post.prob <- data.frame(factor=subBatch[!naInds & !zeroInds], value=fit$posterior)
  prop0 <-  data.frame(factor=subBatch, value=zeroInds)
  prop0.sep <- with(prop0, tapply(value, factor, mean, na.rm=TRUE))

  param.table <- list()
  for(i in 1:length(u.region)){
    tab.region <- rbind((1-prop0.sep[i])*(colMeans(post.prob[post.prob[,1]==u.region[i],2:3])), fit$gamma.pars)
    tab.region <- cbind(comp0=c(prop0.sep[i], NA, NA), tab.region)
    dimnames(tab.region) <- list(c("lambda", "alpha", "beta"), c(paste0("comp", 0:2)))
    param.table[[u.region[i]]] <- tab.region
  }

  post0 <- data.frame(comp1=rep(NA, n), comp2 = rep(NA, n))
  post0[ which(zeroInds),] = 0
  # return probabilities for posteriors
  if(ifelse(is.na(fit$convergence), FALSE, fit$convergence)){
    post0[!naInds & !zeroInds,] <- fit$posterior
  } else {
    param.table = lapply(param.table, function(x){ x[c('alpha', 'beta'), ] = NA; x['lambda', c('comp1', 'comp2')]=NA ; x} )
  }
  post0$x <- y
  # delete large repetitive objects from fit
  fit$x <- fit$posterior <- NULL
  # compute expressionW by subBatch
  w = unsplit(mapply(function(x, pars){pgamma(x, shape=pars[2,3],scale=pars[3,3])}, x=split(post0$x, subBatch), pars=param.table, SIMPLIFY = FALSE), subBatch)
  return(list(fit=fit, param.table=param.table, postZ=post0$comp2, x=post0$x, w=w))
}
