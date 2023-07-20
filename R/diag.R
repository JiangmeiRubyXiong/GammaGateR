#' Checks convergence of GammaGateR/groupGammaGateR object
#'
#' This needs to be run before any plot of fitted result being made.
#'
#' @param GammaGateR.fit GammaGateR or groupGammaGateR object.
#' @param out.all Whether to print out the full convergence result for all markers and slides, or just summary of non-convergence.
#' @param return.result Whether to return the convergence result into an object.
#' @return Returns simple printout or dataframe of convergence result, depending on argument input.
#' @importFrom reactable reactable
#' @export
#' @details convergence check to GammaGateR or batch GammaGateR fits
convCheck <- function(GammaGateR.fit, out.all=FALSE, return.result=FALSE){
  class.fit <- class(GammaGateR.fit)[[1]]
  if(class.fit=="GammaGateR"){
    conv.vec <- sapply(GammaGateR.fit$fit, function (x) x$convergence)
    if(out.all){
      print(conv.vec)
    } else {
      if(!all(conv.vec)){
        print(conv.vec[which(!conv.vec)])
      } else {print("All converged")}
    }
  } else if (class.fit=="groupGammaGateR"){
    conv.tab <- sapply(GammaGateR.fit, function(x) sapply(x$fit, function(y){ y$convergence}))
    if(out.all){
      reactable(conv.tab)
    } else{
      if(any(is.na(conv.tab))){
        na.ind <- which(is.na(conv.tab))
        na.names <- expand.grid(dimnames(conv.tab))[na.ind,]
        colnames(na.names) <- c("marker", "slide")
        print("Convergence NA:")
        print(na.names, row.names = FALSE)
        conv.tab[na.ind] <- TRUE
        conv.msg <- "All converged for other slides."} else {conv.msg <- "All converged."}
      if(!all(conv.tab)){
        uncov.ind <- which(!conv.tab)
        uncov.names <- expand.grid(dimnames(conv.tab))[uncov.ind,]
        colnames(uncov.names) <- c("marker", "slide")
        print("Did not converge:")
        print(uncov.names, row.names = FALSE)
      } else {print(conv.msg)}
    }
  } else {print("Wrong input type!!")}
  if(return.result){return(conv.tab)}
}

#' Diagnostic plot - contrast
#'
#' This can be seen as an overlay of the scatter plot diagnostic plot. While side-by-side comparison of the plots
#' might have scale problem, "overlay" the scatter plot resolves this issue. The corresponding slides are also connected
#' with dashed line, making the changes easier to spot.
#'
#' @param fit1 A groupGammaGateR object.
#' @param fit2 Another groupGammaGateR object, possibly that is a refit of fit1 with updated constraints.
#' @param marker Marker to be plotted, defaults to the first one.
#' @param component Integer specifying which component to plot, 1 is unexpressed nonzero cells, 2 is expressed cells.
#' @param title Title of plot.
#' @param fit.names Names of fit objects. Default to c("fit1", "fit2")
#' @param boundaries Add vertical boundary lines to plot.
#' @param ... Arguments passed to plot.GammaGateR
#' @importFrom ggplot2 aes ggplot ggtitle geom_histogram after_stat stat_function geom_vline unit annotation_custom geom_segment scale_color_manual
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom gridExtra tableGrob
#' @import patchwork
#' @export
#' @details Plot a histogram of a GammaGateR model. Plots one model.
#' Takes a GammaGateR object and plots the histogram with the fitted model and parameter values.
diag_contrast <- function(fit1, fit2, marker=1, component=2, title="Contrastive Diagnostic Plot", fit.names=c("fit1", "fit2"),boundaries=NULL,...){
  plot.df.fit1 = as.data.frame(do.call(rbind, lapply(fit1, function(y){
    do.call(rbind, lapply(y$params[[marker]], function(z) z[,paste0('comp',component)]))})))
  plot.df.fit1$mode <- (plot.df.fit1$alpha-1)*plot.df.fit1$beta
  plot.df.fit1$fit <- fit.names[1]
  plot.df.fit2 = as.data.frame(do.call(rbind, lapply(fit2, function(y){
    do.call(rbind, lapply(y$params[[marker]], function(z) z[,paste0('comp',component)]))})))
  plot.df.fit2$mode <- (plot.df.fit2$alpha-1)*plot.df.fit2$beta
  plot.df.fit2$fit <- fit.names[2]
  plot.df <- rbind(plot.df.fit1, plot.df.fit2)
  connect.df <- data.frame(x1=plot.df.fit1[,"mode"],
                           y1=plot.df.fit1[,"lambda"],
                           x2=plot.df.fit2[,"mode"],
                           y2=plot.df.fit2[,"lambda"])

  p <-  ggplot()+
    geom_point(data=plot.df, aes(x=mode, y=lambda,color=fit),alpha=0.3, size=2.5)+
    theme_ipsum()+ scale_color_manual(values=c("red","forestgreen","burlywood4"))+
    geom_segment(data=connect.df, aes(x = x1, y = y1, xend = x2, yend = y2), linetype="dotted")
  if(!is.null(boundaries)){p <- geom_vline(xintercept = boundaries, color="red", linetype="dashed")}
  if(!is.null(title)){p <- p + ggtitle(title)}
  return(p)
}
