#' Fits a gamma mixture model to mIF imaging data with multiple slides at one time
#'
#' Takes a cfGMM model and forces the right-most component
#' probabilities to be monotonic increasing with respect
#' to the input values. Apply this model to multiple slides with multiple channels.
#'
#' @param expressionMarkers An ncell by nmarker matrix of expression values. Assumed to be coarsely normalized and transformed.
#' @param slide Single vector of length n. Identifies slide number for each cell.
#' @param subBatch If there are multiple subBatch on a slide, subBatch can be used to return probability estimates independently for each region.
#' @param boundaryMarkers A nmarker list of 4x4 matrices giving the boundaries for the modes of the unexpressed and expressed cell distributions.
#' @param qboundaryMarkers A nmarker list of 4x4 matrices giving the qauntile boundaries for the modes for the unexpressed and expressed cell distributions.
#' @param n.cores Number of cores. Default to one less than the number of cores available on the machine.
#' @param ... Arguments passed to GammaGateR.
#' @return Analogous to GammaGateR, groupGammaGateR simpky returns lists in the name of slides, each list being corresponding GammaGateR object.
#' @importFrom cfGMM cfGMM
#' @importFrom stats quantile
#' @importFrom reactable reactable
#' @importFrom parallel mclapply detectCores
#' @export
#' @details Fits cfGMM models to each marker channel in a matrix of marker channels for multiple slides
groupGammaGateR <- function(expressionMarkers, slide, boundaryMarkers=NULL, qboundaryMarkers=NULL, subBatch=NULL, n.cores=NULL, ...){
  if(is.null(n.cores)){n.cores = detectCores()-1}
  if(any(is.na(slide))){expressionMarkers = expressionMarkers[!is.na(slide)]; slide = slide[!is.na(slide)]}
  expressionMarkers = split(expressionMarkers, slide)
  constrCfGMMbunch = mclapply(expressionMarkers, function(x) GammaGateR(x, boundaryMarkers, qboundaryMarkers, subBatch, ...),   mc.cores = n.cores)
  class(constrCfGMMbunch) = c('groupGammaGateR', class(constrCfGMMbunch))
  return(constrCfGMMbunch)
}

#' Inner: Plot function of fitted model in groupGammaGateR object
#'
#' Takes a slide and marker input, and returns the corresponding histogram of expression value,
#' along with fitted density curve.
#'
#' @param x A groupGammaGateR object
#' @param marker Select which markers to plot. Can be a vector of character or integer. If not specified, will print all.
#' @param slide Select which slides to plot. Can be a vector of character or integer. If not specified, will print all.
#' @param component Integer specifying which component to plot, 1 is unexpressed nonzero cells, 2 is expressed cells.
#' @param diagnostic logical indicating whether to create the diagnostic plot. Default value is TRUE.
#' @param interactive logical indicating whether diagnostic plot should be interactive.
#' @param histogram logical indicating whether to create the slide histograms.
#' @param title Title for the plot. Default is the marker name.
#' @param boundary Boundary (vertial dashed line) to be plotted on the histogram.
#' @param color color for points.
#' @param print logical whether to display the plot. Default value TRUE.
#' @param tabl logical. whether to include parameter table in the plot. Default value FALSE.
#' @param ... Arguments passed to XX
#' @importFrom ggplot2 aes ggplot ggtitle stat_function geom_vline unit annotation_custom geom_point xlab ylab
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom stats na.omit
#' @importFrom utils capture.output
#' @importFrom plotly plot_ly layout
#' @details Various diagnostic and QC plots for groupGammaGateR fits.
plotGroupGammaGateR <- function(x, marker=1, slide=1, component=2, diagnostic=TRUE, interactive=FALSE,
                             histogram=FALSE, title=NULL, boundary = NULL,color='grey', print=TRUE, tabl=FALSE,  ...){
  markerind = marker
  if(is.numeric(marker)){marker=colnames(x[[1]][["expressionX"]])[marker] }
  if(diagnostic){
    plot.df = as.data.frame(do.call(rbind, lapply(x, function(y){
      do.call(rbind, lapply(y$params[[marker]], function(z) z[,paste0('comp',component)]) )
    } ) ))
    plot.df$mode <- (plot.df$alpha-1)*plot.df$beta
    if(interactive){
      plot.df$slide_id <- names(x)
      p1 <- plot_ly( plot.df, x = ~mode, y = ~lambda,
                     marker = list(size = 10,color = 'rgba(255, 182, 193, .9)',
                                   line = list(color = 'rgba(152, 0, 0, .8)',width = 2)),
                     text = ~slide_id, type = "scatter", mode = 'markers')
      p1 <- layout(p=p1, title = marker,
                   yaxis = list(zeroline = FALSE),
                   xaxis = list(zeroline = FALSE))
    }else{
      if(is.null(title)){title <- "Diagnostic Plot"}
      p1 <- ggplot(plot.df)+  geom_point(aes(x=mode,y=lambda), alpha=0.4, color=color, size=3) +
        xlab("Mode of Expressed Component")+ylab("Lambda of Expressed")+
        theme_ipsum(axis_title_size = 12) + ggtitle(title)
    }

    if(print){
      print(p1)
    } else {
      return(p1)
    }
  }
  if(histogram){
    p <- plot(x[[slide]], marker = markerind, title = title, boundary=boundary, tabl=tabl, print=FALSE)
    if(print){
      print(p)
    } else {
      return(p)
    }
  }
}

#' Plot function of fitted model in groupGammaGateR object
#'
#' Takes a slide and marker input, and returns the corresponding histogram of expression value,
#' along with fitted density curve.
#'
#' @param x A groupGammaGateR object
#' @param marker Select which markers to plot. Can be a vector of character or integer. If not specified, will print all.
#' @param slide Select which slides to plot. Can be a vector of character or integer. If not specified, will print all.
#' @param component Integer specifying which component to plot, 1 is unexpressed nonzero cells, 2 is expressed cells.
#' @param diagnostic logical indicating whether to create the diagnostic plot. Default value is TRUE.
#' @param interactive logical indicating whether diagnostic plot should be interactive.
#' @param histogram logical indicating whether to create the slide histograms.
#' @param title Title for the plot. Default is the marker name.
#' @param boundary Boundary (vertial dashed line) to be plotted on the histogram.
#' @param color color for points.
#' @param print logical whether to display the plot. Default value TRUE.
#' @param tabl logical. whether to include parameter table in the plot. Default value FALSE.
#' @param ... Arguments passed to XX
#' @importFrom ggplot2 aes ggplot ggtitle stat_function geom_vline unit annotation_custom geom_point xlab ylab
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom stats na.omit
#' @importFrom utils capture.output
#' @importFrom plotly plot_ly layout
#' @export
#' @details Various diagnostic and QC plots for groupGammaGateR fits.
plot.groupGammaGateR <- function(x, marker=NULL, slide=NULL, component=2, diagnostic=TRUE, interactive=FALSE,
                              histogram=FALSE, title=NULL, boundary = NULL,color='grey', print=TRUE, tabl=FALSE, ...){
  if(is.numeric(marker)){marker=colnames(x[[1]][["expressionX"]])[marker] }
  if(diagnostic){
    if(is.null(marker)){markern=names(x[[1]][[1]])}else{markern=marker}
    lapply(as.list(markern), function(mkr)plotGroupGammaGateR(x, marker=mkr, diagnostic=TRUE, interactive = interactive, histogram=FALSE,print=TRUE, title=mkr))
    if(!print){
      return(lapply(as.list(markern), function(mkr)plotGroupGammaGateR(x, marker=mkr, diagnostic=TRUE, interactive = interactive, histogram=FALSE,print=FALSE, title=mkr)))
    }
  }
  if(histogram){
    if(is.null(marker) | is.null(slide)){
      if(is.null(marker))marker=names(x[[1]][[1]])
      if(is.null(slide))slide=names(x)
      ms=expand.grid(marker, slide)
      marker=ms[,1];slide=ms[,2]
      if(is.null(title))title=paste0(slide, "\n", marker)}
    ps <- list()
    for(i in 1:length(marker)){
      marker.i <- as.character(marker[i])
      slide.i <- as.character(slide[i])
      title.i <- title[i]
      if(is.null(boundary)){
        ps[[i]] <- plotGroupGammaGateR(x, marker=marker.i, slide=slide.i, diagnostic = FALSE, histogram = TRUE, title = title.i, boundary = NULL, tabl=tabl, print = FALSE)
      } else {
        ps[[i]] <- plotGroupGammaGateR(x, marker=marker.i, slide=slide.i, diagnostic = FALSE, histogram = TRUE, title = title.i, boundary = boundary[i], tabl=tabl, print = FALSE)
      }
    }
    if(print){
      lapply(ps, print)
    } else {return(ps)}
  }
}

#' Binds groupGammaGateR
#'
#'
#' Extends the c() for groupGammaGateR object
#'
#' @param ... groupGammaGateR objects.
#' @return Analogous to GammaGateR, groupGammaGateR simpky returns lists in the name of slides, each list being corresponding GammaGateR object.
#' @importFrom cfGMM cfGMM
#' @importFrom stats quantile
#' @importFrom reactable reactable
#' @importFrom parallel mclapply detectCores
#' @export
#' @details Fits cfGMM models to each marker channel in a matrix of marker channels for multiple slides
collateGroupGammaGateR <- function(...){
  input.list = list(...)
  out.list <- list()
  for(i in 1:length(input.list)){
    obj.i <- input.list[[i]]
    for(j in 1:length(obj.i)){
      out.list <- append(out.list, list(obj.i[[j]]))
      names(out.list)[[length(out.list)]] <- names(obj.i)[j]
    }
  }
  class(out.list) = c('groupGammaGateR', 'list')
  return(out.list)
}
