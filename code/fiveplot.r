library(patchwork)

# referenceScan <- "K013sA1"
# targetScan <- "K013sA2"
#
# fiveplot(referenceScan, targetScan, "2, 1")
# fiveplot(referenceScan, targetScan, "2, 7")
# fiveplot(referenceScan, targetScan, "2, 6")


fiveplot <- function(comparisonResults, referenceScan, targetScan, cell) {

  referenceCell <- comparisonResults %>%
    filter(comparisonName == sprintf("%s vs. %s", referenceScan, targetScan) &
             cellIndex == cell) %>%
    pull(cellHeightValues) %>%
    .[[1]]

  targetCell <- comparisonResults %>%
    filter(comparisonName == sprintf("%s vs. %s", referenceScan, targetScan) &
             cellIndex == cell) %>%
    pull(alignedTargetCell) %>%
    .[[1]]

  referenceCell$surface.matrix <- (referenceCell$surface.matrix*referenceCell$cmcR.info$scaleByVal + referenceCell$cmcR.info$centerByVal)*1e6
  targetCell$surface.matrix <- (targetCell$surface.matrix*targetCell$cmcR.info$scaleByVal + targetCell$cmcR.info$centerByVal)*1e6


  patchComparisonPlts <- impressions::x3pComparisonPlot(
    reference = referenceCell,
    target = targetCell,
    cutoffThresh = sd(c(c(referenceCell$surface.matrix),c(targetCell$surface.matrix)),na.rm = TRUE),
    plotNames = c(sprintf("%s Cell %s", referenceScan, cell),
                  sprintf("%s Aligned Cell", targetScan),
                  "Filtered Element-wise Average",
                  sprintf("%s Cell %s\nFiltered Differences", referenceScan, cell),
                  sprintf("%s Aligned Cell\nFiltered Differences", targetScan))
  )

  patchComparisonLegend_match <-
    cowplot::plot_grid(patchComparisonPlts$legend$grobs[[1]])

  combinedValues <-  referenceCell %>%
    impressions::x3pToDF() %>%
    rename(refValue = value) %>%
    left_join(targetCell %>%
                impressions::x3pToDF() %>%
                rename(targValue = value),
              by = c("x","y"))

  blobBoundaries <- comparisonResults %>%
    filter((comparisonName == sprintf("%s vs. %s", referenceScan, targetScan) &
              cellIndex == cell)) %>%
    as.data.frame() %>%
    dplyr::select(comparisonName,cellIndex,cellHeightValues,alignedTargetCell) %>%
    pmap(~ {

      reference <-  ..3

      target <- ..4

      reference$surface.matrix <- (reference$surface.matrix*reference$cmcR.info$scaleByVal + reference$cmcR.info$centerByVal)*1e6
      target$surface.matrix <- (target$surface.matrix*target$cmcR.info$scaleByVal + target$cmcR.info$centerByVal)*1e6

      averageBinarized <- bind_rows(reference %>%
                                      impressions::x3pToDF() %>%
                                      mutate(value = value),
                                    target %>%
                                      impressions::x3pToDF() %>%
                                      mutate(value = value)) %>%
        group_by(x,y) %>%
        summarise(difference = diff(value),
                  absDifference = abs(diff(value)),
                  average = mean(value),
                  .groups = "drop")  %>%
        mutate(comparisonName = ..1,
               cellIndex = ..2)%>%
        mutate(value = ifelse(absDifference > sd(c(c(reference$surface.matrix),c(target$surface.matrix)),na.rm = TRUE),TRUE,FALSE))

      suppressWarnings({

        averageMat <- averageBinarized %>%
          mutate(x = x+1,
                 y=y+1) %>%
          as.data.frame() %>%
          dplyr::select(x,y,value) %>%
          imager::as.cimg() %>%
          as.matrix()

      })

      averageMat[is.na(averageMat)] <- 0

      # we pad the matrix so that the contours one the edge blobs are properly
      # identified. the padding is removed in the last lines of the creation of
      # the outline object below
      averageMat  <- averageMat %>%
        imager::as.cimg() %>%
        imager::pad(nPix = 10,axes = "xy",val = 0)

      labels <- imager::label(averageMat)

      bounds <- map(unique(labels[labels > 0]),
                    function(lab){

                      imager::boundary(labels == lab)

                    })

      return(list(bounds,labels))

    })

  # combine all labeled blobs into one image
  boundaryPx <- Reduce("+",blobBoundaries[[1]][[1]] %>%
                         map(as.matrix)) %>%
    imager::as.cimg()

  # the mask used to dilate the blobs will grow them towards the bottom-right of
  # the matrix
  dilatedPx <- imager::dilate_rect(boundaryPx,sx = 2,sy = 2)
  dilatedPx_labels <- imager::dilate_rect(blobBoundaries[[1]][[2]],sx = 2,sy = 2)

  # flip the image and re-apply the dilation to grow the borders to the other
  # corners. flip back after dilation
  dilatedPx_mirrorx <- imager::mirror(imager::dilate_rect(imager::mirror(boundaryPx,axis="x"),sx = 2,sy = 2),axis="x")
  dilatedPx_mirrorx_labels <- imager::mirror(imager::dilate_rect(imager::mirror(blobBoundaries[[1]][[2]],axis="x"),sx = 2,sy = 2),axis="x")

  dilatedPx_mirrory <- imager::mirror(imager::dilate_rect(imager::mirror(boundaryPx,axis="y"),sx = 2,sy = 2),"y")
  dilatedPx_mirrory_labels <- imager::mirror(imager::dilate_rect(imager::mirror(blobBoundaries[[1]][[2]],axis="y"),sx = 2,sy = 2),"y")

  dilatedPx_mirrorxy <- imager::mirror(imager::dilate_rect(imager::mirror(boundaryPx,axis="xy"),sx = 3,sy = 3),"xy")
  dilatedPx_mirrorxy_labels <- imager::mirror(imager::dilate_rect(imager::mirror(blobBoundaries[[1]][[2]],axis="xy"),sx = 3,sy = 3),"xy")

  # combine all of the dilated images together into one image
  dilatedPx_comb <- dilatedPx + dilatedPx_mirrorx + dilatedPx_mirrory + dilatedPx_mirrorxy

  # we just want a binary labeling
  dilatedPx_comb[dilatedPx_comb > 0] <- 1

  # the dilated boundaries will have also grown into the blobs, so we take those
  # pixels out
  dilatedPx_comb[blobBoundaries[[1]][[2]] > 0] <- 0

  # from: https://stackoverflow.com/questions/34756755/plot-outline-around-raster-cells
  outline <- dilatedPx_comb %>%
    as.data.frame() %>%
    filter(value > 0) %>%
    mutate(x = x-1,
           y = y-1) %>%
    raster::rasterFromXYZ() %>%
    raster::rasterToPolygons(dissolve = TRUE) %>%
    fortify() %>%
    #the boundaries around the filtered blobs all share a common value in the
    #"hole" column of TRUE
    filter(hole) %>%
    # remove padding used previously
    mutate(lat = lat-5,
           long = long-5)


  `-.gg` <- function(plot, layer) {
    if (missing(layer)) {
      stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
    }
    if (!is.ggplot(plot)) {
      stop('Need a plot on the left side')
    }
    plot$layers = c(layer, plot$layers)
    plot
  }

  topLeft <- patchComparisonPlts[[1]] +
    cowplot::theme_nothing() +
    labs(subtitle = sprintf("%s Cell %s", referenceScan, cell)) +
    theme(plot.margin = margin(0,0,5,0),
          plot.subtitle = element_text(hjust = .5,size = 8,vjust = -1)) +
    geom_raster(data = combinedValues %>%
                  filter(is.na(refValue) & !is.na(targValue)),
                fill = "gray40")

  bottomLeft <-patchComparisonPlts[[2]] +
    cowplot::theme_nothing() +
    labs(subtitle = paste0(targetScan," Aligned Cell\nat ",unique(comparisonResults$theta),"°")) +
    theme(plot.margin = margin(-20,-100,30,-100),
          plot.subtitle = element_text(hjust = .5,vjust = -78,size = 8)) +
    geom_raster(data = combinedValues %>%
                  filter(!is.na(refValue) & is.na(targValue)),
                fill = "gray40")

  middle <- patchComparisonPlts[[3]] +
    cowplot::theme_nothing() +
    labs(subtitle = "Filtered Element-wise Average\nAbs. Differences at Most 1") +
    theme(plot.margin = margin(0,25,0,25),
          plot.subtitle = element_text(hjust = .5,size = 8,vjust = -5)) -
    geom_raster(fill = "gray80") +
    geom_path(data = outline,  color = "grey40",
              aes(x=long,y=lat,group=group),
              colour = "gray40",
              inherit.aes = FALSE,
              size = .2)

  topRight <- patchComparisonPlts[[4]] +
    cowplot::theme_nothing() +
    labs(subtitle = sprintf("Filtered %s Cell %s\nAbs. Differences Greater Than 1",referenceScan, cell)) +
    theme(plot.margin = margin(0,0,5,0),
          plot.subtitle = element_text(hjust = .5,size = 8)) -
    geom_raster(fill = "gray80") +
    geom_path(data = outline,  color = "grey40",
              aes(x=long,y=lat,group=group),
              colour = "gray40",
              inherit.aes = FALSE,
              size = .1)

  bottomRight <- patchComparisonPlts[[5]] +
    cowplot::theme_nothing() +
    labs(subtitle = sprintf("Filtered %s Aligned Cell\nAbs. Differences Greater Than 1",targetScan)) +
    theme(plot.margin = margin(-20,-100,30,-100),
          plot.subtitle = element_text(hjust = .5,vjust = -78,size = 8)) -
    geom_raster(fill = "gray80") +
    geom_path(data = outline, color = "grey40",
              aes(x=long,y=lat,group=group),
              colour = "gray40",
              inherit.aes = FALSE,
              size = .1)

  design <- "ACCD\nBCCE"

  patchwork::wrap_plots(topLeft,bottomLeft,middle,topRight,bottomRight,design = design) +
    inset_element(patchComparisonLegend_match,left = -2.15,bottom = 0,right = -2.15,top = 0,on_top = FALSE,align_to = 'full')
}

