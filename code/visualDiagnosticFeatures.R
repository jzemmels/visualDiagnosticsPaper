# Requires:
## a directory called visualDiagnostics_entireScanPairs
## the topMatchData_comb data file
## a 4/7/22 version of cmcR (available on Github)
## a 4/7/22 version of impressions (available on Github)

library(cmcR)
library(impressions)
library(tidyverse)

future:::ClusterRegistry("stop")

future::plan(future::multisession(workers = future::availableCores() - 9))

furrr::future_walk(
  # 1,
  1:(nrow(topMatchData_comb) - 1),
  function(refInd){

    walk(
      # 2,
      (refInd + 1):nrow(topMatchData_comb),
      function(targInd){

        refName <- topMatchData_comb[refInd,]$scanName

        targName <- topMatchData_comb[targInd,]$scanName

        if(file.exists(paste0("visualDiagnostics_entireScanPairs/",refName,"_vs_",targName,".RData"))){

          # print(paste0("visualDiagnostics_entireScanPairs/",refName,"_vs_",targName,".RData"))

          return(NULL)
        }

        refScan <- topMatchData_comb[refInd,]$processedScan[[1]]

        targScan <- topMatchData_comb[targInd,]$processedScan[[1]]

        # browser()

        if(!isTRUE(all.equal(refScan$header.info$incrementX,targScan$header.info$incrementX))){

          if(refScan$header.info$incrementX > targScan$header.info$incrementX){

            targScan <- x3ptools::x3p_interpolate(targScan,resx = refScan$header.info$incrementX)

          }
          else{

            refScan <- x3ptools::x3p_interpolate(refScan,resx = targScan$header.info$incrementX)

          }
        }

        print(paste0(refName,"_vs_",targName))

        compData1 <-
          map_dfr(seq(-30,30,by = 3),#seq(-180,180,by = 3),
                  function(theta){

                    cmcR::comparison_allTogether(refScan,targScan,theta = theta,
                                                 numCells = c(1,1),
                                                 maxMissingProp = .99,
                                                 sideLengthMultiplier = 1,
                                                 returnX3Ps = TRUE)
                    #
                  }) %>%
          filter(fft_ccf == max(fft_ccf)) %>%
          mutate(originalMethod = "CMC",
                 direction = "reference_vs_target",
                 comparisonName = paste0(refName,"_vs_",targName))

        compData2 <-
          map_dfr(seq(-30,30,by = 3),#seq(-180,180,by = 3),
                  function(theta){

                    cmcR::comparison_allTogether(targScan,refScan,theta = theta,
                                                 numCells = c(1,1),
                                                 maxMissingProp = .99,
                                                 sideLengthMultiplier = 1,
                                                 returnX3Ps = TRUE)
                    #
                  }) %>%
          filter(fft_ccf == max(fft_ccf)) %>%
          mutate(originalMethod = "CMC",
                 direction = "target_to_reference",
                 comparisonName = paste0(refName,"_vs_",targName))

        compData <- bind_rows(compData1,compData2)

        visualDiagnosticFeatures <- compData %>%
          select(comparisonName,cellIndex,direction,cellHeightValues,alignedTargetCell) %>%
          pmap_dfr(~ {

            reference <- ..4
            target <- ..5

            reference$surface.matrix <- (reference$surface.matrix*reference$cmcR.info$scaleByVal + reference$cmcR.info$centerByVal)*1e6
            target$surface.matrix <- (target$surface.matrix*target$cmcR.info$scaleByVal + target$cmcR.info$centerByVal)*1e6

            # pixelwise average and absolute difference between the two scans
            refTargAverage <- bind_rows(reference %>%
                                          impressions::x3pToDF() %>%
                                          mutate(value = value),
                                        target %>%
                                          impressions::x3pToDF() %>%
                                          mutate(value = value)) %>%
              group_by(x,y) %>%
              summarise(valueDiff = abs(diff(value)),
                        value = mean(value),
                        .groups = "drop")

            ret <- map_dfr(c(.75,1,1.5,2),
                           function(mult){

                             labeledBlobs <- impressions::labelBlobs(reference,target,
                                                                     filterCutoff = mult*sd(c(c(reference$surface.matrix),
                                                                                              c(target$surface.matrix)),
                                                                                            na.rm = TRUE))

                             labeledBlobs <- labeledBlobs %>%
                               #   mutate(value = factor(value,labels = c("Non-filtered",1:max(labeledBlobs$value,na.rm = TRUE)))) %>%
                               rename(blobLabel = value)

                             nonFilteredObs <- labeledBlobs %>%
                               filter(!is.na(blobLabel) & blobLabel > -1)#"Non-filtered")

                             allObs <- labeledBlobs %>%
                               filter(!is.na(blobLabel))

                             # add to the df created above the pixelwise difference between the reference
                             # and average scans. alpha-blend pixels
                             refDifference <- reference %>%
                               impressions::x3pToDF() %>%
                               left_join(refTargAverage %>%
                                           rename(aveValue = value),
                                         by = c("x","y")) %>%
                               mutate(value = ifelse(is.na(aveValue),NA,value)) %>%
                               mutate(value = ifelse(value > mult*sd(c(c(reference$surface.matrix),
                                                                       c(target$surface.matrix)),
                                                                     na.rm = TRUE),
                                                     value,NA))

                             targDifference <- target %>%
                               impressions::x3pToDF() %>%
                               left_join(
                                 refTargAverage %>%
                                   rename(aveValue = value),
                                 by = c("x","y")) %>%
                               mutate(value = ifelse(is.na(aveValue),NA,value)) %>%
                               mutate(value = ifelse(value > mult*sd(c(c(reference$surface.matrix),
                                                                       c(target$surface.matrix)),
                                                                     na.rm = TRUE),
                                                     value,NA))

                             data.frame(comparisonName = ..1,
                                        cellIndex = ..2,
                                        direction = ..3,
                                        sdMultiplier = mult,
                                        sdValue = sd(c(c(reference$surface.matrix),
                                                       c(target$surface.matrix)),
                                                     na.rm = TRUE)) %>%
                               mutate(differenceCor_nonStandardized = cor(refDifference$value,
                                                                          targDifference$value,
                                                                          use = "pairwise.complete.obs"),
                                      blobSize_ave_nonStandardized = nonFilteredObs %>%
                                        group_by(blobLabel) %>%
                                        tally() %>%
                                        pull(n) %>%
                                        mean(na.rm = TRUE),
                                      blobSize_sd_nonStandardized = nonFilteredObs %>%
                                        group_by(blobLabel) %>%
                                        tally() %>%
                                        pull(n) %>%
                                        sd(na.rm = TRUE))

                           })

            return(ret)

          })

        visualDiagnostics_entireScanPairs <- visualDiagnosticFeatures %>%
          left_join(compData %>%
                      select(cellIndex,direction,x,y,theta,fft_ccf,pairwiseCompCor),
                    by = c("cellIndex","direction"))

        save(visualDiagnostics_entireScanPairs,file = paste0("visualDiagnostics_entireScanPairs/",refName,"_vs_",targName,".RData"))

      })

  })


#Helper function to pad the rows/columns of a matrix
padByDim <- function(reference,
                     target,
                     dimToPad,
                     side = "pre"){

  #This function assumes that dimToPad represents the difference in dimension
  #between the reference scan and target scan. It will pad whichever scan has
  #the smaller dimension (represented by a negative value in dimToPad)

  if(side == "pre"){
    if(dimToPad[1] < 0){
      reference$surface.matrix <- rbind(matrix(NA,
                                               nrow = abs(dimToPad[1]),
                                               ncol = ncol(reference$surface.matrix)),
                                        reference$surface.matrix)
      reference$header.info$sizeY <- ncol(reference$surface.matrix)
      reference$header.info$sizeX <- nrow(reference$surface.matrix)
    }
    if(dimToPad[1] > 0){
      target$surface.matrix <- rbind(matrix(NA,
                                            nrow = abs(dimToPad[1]),
                                            ncol = ncol(target$surface.matrix)),
                                     target$surface.matrix)
      target$header.info$sizeY <- ncol(target$surface.matrix)
      target$header.info$sizeX <- nrow(target$surface.matrix)
    }

    if(dimToPad[2] < 0){
      reference$surface.matrix <- cbind(matrix(NA,
                                               ncol = abs(dimToPad[2]),
                                               nrow = nrow(reference$surface.matrix)),
                                        reference$surface.matrix)
      reference$header.info$sizeY <- ncol(reference$surface.matrix)
      reference$header.info$sizeX <- nrow(reference$surface.matrix)
    }
    if(dimToPad[2] > 0){
      target$surface.matrix <- cbind(matrix(NA,
                                            ncol = abs(dimToPad[2]),
                                            nrow = nrow(target$surface.matrix)),
                                     target$surface.matrix)
      target$header.info$sizeY <- ncol(target$surface.matrix)
      target$header.info$sizeX <- nrow(target$surface.matrix)
    }
  }
  if(side == "post"){
    if(dimToPad[1] < 0){
      reference$surface.matrix <- rbind(reference$surface.matrix,
                                        matrix(NA,
                                               nrow = abs(dimToPad[1]),
                                               ncol = ncol(reference$surface.matrix)))
      reference$header.info$sizeY <- ncol(reference$surface.matrix)
      reference$header.info$sizeX <- nrow(reference$surface.matrix)
    }
    if(dimToPad[1] > 0){
      target$surface.matrix <- rbind(target$surface.matrix,
                                     matrix(NA,
                                            nrow = abs(dimToPad[1]),
                                            ncol = ncol(target$surface.matrix)))
      target$header.info$sizeY <- ncol(target$surface.matrix)
      target$header.info$sizeX <- nrow(target$surface.matrix)
    }

    if(dimToPad[2] < 0){
      reference$surface.matrix <- cbind(reference$surface.matrix,
                                        matrix(NA,
                                               ncol = abs(dimToPad[2]),
                                               nrow = nrow(reference$surface.matrix)))
      reference$header.info$sizeY <- ncol(reference$surface.matrix)
      reference$header.info$sizeX <- nrow(reference$surface.matrix)
    }
    if(dimToPad[2] > 0){
      target$surface.matrix <- cbind(target$surface.matrix,
                                     matrix(NA,
                                            ncol = abs(dimToPad[2]),
                                            nrow = nrow(target$surface.matrix)))
      target$header.info$sizeY <- ncol(target$surface.matrix)
      target$header.info$sizeX <- nrow(target$surface.matrix)
    }
  }

  return(list(reference,target))
}

#Pads two scans such that (1) their central pixels coincide with their
#respective (estimated) center of the firing pin hole and (2) the matrices are
#the same size. If preProcess_rotateScan is used before this function, we can
#assume that the two scans are translationally & rotationally aligned. The
#differences argument then indicates whether only the intersection of the two
#surfaces should be considered.
alignScans <- function(reference,
                       target,
                       reference_center,
                       target_center){

  if(any((reference_center - target_center) != 0)){

    #Pad the matrices so that the centers of their respective firing pin impression
    #holes occupy the same index (shoves both matrices into the top-left corner)
    centerPadded <- padByDim(reference,
                             target,
                             dimToPad = reference_center - target_center,
                             side = "pre")
  }
  else{

    centerPadded <- list(reference,target)

  }

  #Then perform extra padding to make both matrices the same dimension (shoves
  #both matrices into the top-left corner)
  sameSizePadded <- padByDim(centerPadded[[1]],
                             centerPadded[[2]],
                             dimToPad = dim(centerPadded[[1]]$surface.matrix) - dim(centerPadded[[2]]$surface.matrix),
                             side = "post")

  reference_padded <- sameSizePadded[[1]]
  target_padded <- sameSizePadded[[2]]


  return(list(reference_padded,
              target_padded))
}

rotateSurfaceMatrix <- function(surfaceMat,
                                theta = 0,
                                interpolation = 0){
  surfaceMatFake <- (surfaceMat*10^5) + 1 #scale and shift all non-NA pixels up 1 (meter)
  # imFakeRotated <- :bilinearInterpolation(imFake,theta)
  surfaceMatFakeRotated <- surfaceMatFake %>%
    imager::as.cimg() %>%
    imager::imrotate(angle = theta,
                     interpolation = interpolation, #linear interpolation,
                     cx = floor(nrow(.)/2), #imager treats the rows as the "x" axis of an image
                     cy = floor(ncol(.)/2),
                     boundary = 0) %>% #pad boundary with 0s (dirichlet condition)
    as.matrix()

  surfaceMatFakeRotated[surfaceMatFakeRotated == 0] <- NA
  #shift all of the legitimate pixels back down by 1:
  surfaceMatRotated <- (surfaceMatFakeRotated - 1)/(10^5)

  return(surfaceMatRotated)
}

scanPair_manualRegistration <- function(reference,target,x,y,theta){

  target$surface.matrix <- rotateSurfaceMatrix(target$surface.matrix,
                                               theta = theta)

  translatedScans <- alignScans(reference = reference,
                                target = target,
                                reference_center = dim(reference$surface.matrix)/2,
                                target_center = dim(target$surface.matrix)/2 - c(y,x))

  return(translatedScans)
}

future:::ClusterRegistry("stop")

future::plan(future::multisession(workers = future::availableCores() - 9))

furrr::future_walk(
  list.files("visualDiagnostics_entireScanPairs/",full.names = TRUE),
  #"visualDiagnostics_entireScanPairs/K002eG1_vs_K002eG3.RData",
  function(fileName){

    load(fileName)

    # print(fileName)

    scanNames <- str_split(str_remove(str_remove(fileName,"visualDiagnostics_entireScanPairs/"),".RData"),"_vs_")[[1]]

    refName <- topMatchData_comb %>%
      filter(scanName == scanNames[1]) %>%
      pull(scanName)

    targName <- topMatchData_comb %>%
      filter(scanName == scanNames[2]) %>%
      pull(scanName)

    reference <- topMatchData_comb %>%
      filter(scanName == scanNames[1]) %>%
      pull(processedScan) %>%
      .[[1]]

    target <- topMatchData_comb %>%
      filter(scanName == scanNames[2]) %>%
      pull(processedScan) %>%
      .[[1]]

    if(!isTRUE(all.equal(reference$header.info$incrementX,target$header.info$incrementX))){

      if(reference$header.info$incrementX > target$header.info$incrementX){

        target <- x3ptools::x3p_interpolate(target,resx = reference$header.info$incrementX)

      }
      else{

        reference <- x3ptools::x3p_interpolate(reference,resx = target$header.info$incrementX)

      }
    }

    estimRegistration_refToTarget <- visualDiagnostics_entireScanPairs %>%
      filter(direction == "reference_vs_target") %>%
      select(x,y,theta) %>%
      distinct()

    estimRegistration_targetToRef <- visualDiagnostics_entireScanPairs %>%
      filter(direction == "target_to_reference") %>%
      select(x,y,theta) %>%
      distinct()

    registeredScans_refToTarget <- scanPair_manualRegistration(reference = reference,
                                                               target = target,
                                                               x = estimRegistration_refToTarget$x,
                                                               y = estimRegistration_refToTarget$y,
                                                               theta = estimRegistration_refToTarget$theta)

    registeredScans_targetToRef <- scanPair_manualRegistration(reference = target,
                                                               target = reference,
                                                               x = estimRegistration_targetToRef$x,
                                                               y = estimRegistration_targetToRef$y,
                                                               theta = estimRegistration_targetToRef$theta)

    compData1 <-
      map_dfr(-2:2,
              function(thet){

                cmcR::comparison_allTogether(reference = registeredScans_refToTarget[[1]],
                                             target = registeredScans_refToTarget[[2]],
                                             theta = thet,
                                             numCells = c(8,8),
                                             maxMissingProp = .99,
                                             sideLengthMultiplier = 1.1,
                                             returnX3Ps = TRUE)

              }) %>%
      group_by(cellIndex) %>%
      filter(fft_ccf == max(fft_ccf)) %>%
      ungroup() %>%
      filter(fft_ccf == max(fft_ccf)) %>%
      mutate(originalMethod = "CMC",
             direction = "reference_vs_target",
             comparisonName = paste0(refName,"_vs_",targName))

    compData2 <-
      map_dfr(-2:2,
              function(thet){

                cmcR::comparison_allTogether(reference = registeredScans_targetToRef[[1]],
                                             target = registeredScans_targetToRef[[2]],
                                             theta = thet,
                                             numCells = c(8,8),
                                             maxMissingProp = .99,
                                             sideLengthMultiplier = 1.1,
                                             returnX3Ps = TRUE)

              }) %>%
      group_by(cellIndex) %>%
      filter(fft_ccf == max(fft_ccf)) %>%
      ungroup() %>%
      mutate(originalMethod = "CMC",
             direction = "target_vs_reference",
             comparisonName = paste0(refName,"_vs_",targName))

    compData <- bind_rows(compData1,compData2)

    visualDiagnosticFeatures <- compData %>%
      select(comparisonName,cellIndex,direction,cellHeightValues,alignedTargetCell) %>%
      pmap_dfr(~ {

        reference <- ..4
        target <- ..5

        reference$surface.matrix <- (reference$surface.matrix*reference$cmcR.info$scaleByVal + reference$cmcR.info$centerByVal)*1e6
        target$surface.matrix <- (target$surface.matrix*target$cmcR.info$scaleByVal + target$cmcR.info$centerByVal)*1e6

        # pixelwise average and absolute difference between the two scans
        refTargAverage <- bind_rows(reference %>%
                                      impressions::x3pToDF() %>%
                                      mutate(value = value),
                                    target %>%
                                      impressions::x3pToDF() %>%
                                      mutate(value = value)) %>%
          group_by(x,y) %>%
          summarise(valueDiff = abs(diff(value)),
                    value = mean(value),
                    .groups = "drop")

        ret <- map_dfr(c(.75,1,1.5,2),
                       function(mult){

                         labeledBlobs <- impressions::labelBlobs(reference,target,
                                                                 filterCutoff = mult*sd(c(c(reference$surface.matrix),
                                                                                          c(target$surface.matrix)),
                                                                                        na.rm = TRUE))

                         labeledBlobs <- labeledBlobs %>%
                           #   mutate(value = factor(value,labels = c("Non-filtered",1:max(labeledBlobs$value,na.rm = TRUE)))) %>%
                           rename(blobLabel = value)

                         nonFilteredObs <- labeledBlobs %>%
                           filter(!is.na(blobLabel) & blobLabel > -1)#"Non-filtered")

                         allObs <- labeledBlobs %>%
                           filter(!is.na(blobLabel))

                         # add to the df created above the pixelwise difference between the reference
                         # and average scans. alpha-blend pixels
                         refDifference <- reference %>%
                           impressions::x3pToDF() %>%
                           left_join(refTargAverage %>%
                                       rename(aveValue = value),
                                     by = c("x","y")) %>%
                           mutate(value = ifelse(is.na(aveValue),NA,value)) %>%
                           mutate(value = ifelse(value > mult*sd(c(c(reference$surface.matrix),
                                                                   c(target$surface.matrix)),
                                                                 na.rm = TRUE),
                                                 value,NA))

                         targDifference <- target %>%
                           impressions::x3pToDF() %>%
                           left_join(
                             refTargAverage %>%
                               rename(aveValue = value),
                             by = c("x","y")) %>%
                           mutate(value = ifelse(is.na(aveValue),NA,value)) %>%
                           mutate(value = ifelse(value > mult*sd(c(c(reference$surface.matrix),
                                                                   c(target$surface.matrix)),
                                                                 na.rm = TRUE),
                                                 value,NA))

                         data.frame(comparisonName = ..1,
                                    cellIndex = ..2,
                                    direction = ..3,
                                    sdMultiplier = mult,
                                    sdValue = sd(c(c(reference$surface.matrix),
                                                   c(target$surface.matrix)),
                                                 na.rm = TRUE)) %>%
                           mutate(differenceCor_nonStandardized = cor(refDifference$value,
                                                                      targDifference$value,
                                                                      use = "pairwise.complete.obs"),
                                  blobSize_ave_nonStandardized = nonFilteredObs %>%
                                    group_by(blobLabel) %>%
                                    tally() %>%
                                    pull(n) %>%
                                    mean(na.rm = TRUE),
                                  blobSize_sd_nonStandardized = nonFilteredObs %>%
                                    group_by(blobLabel) %>%
                                    tally() %>%
                                    pull(n) %>%
                                    sd(na.rm = TRUE))

                       })

        return(ret)

      })

    visualDiagnosticsFeatures <- visualDiagnosticFeatures %>%
      left_join(compData %>%
                  select(cellIndex,direction,x,y,theta,fft_ccf,pairwiseCompCor),
                by = c("cellIndex","direction"))

    save(visualDiagnosticsFeatures,file = paste0("topMatch_preRegisteredCellFeatures/",refName,"_vs_",targName,".RData"))

  })