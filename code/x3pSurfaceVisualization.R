if(!file.exists("figures/topMatchAnnotatedreference.png")){

  data.frame(x=1:2,y=1:2,value=c(1,1)) %>%
    x3ptools::df_to_x3p() %>%
    x3ptools::x3p_image()

  reference <- x3ptools::read_x3p("data/K013sA1.x3p") %>%
    # x3ptools::x3p_sample(m=4) %>%
    x3ptools::x3p_rotate() %>%
    x3ptools::x3p_flip_y() %>%
    x3ptools::x3p_rotate(angle = 270)

  reference$mask[reference$mask == "#000000FF"] <- "#cd7f32"
    reference$mask[reference$mask == "#FFFFFFFF"] <- "#FF0000"

    reference %>%
      x3ptools::x3p_image(size = c(750,750),zoom=.75)
    x3ptools::x3p_snapshot(file = "figures/topMatchAnnotatedreference.png")

    knitr::plot_crop("figures/topMatchAnnotatedreference.png")

}

if(!file.exists("figures/topMatchAnnotatedTarget.png")){

  data.frame(x=1:2,y=1:2,value=c(1,1)) %>%
    x3ptools::df_to_x3p() %>%
    x3ptools::x3p_image()

  target <- x3ptools::read_x3p("data/K013sA2-2.x3p") %>%
    # x3ptools::x3p_sample(m=4) %>%
    x3ptools::x3p_rotate() %>%
    x3ptools::x3p_flip_y() %>%
    x3ptools::x3p_rotate(angle = 270)

  target$mask[target$mask == "#000000FF"] <- "#cd7f32"
    target$mask[target$mask == "#FFFFFFFF"] <- "#FF0000"

    target %>%
      x3ptools::x3p_image(size = c(750,750),zoom=.75)
    x3ptools::x3p_snapshot(file = "figures/topMatchAnnotatedTarget.png")

    knitr::plot_crop("figures/topMatchAnnotatedTarget.png")


}

if(!file.exists("figures/referenceTargetSideBySide.png")){

  x3pPlts <- cmcR::x3pListPlot(list("K013sA1" = reference_processed,
                                    "K013sA2" = target_processed)
                               ,type = "list",na.value = "gray65") %>%
    map(~ {

      . + theme(legend.position = "none",plot.title = element_text(hjust = .5))

    })


  x3pCombinedPlt <- (x3pPlts[[1]]+
                       theme(plot.title = element_text(vjust = -32,
                                                       hjust = .6,
                                                       size = 18))) /
    (x3pPlts[[2]] +
       theme(plot.title = element_text(vjust = -32,
                                       hjust = .6,
                                       size = 18)))

  ggsave(plot = x3pCombinedPlt,filename = "figures/referenceTargetSideBySide.png",bg ="white")

  knitr::plot_crop("figures/referenceTargetSideBySide.png")

}
