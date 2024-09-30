library(officer)
library(rvg)
library(ggpubr)

to_pptx <- function(draft, saveFile = "draft.pptx"){
  p <- draft
  p_dml <- rvg::dml(ggobj = p)
  # initialize PowerPoint slide ----
  officer::read_pptx() %>%
    # add slide ----
  officer::add_slide() %>%
    # specify object and location of object ----
  officer::ph_with(p_dml, ph_location()) %>%
    # export slide -----
  base::print(
    target = here::here(
      saveFile
    )
  )
}
