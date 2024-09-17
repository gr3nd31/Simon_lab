library(officer)
library(rvg)
library(ggpubr)

p <- ggplot(data = iris, aes(x = Species,
                             y = Petal.Width))+
  geom_boxplot()+
  geom_jitter(width = 0.3)+
  theme_bw()
p

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
    "demo_one.pptx"
  )
)
