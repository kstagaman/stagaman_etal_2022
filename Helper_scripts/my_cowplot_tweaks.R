# my_cowplot_tweaks.R

my_theme <- theme_update(
  plot.background = element_rect(fill = "white", color = "white"),
  legend.position = "top",
  legend.box = "vertical",
  legend.box.just = "left",
  legend.title = element_text(size = 9),
  legend.text = element_text(size = 8),
  strip.text = element_text(size = 9),
  plot.caption = element_text(hjust = 0, size = 10),
  axis.text = element_text(size = 8),
  axis.title = element_text(size = 10)
)
