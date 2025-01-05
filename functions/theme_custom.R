
 # Description:  Custom ggplot2 style

 # Arguments:
 #   base_size   - size of the fonts
 #   base_family - main font family
 
 function (base_size = 12, base_family = "")
 {
 theme_bw(base_size = base_size, base_family = base_family) %+replace% 
 theme(
 panel.grid.major = element_line(colour = "grey60", 
 size = 0.5, linetype = "longdash"),
 panel.grid.minor = element_blank(),
 # panel.background = element_rect(fill="grey90",colour = "black", size = 1, linetype="solid"),
 axis.title   = element_text(colour = "black"),
 axis.title.x = element_text(hjust = 0),
 plot.title   = element_text(size = rel(1.1), hjust = 0.5, 
 vjust = 1, margin = margin(b = base_size/2)),
 legend.position = "bottom", # c(0.15, 0.8),
 legend.title = element_blank(),
 aspect.ratio = 0.5
 )
 }
