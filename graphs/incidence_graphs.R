library(readxl)
library(dplyr) #data manipulation
library(readxl) #read in excel files 
library(purrr) #pmap to run over functions
library(tibble) #make tibbles
library(ggplot2) #plot
library(ggpubr) 

#function wrapper around graph code originally created by Carlos Sillero Rejon

time_series_graph <- function(data,
                              cohort,
                              cohort_title,
                              var,
                              groups_to_include,
                              steps,
                              max_value,
                              colours){
  
  #unlist the groups_to_include and colours
  if (is.list(colours)) {
    colours <- colours[[1]]
  }
  
  if (is.list(groups_to_include)){
    groups_to_include <- groups_to_include[[1]]
  }
  
  #construct the outcome names
  outcome <- paste0(var, "_", cohort, "_", groups_to_include)
  
  
  #basically splitting between the somali/white/other minority and weston and non weston ones
  if (length(groups_to_include) > 2){
  
    #this may be a gross implementation but we are relying on positioning in the set up phase to pick outcome and colours
  p <- ggplot(data, aes(x=Month, y=.data[[outcome[1]]])) +
    #intervention and outbreak lines
     geom_vline(xintercept = as.Date("2023-06-01"),
               linetype = "dashed", 
               color = "firebrick", 
               linewidth = 0.7)+
    geom_vline(xintercept = as.Date("2024-03-01"),
               linetype = "dashed", 
               color = "firebrick", 
               linewidth = 0.7)+
    annotate(
      "text",
      x = as.Date("2023-11-01"),
      y = 375,
      label = "MMR escalation response",
      color = "firebrick",
      size=3,
      angle = 0)+
    #Somali, in theory
  geom_point (size=1.5,
              color = colours[1])+
    geom_line (linewidth=0.7,
               color = colours[1])  +
    geom_smooth(method=NULL , se=FALSE,
                linewidth = 0.7,
                linetype =3,
                color = colours[1]) +
    #White, in theory
  geom_point (aes(x=Month, y=.data[[outcome[2]]]),
              size=1.5,
              color = colours[2])+
    geom_line (aes(x=Month, y=.data[[outcome[2]]]),
               linewidth=0.7,
               color = colours[2]) +
    geom_smooth(aes(x=Month, y=.data[[outcome[2]]]),
                method=NULL , se=FALSE,
                linewidth = 0.7,
                linetype =3,
                color = colours [2]) +
    #Other minority, in theory
  geom_point (aes(x=Month, y=.data[[outcome[3]]]),
              size=1.5,
              color = colours[3])+
    geom_line (aes(x=Month, y=.data[[outcome[3]]]),
               linewidth=0.7,
               color = colours[3]) +
    geom_smooth(aes(x=Month, y=IncidencePct_C1_OtherMinority),
                method=NULL , se=FALSE,
                linewidth = 0.7,
                linetype =3,
                color = colours[3]) +
    #y axis and x axis scales
  scale_y_continuous(
    breaks = seq (0,max_value, steps),
    limits = c(0, max_value),
    expand = c(0, 0))  +
    scale_x_date (
      date_breaks = "3 month",       # breaks every month
      date_labels = "%b %Y")  +    # labels like "Jan 2024", "Feb 2024"
    labs(title = cohort_title, 
         #subtitle = "",
         x = "Time - Month",
         y =  "Incidence Rate per 1000 individuals") +
    #theme
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=.1, color="grey90"),
        panel.background = element_rect(fill="white"),
        legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        title = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 13,
                                   angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.line = element_line(linewidth = 0.5, 
                                 colour = "black", 
                                 linetype=1))

  
  } else { #implement with two lines, for Weston/Non-Weston
    p <- ggplot(data, aes(x=Month, y=.data[[outcome[1]]])) +
      #intervention and outbreak lines
      geom_vline(xintercept = as.Date("2023-06-01"),
                 linetype = "dashed", 
                 color = "firebrick", 
                 linewidth = 0.7)+
      geom_vline(xintercept = as.Date("2024-03-01"),
                 linetype = "dashed", 
                 color = "firebrick", 
                 linewidth = 0.7)+
      annotate(
        "text",
        x = as.Date("2023-11-01"),
        y = 375,
        label = "MMR escalation response",
        color = "firebrick",
        size=3,
        angle = 0)+
      #Weston deprived, in theory
      geom_point (size=1.5,
                  color = colours[1])+
      geom_line (linewidth=0.7,
                 color = colours[1])  +
      geom_smooth(method=NULL , se=FALSE,
                  linewidth = 0.7,
                  linetype =3,
                  color = colours[1]) +
      #Non-Weston Deprived, in theory
      geom_point (aes(x=Month, y=.data[[outcome[2]]]),
                  size=1.5,
                  color = colours[2])+
      geom_line (aes(x=Month, y=.data[[outcome[2]]]),
                 linewidth=0.7,
                 color = colours[2]) +
      geom_smooth(aes(x=Month, y=.data[[outcome[2]]]),
                  method=NULL , se=FALSE,
                  linewidth = 0.7,
                  linetype =3,
                  color = colours [2]) +
      #scale y/x axes
      scale_y_continuous(
        breaks = seq (0, max_value, steps),
        limits = c(0, max_value),
        expand = c(0, 0))  +
      scale_x_date (
        date_breaks = "3 month",       # breaks every month
        date_labels = "%b %Y")  +    # labels like "Jan 2024", "Feb 2024"
      labs(title = cohort_title, 
           #subtitle = "",
           x = "Time - Month",
           y =  "Incidence Rate per 1000 individuals") +
      #theme
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size=.1, color="grey90"),
          panel.background = element_rect(fill="white"),
          legend.position = "right",
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          title = element_text(size = 13),
          axis.title = element_text(size = 13),
          axis.text.x = element_text(size = 13,
                                     angle = 45,
                                     hjust = 1),
          axis.text.y = element_text(size = 13),
          axis.line = element_line(linewidth = 0.5, 
                                   colour = "black", 
                                   linetype=1))

  }
  
  return(p)
}


#get working directory from paths.R
setwd(wd)

#read in data, fix dates
vax <- read_excel("Maindata_13_08_26.xlsx", sheet = "BNSSG_Monthly")
vax$Month <- as.Date(as.character(vax$Month))
vax$month_start <- as.Date(vax$month_start)

cohort_spec <- tribble(
  ~cohort, ~groups_to_include, ~cohort_title, ~colours, ~var, ~steps, ~max_value,
  
  "C1",
  list(c("Somali", "White", "OtherMinority")),
  "Between 12 and 47 months (routine)",
  list(c("#5D8AA8",
         "#D99058",
         "#B284BE")),
  "IncidencePct",
  25,
  500,
  

  "C2",
  list(c("Somali", "White", "OtherMinority")),
  "Between 4 and 16 years (children catch-up)",
  list(c("#5D8AA8",
         "#D99058",
         "#B284BE")),
  "IncidencePct",
  10,
  50,

  
  "C3",
  list(c("Somali", "White", "OtherMinority")),
  "Between 17 and 34 years (adult catch-up)",
  list(c("#5D8AA8",
         "#D99058",
         "#B284BE")),
  "IncidencePct",
  2,
  14,
  
  "C1",
  list(c("WestonDeprived", "NonWestonDeprived")),
  "Between 12 and 47 months (routine)",
  list(c("#006D6F",
         "#702963")),
  "IncidencePct",
  25,
  500,
  
  "C2",
  c("WestonDeprived", "NonWestonDeprived"),
  "Between 4 and 16 years (children catch-up)",
  list(c("#006D6F",
         "#702963")),
  "IncidencePct",
  10,
  50,
  
  "C3",
  c("WestonDeprived", "NonWestonDeprived"),
  "Between 17 and 34 years (adult catch-up)",
  list(c("#006D6F",
         "#702963")),
  "IncidencePct",
  2,
  14
)


plot_output <- pmap(
  cohort_spec,
  function(cohort,
           groups_to_include,
           cohort_title,
           colours,
           var,
           steps,
           max_value){
    
    time_series_graph(
      data = vax,
      cohort = cohort,
      cohort_title = cohort_title,
      var = var,
      groups_to_include = groups_to_include,
      steps = steps,
      max_value = max_value,
      colours = colours
    )
    
  }
)

#arrange plots stacked

somali_plot <- ggarrange(
  plotlist = c(plot_output[1],plot_output[2],plot_output[3]),
  ncol = 1,
  nrow = 3
) +
  theme(
    plot.background = element_rect(fill = "white")
  )

weston_plot <- ggarrange(
  plotlist = c(plot_output[4],plot_output[5],plot_output[6]),
  ncol = 1,
  nrow = 3
) +
  theme(
    plot.background = element_rect(fill = "white")
  )

#save plots 

ggsave(
  filename = "somali.png",
  plot = somali_plot,
  width = 4000,
  height = 4000,
  units = "px",
  dpi = 600
)

ggsave(filename = "weston.png",
        plot = weston_plot,
        width = 4000,
        height = 4000,
        units = "px",
        dpi = 600
)
  

