library(sf)
library(dplyr)
library(readxl)
library(purrr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(ggspatial)

#function to make a single map, taking arguments for breaks etc
#focused on incidence
#we use the wrapper function below to call this and make the arguments
make_map <- function(data,
                     var,
                     title,
                     low_col,
                     high_col,
                     max_value,
                     step,
                     map_elements = FALSE) {
  
  p <- ggplot(data) +
    geom_sf(
      aes(fill = .data[[var]]),
      colour = NA
    ) +
    scale_fill_gradient(
      low = low_col,
      high = high_col,
      name = "Rate per 1,000",
      breaks = seq(from = 0, to = max_value, by = step),
      limits = c(0, max_value),
      na.value = "grey90",
      guide = guide_colourbar(
        title.position = "left",
        title.hjust = 0.5,
        barwidth = unit(12, "cm"),
        barheight = unit(0.25, "cm")
      )
    )+
    theme_void() +
    theme(
      plot.title = element_text(
        size = 10,
        face = "plain",
        hjust = 0.5
      ),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    ) +
    labs(title = title)
  
  if (map_elements) {
    
    p <- p +
      annotation_scale(
        location = "br",
        width_hint = 0.3
      ) +
      annotation_north_arrow(
        location = "tr",
        pad_x = unit(0.2, "cm"),
        pad_y = unit(3, "cm"),
        width = unit(0.4, "cm"),
        height = unit(0.4, "cm"),
        which_north = "true",
        style = north_arrow_orienteering
      )
    
  }
  
  return(p)
  
}

#function to go over the cohorts and ethnicites and make maps wtih make_map
#this makes a row of maps by cohort
make_map_row <- function(data,
                         cohort,
                         cohort_title,
                         low_col,
                         high_col,
                         max_value,
                         step,
                         ethnicity_spec) {
  
  map_list <- ethnicity_spec %>%
    mutate(
      var = paste0("mean_incidence_", cohort, "_", ethnicity),
      map_elements = row_number() == 1 #this means the map elements (north arrow/scale bar) will only be on the first map
    ) %>%
    pmap(
      \(ethnicity, plot_title, var, map_elements) {
        
        make_map(
          data = data,
          var = var,
          title = plot_title,
          low_col = low_col,
          high_col = high_col,
          max_value = max_value,
          step = step,
          map_elements = map_elements
        )
        
      }
    )
  
  ggarrange(
    plotlist = map_list,
    ncol = length(map_list),
    nrow = 1,
    common.legend = TRUE,
    legend = "bottom"
  ) +
    labs(title = cohort_title) +
    theme(
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0
      )
    )
}

#set wd from paths.R
setwd(wd)

#load msoa as sf
msoa <- read_sf("Middle_layer_Super_Output_Areas_December_2021_Boundaries.geojson")

#filter to bnssg
bnssg_msoa <- msoa |>
  filter(
    grepl("Bristol", MSOA21NM) |
      grepl("North Somerset", MSOA21NM) |
      grepl("South Gloucestershire", MSOA21NM)
  )

#cleanup
rm(msoa)

#locality look up thanks to Carlos
locality_lookup <- tribble(
  ~locality, ~MSOA21CD,
  "Inner City and East", "E02003039",
  "Inner City and East", "E02003024",
  "Inner City and East", "E02003043",
  "Inner City and East", "E02006888",
  "Inner City and East", "E02006890",
  "Inner City and East", "E02003036",
  "Inner City and East", "E02006889",
  "Inner City and East", "E02003034",
  
  "North and West Bristol", "E02003022",
  "North and West Bristol", "E02003021",
  "North and West Bristol", "E02003013",
  "North and West Bristol", "E02003016",
  "North and West Bristol", "E02003041",
  "North and West Bristol", "E02003012",
  "North and West Bristol", "E02003019",
  "North and West Bristol", "E02003033",
  "North and West Bristol", "E02003017",
  "North and West Bristol", "E02003037",
  
  "South Bristol", "E02003063",
  "South Bristol", "E02003050",
  "South Bristol", "E02003046",
  "South Bristol", "E02003056",
  "South Bristol", "E02003059",
  "South Bristol", "E02003064",
  "South Bristol", "E02003057",
  "South Bristol", "E02003053",
  "South Bristol", "E02003060",
  
  "South Gloucestershire", "E02003119",
  "South Gloucestershire", "E02003104",
  "South Gloucestershire", "E02003100",
  "South Gloucestershire", "E02003095",
  "South Gloucestershire", "E02003111",
  "South Gloucestershire", "E02003110",
  "South Gloucestershire", "E02003105",
  "South Gloucestershire", "E02003094",
  "South Gloucestershire", "E02003120",
  "South Gloucestershire", "E02003116",
  "South Gloucestershire", "E02003112",
  "South Gloucestershire", "E02003015",
  "South Gloucestershire", "E02003117",
  "South Gloucestershire", "E02003018",
  "South Gloucestershire", "E02003090",
  "South Gloucestershire", "E02003091",
  "South Gloucestershire", "E02003113",
  "South Gloucestershire", "E02003096",
  
  "Weston, Worle and Villages", "E02003080",
  "Weston, Worle and Villages", "E02003084",
  "Weston, Worle and Villages", "E02003085",
  "Weston, Worle and Villages", "E02006845",
  "Weston, Worle and Villages", "E02003081",
  "Weston, Worle and Villages", "E02003086",
  "Weston, Worle and Villages", "E02003088",
  
  "Woodspring", "E02003069",
  "Woodspring", "E02003067",
  "Woodspring", "E02003068",
  "Woodspring", "E02003065",
  "Woodspring", "E02003072"
)

#join locality lookup to MSOAs
bnssg_msoa <- bnssg_msoa %>%
  left_join(locality_lookup, by = "MSOA21CD")

#dissolve boundaries to locality (sf format enables this)
bnssg_localities <- bnssg_msoa %>%
    group_by(locality) %>%
    summarise()

#cleanup
rm(bnssg_msoa)

#add in data and join to localities
incidence <- read_excel("somali_comparison_with_%.xlsx", 
                        sheet = "geo")

bnssg_localities <- bnssg_localities %>%
  left_join(incidence, by = "locality")

#set specifications to run functions
#ethnicity one is simple
ethnicity_spec <- tribble(
  ~ethnicity,      ~plot_title,
  "Somali",        "Somali",
  "White",         "White",
  "OtherMinority", "Other Minority"
)

#cohort specification includes colour for colour ramp, low/high colours and breaks
cohort_spec <- tribble(
  ~cohort, ~cohort_title, ~low_col,   ~high_col, ~max_value, ~step,
  
  "C1",
  "Between 12 and 47 months (routine)",
  "#87CEFA",
  "#00009C",
  350,
  50,
  
  "C2",
  "Between 4 and 16 years (children catch-up)",
  "#FFE4E1",
  "#FE0000",
  10,
  2,
  
  "C3",
  "Between 17 and 34 years (adult catch-up)",
  "#FFFDD0",
  "#EF9B0F",
  2,
  0.5
)

#generate cohort maps using the spec
#using pmap from purrr which conveniently loops over functions
#feeding in over cohort_spec
#no need for map element as they're created in the make_map_row wrapper function 
cohort_maps <- cohort_spec %>%
  pmap(
    \(cohort,
      cohort_title,
      low_col,
      high_col,
      max_value,
      step) {
      
      make_map_row(
        data = bnssg_localities,
        cohort = cohort,
        cohort_title = cohort_title,
        low_col = low_col,
        high_col = high_col,
        max_value = max_value,
        step = step,
        ethnicity_spec = ethnicity_spec
      )
      
    }
  )

#arrange all maps into one plot
map_plot <- ggarrange(
  plotlist = cohort_maps,
  ncol = 1,
  nrow = length(cohort_maps)
) +
  theme(
    plot.background = element_rect(fill = "white")
  )

#save map

ggsave(
  filename = "locality_cohort_map.png",
  plot = map_plot,
  width = 4000,
  height = 4000,
  units = "px",
  dpi = 320
)
