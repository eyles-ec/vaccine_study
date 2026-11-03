library(sf) #spatial framework
library(dplyr) #data manipulation
library(readxl) #read in excel files 
library(purrr) #pmap to run over functions
library(tibble) #make tibbles
library(ggplot2) #plot
library(ggpubr) #map publication stuff
library(ggspatial) #use ggplot framework with spatial

#function to make a single map, taking arguments for breaks etc
#focused on incidence
#we use the wrapper function below to call this and make the arguments
make_map <- function(data,
                     var,
                     title,
                     colours,
                     max_value,
                     step,
                     na_colour = "grey70",
                     discrete_scale = FALSE,
                     map_elements = FALSE) {
  
  #because some data have different distribution, may want a discrete scale
  #first coerce colours out of a list into a vector just in case, since the way we define
  #colours in the cohort_spec is in list form
  
  if (is.list(colours)) {
    colours <- colours[[1]]
  }
  
  fill_scale <- if (discrete_scale) {
    
      #will go along with discrete colours rather than a gradient 
      #this means you can choose which colours for each step
      scale_fill_stepsn(
      colours = colours,
      name = "Rate per 1,000",
      breaks = seq(0, max_value, step),
      limits = c(0, max_value),
      na.value = na_colour,
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(12, "cm"),
        barheight = unit(0.25, "cm")
      )
    )
    
  } else {
    
    #what was implemented before with low/high colour 
    #automatically fills in the steps
    scale_fill_gradient(
      low = colours[1],
      high = colours[2],
      name = "Rate per 1,000",
      breaks = seq(0, max_value, step),
      limits = c(0, max_value),
      na.value = na_colour,
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(12, "cm"),
        barheight = unit(0.25, "cm")
      )
    )
    
  }
  
  
  p <- ggplot(data) +
    geom_sf(
      aes(fill = .data[[var]]),
      colour = "grey40"
    ) +
    fill_scale   +
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
        pad_y = unit(2.3, "cm"),
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
                         colours,  
                         max_value,
                         step,
                         discrete_scale,
                         ethnicity_spec 
                         #na_colour
                         ) {
  
  map_list <- ethnicity_spec %>%
    mutate(
      var = paste0("mean_incidence_", cohort, "_", ethnicity),
      map_elements = row_number() == 1 & cohort == "C1" #this means the map elements (north arrow/scale bar) will only be on the first map of C1
    ) %>%
    pmap(
      \(ethnicity, plot_title, var, map_elements) {
        
        make_map(
          data = data,
          var = var,
          title = plot_title,
          colours = colours,
          max_value = max_value,
          step = step,
          #na_colour = na_colour,
          discrete_scale = discrete_scale,
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
        size = 10,
        face = "bold",
        hjust = 0.5
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

#locality lookup, added missing areas to make locality map work
#previous lookup was only where GPs existed 
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
  "Inner City and East", "E02003023",
  "Inner City and East", "E02003025",
  "Inner City and East", "E02003027",
  "Inner City and East", "E02003028",
  "Inner City and East", "E02003029",
  "Inner City and East", "E02003030",
  "Inner City and East", "E02003031",
  "Inner City and East", "E02003032",
  "Inner City and East", "E02003038",
  "Inner City and East", "E02003040",
  "Inner City and East", "E02003044",
  "Inner City and East", "E02003045",
  "Inner City and East", "E02006887",
  "Inner City and East", "E02006949",
  "Inner City and East", "E02006950",
  "Inner City and East", "E02006951",
  "Inner City and East", "E02006952",
  
  "North and West Bristol", "E02003022",
  "North and West Bristol", "E02003021",
  "North and West Bristol", "E02003026",
  "North and West Bristol", "E02003015",
  "North and West Bristol", "E02003013",
  "North and West Bristol", "E02003016",
  "North and West Bristol", "E02003041",
  "North and West Bristol", "E02003012",
  "North and West Bristol", "E02003019",
  "North and West Bristol", "E02003033",
  "North and West Bristol", "E02003017",
  "North and West Bristol", "E02003037",
  "North and West Bristol", "E02003018",
  "North and West Bristol", "E02003014",
  "North and West Bristol", "E02003026",
  "North and West Bristol", "E02003020",

  "South Bristol", "E02003063",
  "South Bristol", "E02003062",
  "South Bristol", "E02003050",
  "South Bristol", "E02003051",
  "South Bristol", "E02003052",
  "South Bristol", "E02003053",
  "South Bristol", "E02003046",
  "South Bristol", "E02003047",
  "South Bristol", "E02003048",
  "South Bristol", "E02003049",
  "South Bristol", "E02003056",
  "South Bristol", "E02003055",
  "South Bristol", "E02003054",
  "South Bristol", "E02003059",
  "South Bristol", "E02003064",
  "South Bristol", "E02003057",
  "South Bristol", "E02003058",
  "South Bristol", "E02003053",
  "South Bristol", "E02003060",
  "South Bristol", "E02003061",
  
  "South Gloucestershire", "E02003119",
  "South Gloucestershire", "E02003104",
  "South Gloucestershire", "E02003100",
  "South Gloucestershire", "E02003095",
  "South Gloucestershire", "E02003111",
  "South Gloucestershire", "E02003110",
  "South Gloucestershire", "E02003094",
  "South Gloucestershire", "E02003120",
  "South Gloucestershire", "E02003116",
  "South Gloucestershire", "E02003112",
  "South Gloucestershire", "E02003117",
  "South Gloucestershire", "E02003090",
  "South Gloucestershire", "E02003091",
  "South Gloucestershire", "E02003113",
  "South Gloucestershire", "E02003096",
  "South Gloucestershire", "E02003092",
  "South Gloucestershire", "E02003093",
  "South Gloucestershire", "E02003097",
  "South Gloucestershire", "E02003098",
  "South Gloucestershire", "E02003099",
  "South Gloucestershire", "E02003101",
  "South Gloucestershire", "E02003102",
  "South Gloucestershire", "E02003103",
  "South Gloucestershire", "E02003107",
  "South Gloucestershire", "E02003108",
  "South Gloucestershire", "E02003109",
  "South Gloucestershire", "E02003114",
  "South Gloucestershire", "E02003115",
  "South Gloucestershire", "E02003118",
  "South Gloucestershire", "E02003121",
  "South Gloucestershire", "E02007064",
  "South Gloucestershire", "E02007065",
  "South Gloucestershire", "E02007066",
  "South Gloucestershire", "E02007067",
  "South Gloucestershire", "E02003105",
  
  "Weston, Worle and Villages", "E02003080",
  "Weston, Worle and Villages", "E02003084",
  "Weston, Worle and Villages", "E02003085",
  "Weston, Worle and Villages", "E02006845",
  "Weston, Worle and Villages", "E02006846",
  "Weston, Worle and Villages", "E02003081",
  "Weston, Worle and Villages", "E02003086",
  "Weston, Worle and Villages", "E02003088", #this MSOA also covers part of Woodspring
  "Weston, Worle and Villages", "E02003089",
  "Weston, Worle and Villages", "E02003087",
  "Weston, Worle and Villages", "E02003078", #this MSOA also covers part of Woodspring
  "Weston, Worle and Villages", "E02003079",
  "Weston, Worle and Villages", "E02003082",
  
  "Woodspring", "E02003069",
  "Woodspring", "E02003067",
  "Woodspring", "E02003068",
  "Woodspring", "E02003065",
  "Woodspring", "E02003072",
  "Woodspring", "E02003066",
  "Woodspring", "E02003070",
  "Woodspring", "E02003071",
  "Woodspring", "E02003073",
  "Woodspring", "E02003074",
  "Woodspring", "E02003075",
  "Woodspring", "E02003076",
  "Woodspring", "E02003077"
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
  ~cohort, ~cohort_title, ~colours, ~max_value, ~step, ~discrete_scale,
  
  "C1",
  "Between 12 and 47 months (routine)",
  list(c("#deebf7",
    "#9ecae1",
    "#6baed6",
    "#3182bd",
    "#08519c")),
  350,
  25,
  TRUE,
  
  "C2",
  "Between 4 and 16 years (children catch-up)",
  list(c("#FFE5d9",
    "#cb181d")),
  10,
  2,
  FALSE,
  
  "C3",
  "Between 17 and 34 years (adult catch-up)",
  list(c("#Fef0d9",
  "#e69f00")),
  2,
  0.5,
  FALSE
)

#generate cohort maps using the spec
#using pmap from purrr which conveniently loops over functions
#feeding in over cohort_spec
#no need for map element as they're created in the make_map_row wrapper function 
cohort_maps <- cohort_spec %>%
  pmap(
    \(cohort,
      cohort_title,
      colours,
      max_value,
      step,
      discrete_scale) {
      
      make_map_row(
        data = bnssg_localities,
        cohort = cohort,
        cohort_title = cohort_title,
        colours = colours,
        max_value = max_value,
        step = step,
        discrete_scale = discrete_scale,
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
  dpi = 600
)
