library(sf)
library(dplyr)
library(readxl)

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



