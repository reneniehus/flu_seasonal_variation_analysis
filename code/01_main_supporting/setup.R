# setup.R
#
# Session bootstrap for every pipeline script: quiet library loading, dplyr/lubridate remasking,
# repo-wide options (incl. the default ggplot discrete palette), and the small set of shared
# helpers that survive actual use (get_in_hms, replace_inf, EU_short + the country name tables).
# Keep this LEAN: helpers belong here only once more than one script uses them.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Libraries (quiet) ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library <- function(...) suppressPackageStartupMessages(base::library(...,quietly=TRUE))

# libraries
start_time <- Sys.time()

library(magrittr) # better pipes
library(tidyverse) 
library(scales)
library(readxl)
library(purrr)
library(tidyr)
library(readr)
#library(rethinking)
library(patchwork)
library(viridis)
library(here)
# library(dagitty)
library(ISOweek)
library(fst)
#library( "EcdcColors" )
#library(cmdstanr)
library(crayon)
library(lubridate)
# library(Hmisc)
# library(data.table)
# library(dtplyr)

# libraries not to include 
# library(tsibble)

# remasking
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
date <- lubridate::date
intersect <- base::intersect
setdiff <- base::setdiff
union <- base::union
#rstudent <- rethinking::rstudent
expand <- tidyr::expand
map <- purrr::map
discard <- purrr::discard
col_factor <- readr::col_factor
area <- patchwork::area
#compare <- rethinking::compare

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Super basic functions ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
get_in_hms <- function(time1, time2) {
  format(as.POSIXct(as.numeric(difftime(time1, time2, units = 'secs')), 
                    origin = '1970-01-01', tz = 'UTC'), '%H:%M:%S')  
}

replace_inf = function(vector_with_inf,replacement=NA){
  vector_with_inf[is.infinite(vector_with_inf)] = replacement
  return(vector_with_inf)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Other Options ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Reset R's most annoying default options
options(stringsAsFactors = FALSE, 
        scipen = 999, 
        dplyr.summarise.inform = FALSE,
        tibble.print_min=4)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Settings for plotting ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Default DISCRETE colour palette (Dark2) for every ggplot in scripts that source setup.R --
# set the idiomatic way (options), NOT by shadowing ggplot(). The old ggplot() wrapper appended
# scale_color_brewer() to every plot, which silently restyled plots and made any script-level
# scale_colour_* call emit 'Adding another scale for colour...'. An options default is
# overridable per plot without noise, and scripts that never source setup.R behave the same.
options(ggplot2.discrete.colour = function(...) ggplot2::scale_colour_brewer(..., palette = "Dark2"))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Medium-complex functions ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
countries <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", 
               "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", 
               "Hungary", "Iceland", "Ireland", "Italy", "Latvia", "Liechtenstein", 
               "Lithuania", "Luxembourg", "Malta", "Netherlands", "Norway", 
               "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", 
               "Sweden",
               # non EU/EEA
               "Switzerland","England","Northern Ireland","Scotland","EU/EEA")
countries_short <- c("AT", "BE", "BG", "HR", "CY", "CZ", 
                     "DK", "EE", "FI", "FR", "DE", "GR", 
                     "HU", "IS", "IE", "IT", "LV", "LI", 
                     "LT", "LU", "MT", "NL", "NO", 
                     "PL", "PT", "RO", "SK", "SI", "ES", 
                     "SE",
                     # non EU/EEA
                     "CH","GB-ENG","GB-NIR","GB-SCT","EU/EEA")


# EL 
#EU_short("Greece") <- "EL"
# EU_short("Greece","EL")
EU_short <- function(name_long,greece="GR" # or "EL
){
  name_short = name_long
  for (i in 1:length(name_long)) {
    name_short[i] <- countries_short[which(countries%in%name_long[i])]
    if (name_long[i]=="Greece"&greece=="GR") name_short[i]<-"GR"
    if (name_long[i]=="Greece"&greece!="GR") name_short[i]<-"EL"
  }
  
  return(name_short)
}
# inverse mapping (short -> long); used by the compartmental model's data builder, whose contact
# matrices and population pyramids are keyed by long country names
EU_long <- function(name_short, greece="GR"){
  name_long = name_short
  for (i in 1:length(name_long)) {
    if (name_short[i]=="EL"&greece=="EL") name_short[i]<-"GR"
    name_long[i] <- countries[which(countries_short%in%name_short[i])]
  }
  return(name_long)
}
# very end: timing
end_time <- Sys.time()
pr=paste("> Setup script run:",round(end_time - start_time,2),"sec \n"); cat(green(pr))
