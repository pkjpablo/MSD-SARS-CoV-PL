# Antibody lanscape, 204 participants
# Pablo 
# Base on roessler_netzl_et_al2022
#https://github.com/acorg/roessler_netzl_et_al2022/tree/main
setwd("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/1. Landscape clean/MSD")

#library 
library(Racmacs)
library(readxl)
library(dplyr)
library(stringr)

# ------------------------------------------ Dataset  --------------------------------------

long <- read_excel("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/MSD geral/MSD/dataset/bases anteriores/2024 06 05 192 long primiera parte.xlsx")
db2 <- long

# Crear una nueva columna de survey+ID+sample_typr
db2$id <- paste(db2$survey, db2$idnova, db2$grupo2, sep = "_")
db2$id

table(db2$lineage)

db2$lineage <- ifelse(db2$lineage == "ancestral" & db2$panel == 32, "ancestral_p32", db2$lineage)
db2$lineage <- ifelse(db2$lineage == "ancestral" & db2$panel == 23, "D614G", db2$lineage)
db2$lineage <- ifelse(db2$lineage == "BA.1" & db2$panel == 23, "BA.1_p23", db2$lineage)
db2$lineage <- ifelse(db2$lineage == "P.1", "P.1.1", db2$lineage)
db2$lineage <- ifelse(db2$lineage == "B.1.617.2.Seq2", "B.1.617.2", db2$lineage)
db2$lineage <- ifelse(db2$lineage == "BA.2.75", "BA.2", db2$lineage)
#db2$lineage <- ifelse(db2$lineage == "BA.5", "BA.5.3.2", db2$lineage)
db2$lineage <- ifelse(db2$lineage == "BQ.1", "BQ.1.3", db2$lineage)

db2 <- db2[db2$lineage != "ancestral_p32", ]
db2 <- db2[db2$lineage != "BA.1_p23", ]
db2 <- db2[db2$lineage != 'AY.4.2', ]
db2 <- db2[db2$lineage != 'BA.2.75.2', ]
db2 <- db2[db2$lineage != 'BA.4.6', ]
db2 <- db2[db2$lineage != 'B.1.617.Seq1', ]
db2 <- db2[db2$lineage != 'BA.5.3.2', ]


db2$lineage[db2$lineage == "P.1"] <- "P.1.1"
db2$lineage[db2$lineage == "B.1.617.Seq2"] <- "B.1.617.2"
db2$lineage[db2$lineage == "BA.2.75"] <- "BA.2"
db2$lineage[db2$lineage == "BA.5"] <- "BA.5.2.1"
db2$lineage[db2$lineage == "BQ.1"] <- "BQ.1.3"
db2$lineage[db2$lineage == "ancestral"] <- "D614G"

table(db2$lineage)

#

db3 <- db2[,c("id","lineage", "value")]

library(tidyr)
db3 <- db3 %>%
  pivot_wider(names_from = lineage, values_from = value)


# eliminar muestras sin resutlados 
head(db3)


# Transformar valores menores o iguales a 0 a 0.1 en todas las variables del dataframe db3

db3[db3 <= 0] <- 0.01

# Transponer el dataframe db3 en el formato para cartgrafia 
db4 <- db3
#db3 <- db4_subset <- subset(db4, grepl("^L45", id))

db3_transpuesto <- t(db3)
colnames(db3_transpuesto) <- db3_transpuesto[1, ]
db3_transpuesto <- db3_transpuesto[-1, ]

a <- print(db3_transpuesto)

# ------------------------------------------ Roessler map --------------------------------------


map <- read.acmap("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/Antigenic maps/roessler_netzl_et_al2023/roessler_netzl_et_al2023-main/data/maps/alignment_map.ace")


## Select the MSD 
map <- subsetMap(map, antigens = c("D614G",
                                   "B.1.1.7",
                                   "B.1.351",
                                   "P.1.1",
                                   "B.1.617.2",
                                   "BA.1",
                                   "BA.2",
                                   "BA.5.2.1",
                                   "BF.7",
                                   "BQ.1.3",
                                   "BQ.1.1",
                                   "XBB.1"))




full_map_p1_adj <- map
mapColors <- c(
  "D614G" = "#393b79",
  "BF.7" = "#647A39",
  "B.1.617.2" = "#d18652",
  "XBB.1" = "#5B004C",
  "BA.2" = "#5B004C",
  "B.1.1.7" = "#637939",
  "BQ.1.3" = "#CD9B1D",
  "P.1.1" = "#b471ab",
  "BA.1" = "#EF3737",
  "BA.5.2.1" = "#bfdaa0",
  "B.1.351" = "#107f97",
  "BQ.1.1" = "#107f97"
)

# Aplicar los colores
agFill(full_map_p1_adj) <- mapColors[agNames(full_map_p1_adj)]

plot(full_map_p1_adj)


#save.acmap(full_map_p1_adj, "maps/Roessler map.ace")

## ver para salvar esto 

# ------------------------------------------ PL pablo map --------------------------------------

set.seed(104)

# 2. Create the acmap object, specifying the titer table ----
library(Racmacs)
library(readr)
library(titertools)
library(ablandscapes) 
library(r3js)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)

map <- acmap(
  titer_table = a
)

mapColors <- read.csv("dataset/Ejemplo 3/map-colors v2.csv", 
                      comment = ";",row.names = 'Antigen', header = TRUE)
mapColors <- read_csv("dataset/Ejemplo 3/map-colors v2.csv")

sr_group_colors <- read_delim("dataset/Ejemplo 3/sr_group_colors.csv", 
                              delim = ";", 
                              escape_double = FALSE, 
                              trim_ws = TRUE)

#sr_group_colors <- sr_group_colors %>%
#  tibble::column_to_rownames(var = "SerumGroup")


# 3. Ejecutar un conjunto de optimizaciones ---- 
# para intentar encontrar el mejor mapa para representar los datos.

map <- optimizeMap(
  map                     = map,
  number_of_dimensions    = 2,
  number_of_optimizations = 1000,
  minimum_column_basis    = "none"
)

plot(map)


mapColors <- c(
  "D614G" = "#393b79",
  "BF.7" = "#647A39",
  "B.1.617.2" = "#d18652",
  "XBB.1" = "#5B004C",
  "BA.2" = "#5B004C",
  "B.1.1.7" = "#637939",
  "BQ.1.3" = "#CD9B1D",
  "P.1.1" = "#b471ab",
  "BA.1" = "#EF3737",
  "BA.5.2.1" = "#bfdaa0",
  "B.1.351" = "#107f97",
  "BQ.1.1" = "#107f97"
)

# Aplicar los colores
agFill(map) <- mapColors[agNames(map)]

agGroups(map) <- agNames(map)

srNames(map)

plot(map)

# Suponiendo que 'srNames(map)' es tu vector

sr_groups <- sub(".+_(.+)$", "\\1", srNames(map))
#sr_groups <- substr(srNames(map), 1, 3)

srGroups(map) <- factor(sr_groups, levels = sr_group_colors$SerumGroup)
srOutline(map) <- mapColors[as.character(srGroups(map)),]

##map <- apply_style(map)
plot(map)

#save.acmap(map, "maps/pau da lima p23 and 32 L45-49 all.ace")

# ------------------------------------------ map opti  --------------------------------------
#full_map <- read.acmap("maps/pau da lima p23 L45-L47.ace")
full_map = map
plot(full_map)
full_map_p1_adj <- optimizeAgReactivity(full_map, fixed_ag_reactivities = c(0, 0, 0,
                                                                            0, 0, 0,
                                                                            0, 0, 0,
                                                                            0, 0,0), reoptimize = F)
full_map_p1_adj <- optimizeMap(full_map_p1_adj, number_of_dimensions = 2, number_of_optimizations = 1000, options = list(ignore_disconnected = TRUE))
full_map_p1_adj <- realignMap(full_map_p1_adj, full_map)

plot(full_map_p1_adj)

map_pablopl <-  full_map_p1_adj
#save.acmap(full_map_p1_adj, "maps/pau da lima p23 and 32 L45-49 all.ace")

#----------------------------- mapa roessler_netzl_et_al2022
#map_jp <- read.acmap("maps/pau da lima p23 and 32 L45-L48 all.ace")
map_jp <- map_pablopl
plot(map_jp)

map <- read.acmap("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/Antigenic maps/roessler_netzl_et_al2023/roessler_netzl_et_al2023-main/data/maps/alignment_map.ace")
map <- read.acmap("maps/Roessler map.ace")

map <- subsetMap(map, antigens = c("D614G",
                                   "B.1.1.7",
                                   "B.1.351",
                                   "P.1.1",
                                   "B.1.617.2",
                                   "BA.1",
                                   "BA.2",
                                   "BA.5.2.1",
                                   "BF.7",
                                   "BQ.1.3",
                                   "BQ.1.1",
                                   "XBB.1"))



full_map_p1_adj <- map
plot(full_map_p1_adj)





# ------------------------------------------ landscape --------------------------------------

### Use the version 4.2.2 

#setup page and load metadata
#rm(list = ls())
library(Racmacs)
library(tidyverse)
library(meantiter)
library(titertools)
library(ablandscapes) 
library(r3js)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)

set.seed(100)

source("./functions/remove_reactivity_bias.R")
source("./functions/map_longinfo.R")
source("./functions/sams_landscape_functions.R")

sr_group_colors <- read_delim("dataset/Ejemplo 3/sr_group_colors.csv", 
                              delim = ";", 
                              escape_double = FALSE, 
                              trim_ws = TRUE)

sr_group_colors <- sr_group_colors %>%
  tibble::column_to_rownames(var = "SerumGroup")

mapColors <- read.csv("dataset/Ejemplo 3/map-colors v2.csv", 
                      comment = ";",row.names = 'Antigen', header = TRUE)

#sr_group_colors <- read_delim("dataset/Ejemplo 3/sr_group_colors.csv", 
#                             delim = ";", 
#                            escape_double = FALSE, 
#                           trim_ws = TRUE)


# read the full map

map_orig = full_map_p1_adj
plot(map_orig)


sr_colors <- read_delim("dataset/Ejemplo 3/sr_group_colors.csv", 
                        delim = ";", 
                        escape_double = FALSE, 
                        trim_ws = TRUE)
sr_colors <- sr_colors %>%
  tibble::column_to_rownames(var = "SerumGroup")



map_jp <- read.acmap("maps/pau da lima p23 and 32 L45-49 all.ace")
plot(map_jp)
map_long <- long_map_info(map_jp)



#map_long <-map_long[grep("^L45", map_long$sr_name), ]
#map_long <-map_long[grep("^L46", map_long$sr_name), ]
#map_long <-map_long[grep("^L47", map_long$sr_name), ]
#map_long <-map_long[grep("^L48", map_long$sr_name), ]

head(map_long)


table(map_long$sr_name)
hist(as.numeric(map_long$titer))

library(dplyr)

## exponencial de base 2
map_long <- map_long %>%
  mutate(titer = 10*2^(as.numeric(titer)*10))
head(map_long)
map_long %>% ggplot(aes(x=titer))+geom_histogram()+scale_x_continuous(trans = "log")
map_long$titer_adjusted <- as.numeric(map_long$titer_adjusted)
hist(map_long$titer_adjusted)


map_long %>%
  select(titer, ag_name, sr_name, sr_group) -> titerdata


titerdata %>%
  group_by(
    sr_group
  ) -> titerdata

#titerdata$sr_group <- "L45"
titerdata$sr_group <- substr(titerdata$sr_name, 1, 3)


titerdata %>%
  group_map(
    get_titertable
  ) -> titertables

lndscp_fits <- lapply(
  titertables,
  function(titertable) {
    
    ablandscape.fit(
      titers = titertable, #[,ags_to_fit_lndscp],
      bandwidth = 1,
      degree = 1,
      method = "cone",  ### cone or loess
      error.sd = 1,
      acmap = full_map_p1_adj,
      control = list(
        optimise.cone.slope = TRUE
      )
    )
  }
)


plot(full_map_p1_adj)

titertables_groups <- group_data(titerdata)

titertables_groups$sr_group
titertables_groups$.rows



# Add impulses


titerdata %>%
  group_by(
    sr_group,
    ag_name
  ) %>%
  dplyr::summarize(gmt = mean(titer, na.rm = T)/1000)-> gmt_data

titerdata %>%
  group_by(
    sr_group,
    ag_name
  ) %>%
  dplyr::summarize(gmt = titertools::gmt(titer, dilution_stepsize = 0)["mean", "estimate"])-> gmt_data


# angle for html page
angle <- list(
  rotation = c(-1.4427, 0.0100, -0.0263), #c(-1.3365, 0.0055, -0.0576),# c(-1.4592, 0.0045, -0.0144)
  translation = c(0, 0,0), #translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.5
  # zoom = 1.1646 # higher is more zoomed out
)


#Modified function for base plot
base_plot_data3js_grs <- function(map, lndscp_fits, highlighted_ags, lims, ag_plot_names, alternative_ba5 = FALSE, opti_nr = 1,
                                  add_border = TRUE, add_axis = TRUE){
  
  x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
  y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
  z_coords <- rep(0.02, length(highlighted_ags))
  ag_point_size <- c(rep(14, length(highlighted_ags))) / 5
  ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags])
  ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags])
  labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags])
  border_col <- "grey50"
  
  z_lims <- c(0,10)
  axis_at <- seq(z_lims[1], z_lims[2],2)
  # Setup plot
  data3js <- ablandscapes:::lndscp3d_setup(
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = z_lims,
    aspect.z = 0.5,
    options = list(
      lwd.grid =  0.05,
      sidegrid.lwd = 1,
      sidegrid.col = border_col,
      sidegrid.at = list("z" = axis_at),
      zaxt = "lin"
    ),
    show.axis = FALSE
  )
  
  if(add_axis){
    
    axis_labels <- seq(0,1,by=0.2)
    
    data3js <- r3js::axis3js(
      data3js,
      side = "z",
      at = axis_at,
      labels = axis_labels,
      # labeloffset = 0.11,
      cornerside = "f",
      size = 20,
      alignment = "right"
    )
  }
  
  # Add basemap
  data3js <- lndscp3d_map(
    data3js = data3js,
    fit = lndscp_fits[[1]],
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = c(0, 10),
    show.map.sera = FALSE,
    options = list(
      opacity.basemap = 0.3
    )
  )
  
  data3js <- r3js::points3js(
    data3js,
    x          = x_coords,
    y          = y_coords,
    z          = z_coords,
    size       = ag_point_size,
    col        = ag_col,
    fill       = ag_fill,
    lwd        = 0.5,
    opacity    = 1,
    highlight  = list(col = "red"),
    label      = labels,
    toggle     = "Basepoints",
    depthWrite = FALSE,
    shape      = "circle filled"
  )
  
  if(add_border){
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[1]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[2],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    # y border
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[1]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[2], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    data3js <- r3js::box3js(
      data3js,
      col   = border_col
    )
    
  }
  
  return(data3js)
}


#Base plot
lndscp_list <- list()
data3js <- base_plot_data3js_grs(full_map_p1_adj, lndscp_fits, agNames(full_map_p1_adj), lims = list(xlim=c(-2,6),ylim=c(-3,2)), agNames(full_map_p1_adj), add_axis = T)
data3js
#str(data3js)

titertables_groups


land_col <- data.frame(
  Color = c("#FF4040", "#1E90FF", "#FFB90F", "#00688B", "#6A5ACD"), 
  row.names = c("L45", "L46", "L47", "L48", "L49")
)

# do plot with L46 and L47 NVSN SS1 SS2 RE V SV
L45_49 <- plot_landscapes_from_list(data3js,
                                    titertables_groups[c(1:5),],
                                    lndscp_fits[c(1:5)],
                                    full_map_p1_adj,
                                    gmt_data, agNames(full_map_p1_adj),
                                    agNames(full_map_p1_adj),
                                    lndscp_colors = land_col,
                                    show_gmts = F)
L45_49


L45 <- plot_landscapes_from_list(data3js, titertables_groups[c(1:1),], lndscp_fits[c(1:1)], map, gmt_data, agNames(map), agNames(map), lndscp_colors = land_col,show_gmts = F)
L45

L46 <- plot_landscapes_from_list(data3js, titertables_groups[c(2:2),], lndscp_fits[c(2:2)], map, gmt_data, agNames(map), agNames(map), lndscp_colors = land_col,show_gmts = F)
L46

L47 <- plot_landscapes_from_list(data3js, titertables_groups[c(3:3),], lndscp_fits[c(3:3)], map, gmt_data, agNames(map), agNames(map), lndscp_colors = land_col,show_gmts = F)
L47

L48 <- plot_landscapes_from_list(data3js, titertables_groups[c(4:4),], lndscp_fits[c(4:4)], map, gmt_data, agNames(map), agNames(map), lndscp_colors = land_col,show_gmts = F)
L48

L49 <- plot_landscapes_from_list(data3js, titertables_groups[c(5:5),], lndscp_fits[c(5:5)], map, gmt_data, agNames(map), agNames(map), lndscp_colors = land_col,show_gmts = F)
L49




panel_plot_grid <- (L45 + L46 + L47 + L48 + L49) + plot_layout(ncol = 2)

print(panel_plot_grid)


### extraer info del landscape
summary(lndscp_fits[c(1:1)])
lndscp_fits[c(1:1)]


result <- sapply(L45_49[["plot"]], function(elem) {
  !is.null(elem$x) && !is.null(elem$y) && !is.null(elem$z)
})

which(result)

lndscp_fits[[1]][["fitted.values"]]
lndscp_fits[[2]][["fitted.values"]]
lndscp_fits[[3]][["fitted.values"]]
lndscp_fits[[4]][["fitted.values"]]
lndscp_fits[[5]][["fitted.values"]]

L45_49[["plot"]][[176]][["z"]]
L45_49[["plot"]][[238]][["z"]]
L45_49[["plot"]][[300]][["z"]]
L45_49[["plot"]][[362]][["z"]]
L45_49[["plot"]][[424]][["z"]]

Racmacs::agCoords(map)

lndscp_fits[[1]][["fitted.values"]]


lndscp_fits[[1]][["cone"]][["cone_heights"]]
lndscp_fits[[2]][["cone"]][["cone_heights"]]
lndscp_fits[[3]][["cone"]][["cone_heights"]]
lndscp_fits[[4]][["cone"]][["cone_heights"]]
lndscp_fits[[5]][["cone"]][["cone_heights"]]

#install.packages("rgl")

library(rgl)

## Cone slope
a <- lndscp_fits[[1]]
lndscp_fits[[1]][["cone"]][["cone_slope"]]
lndscp_fits[[2]][["cone"]][["cone_slope"]]
lndscp_fits[[3]][["cone"]][["cone_slope"]]
lndscp_fits[[4]][["cone"]][["cone_slope"]]
lndscp_fits[[5]][["cone"]][["cone_slope"]]


library(gtsummary)

# Crear el data frame con los valores extraídos
df_slope <- data.frame(
  Survey = paste("Survey", 1:5),
  Cone_Slope = c(
    lndscp_fits[[1]][["cone"]][["cone_slope"]],
    lndscp_fits[[2]][["cone"]][["cone_slope"]],
    lndscp_fits[[3]][["cone"]][["cone_slope"]],
    lndscp_fits[[4]][["cone"]][["cone_slope"]],
    lndscp_fits[[5]][["cone"]][["cone_slope"]]
  )
)

library(knitr)

# Print df_slope as a table
kable(df_slope, col.names = c("Survey", "Cone Slope"), caption = "Cone Slope by Survey")


L45_49[["plot"]][[176]][["z"]]
L45_49[["plot"]][[238]][["z"]]
L45_49[["plot"]][[300]][["z"]]
L45_49[["plot"]][[362]][["z"]]
L45_49[["plot"]][[424]][["z"]]


### Calculando el vertice 
# Suponiendo que `cone$cone_coords` tiene las columnas x e y, y `cone$cone_heights` tiene los valores z:
x <- L45_49[["plot"]][[176]][["x"]]  # Primera columna (X)
y <- L45_49[["plot"]][[176]][["y"]]  # Segunda columna (Y)
z1 <- L45_49[["plot"]][[176]][["z"]]     # Alturas (Z)
z2 <- L45_49[["plot"]][[238]][["z"]]
z3 <- L45_49[["plot"]][[300]][["z"]]
z4 <- L45_49[["plot"]][[362]][["z"]]
z5 <- L45_49[["plot"]][[424]][["z"]]

# Matrices de entrada: x, y, z1
# Calculamos el tamaño de las celdas en las direcciones x e y
dx <- diff(x[, 1])  # Tamaño de celda en x
dy <- diff(y[1, ])  # Tamaño de celda en y



# Encontrar los vértices de los 5 conos
vertices <- data.frame(
  Antigen = paste0("vertice", 1:5),
  x = c(x[which.max(z1)], x[which.max(z2)], x[which.max(z3)], x[which.max(z4)], x[which.max(z5)]),
  y = c(y[which.max(z1)], y[which.max(z2)], y[which.max(z3)], y[which.max(z4)], y[which.max(z5)]),
  Color = "Black"
)

print(vertices)

# Convertir las coordenadas de los antígenos en un data frame
antigen_coords <- as.data.frame(Racmacs::agCoords(full_map_p1_adj))
antigen_coords$Antigen <- rownames(antigen_coords)

# Convertir los colores en un data frame
mapColors2 <- as.data.frame(mapColors)
mapColors2$Antigen <- rownames(mapColors2)

mapColors2 <- mapColors2 %>%
  mutate(Color = case_when(
    Antigen == "BA.5.2.1" ~ "#bfdaa0",
    Antigen %in% c("BQ.1.3", "BQ.1.1") ~ "#107f97",
    TRUE ~ Color
  ))


# Merge de datos base
merged_df <- merge(antigen_coords, mapColors2, by = "Antigen", all.x = TRUE)
colnames(merged_df)[2:3] <- c("x", "y")

vertices$x[1]

# Lista de vértices con coordenadas
vertices <- list(
  vertice1 = c(vertices$x[1], vertices$y[1]),
  vertice2 = c(vertices$x[2], vertices$y[2]),
  vertice3 = c(vertices$x[3], vertices$y[3]),
  vertice4 = c(vertices$x[4], vertices$y[4]),
  vertice5 = c(vertices$x[5], vertices$y[5])
)

# Crear gráficos individuales
plots <- list()

for (vertice in names(vertices)) {
  # Agregar el vértice al dataframe
  temp_df <- rbind(merged_df, data.frame(Antigen = vertice, x = vertices[[vertice]][1], y = vertices[[vertice]][2], Color = "Black"))
  
  # Definir la forma del punto
  temp_df$shape <- ifelse(temp_df$Antigen == vertice, "x", "circle")
  
  # Generar el gráfico
  p <- ggplot(temp_df, aes(x = x, y = y, fill = Color)) +
    geom_point(aes(shape = shape), color = "black", size = 12, stroke = 2) +
    geom_text(aes(label = Antigen), vjust = -3.2, size = 3.5) +
    scale_shape_manual(values = c("circle" = 21, "x" = 4)) +
    scale_fill_identity() +
    theme_bw(base_line_size = 1,base_rect_size = 2) +
    theme(legend.position = "none") +
    labs(title = vertice, x = "", y = "")
  
  # Guardar en la lista de gráficos
  plots[[vertice]] <- p
}

# Mostrar gráficos
plots$vertice1
plots$vertice2
plots$vertice3
plots$vertice4
plots$vertice5


cowplot::plot_grid(plots$vertice1,
                   plots$vertice2,
                   plots$vertice3,
                   plots$vertice4,
                   plots$vertice5, nrow = 1)




######## Calcualr volumenes e slope
result <- sapply(L45_49[["plot"]], function(elem) {
  !is.null(elem$x) && !is.null(elem$y) && !is.null(elem$z)
})
unique_indices <- which(result)[seq(1, length(which(result)), by = 2)]


result <- sapply(L45_49[["plot"]], function(elem) {
  !is.null(elem$x) && !is.null(elem$y) && !is.null(elem$z)
})
unique_indices <- which(result)[seq(1, length(which(result)), by = 2)]



# Definir la función para calcular el volumen, considerando z negativo como 0
calcular_volumen <- function(index) {
  x <- L45_49[["plot"]][[unique_indices[index]]][["x"]]
  y <- L45_49[["plot"]][[unique_indices[index]]][["y"]]
  z <- L45_49[["plot"]][[unique_indices[index]]][["z"]]
  
  # Establecer a 0 todos los valores negativos de z
  #z[z < 0] <- 0
  
  # Convertir x, y, z a matrices de 13x10
  x_matrix <- matrix(x, nrow = 10, ncol = 13, byrow = TRUE)
  y_matrix <- matrix(y, nrow = 10, ncol = 13, byrow = TRUE)
  z_matrix <- matrix(z, nrow = 10, ncol = 13, byrow = TRUE)
  
  # Calcular el área de cada celda
  dx <- abs(x_matrix[1, 1] - x_matrix[1, 2])  # Distancia en x
  dy <- abs(y_matrix[1, 1] - y_matrix[2, 1])  # Distancia en y
  
  # Calcular el volumen usando la regla del trapecio
  volumen <- 0
  
  for (i in 1:(nrow(z_matrix) - 1)) {
    for (j in 1:(ncol(z_matrix) - 1)) {
      # Calcular el volumen de cada celda
      volume_cell <- (z_matrix[i, j] + z_matrix[i + 1, j] + z_matrix[i, j + 1] + z_matrix[i + 1, j + 1]) / 4 * dx * dy
      volumen <- volumen + volume_cell
    }
  }
  
  return(volumen)
}

# Aplicar la función a todos los índices en unique_indices, con valor absoluto si se desea
volumenes <- sapply(1:length(unique_indices), function(i) calcular_volumen(i))

# Imprimir todos los volúmenes
print(volumenes) ### 270 es el 100% del volumen

# Volúmenes correspondientes a cada encuesta
volumenes <- c(21.67979 , 65.26461, 175.38968, 191.66465, 205.88128)

# Agregar los volúmenes a la tabla df_slope
df_slope$Volume <- volumenes
df_slope$Volume100 <- volumenes/270

# If df_slope now has 3 columns (e.g., Survey, Cone Slope, Volume), update the col.names argument
kable(df_slope, col.names = c("Survey", "Cone Slope", "Volume"), caption = "Cone Slope and Volume by Survey")

dates <- as.Date(c("2021-01-07", "2021-08-26", "2022-06-12", "2023-02-11", "2024-01-08"))
df_slope$dates <- dates

df_slope
#surveys dates 
vertical_lines <- as.Date(c("2019-09-09", "2019-11-09", "2020-11-18", "2021-02-26", 
                            "2021-07-14", "2021-10-09", "2022-03-25", "2022-08-31", 
                            "2022-11-16", "2023-05-09", "2023-10-25","2024-03-24"))

# Crear la figura con colores pastel y líneas verticales
Fig3 <- ggplot(df_slope) +
  # Eje para el volumen con color pastel
  geom_line(aes(x = dates, y = Volume100*100), color = "#A3C4BC", size = 1) +  # Pastel verde-azulado
  geom_point(aes(x = dates, y = Volume100*100), color = "#A3C4BC", size = 3) +
  # Eje para el slope con color pastel
  geom_line(aes(x = dates, y = Cone_Slope  * 100), color = "#F7B7A3", size = 1) +  # Pastel rosado
  geom_point(aes(x = dates, y = Cone_Slope  * 100), color = "#F7B7A3", size = 3) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%m-%Y",
               limits = as.Date(c('2020-03-01', '2024-04-01')))+
  # Añadir líneas verticales
  lapply(vertical_lines, function(x) {
    geom_vline(xintercept = as.numeric(x), linetype = "dashed", color = "#82b0d2", size = 1)
  }) +
  scale_y_continuous(
    name = "Volumen (%)",
    sec.axis = sec_axis(~ . / 100, name = "Slope")
  ) +
  # Etiquetas y título
  labs(
    x = "",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )+theme_bw()+theme_classic()

Fig3



# Crear la figura con colores pastel y áreas de fondo marcando los periodos
Fig3a <- ggplot(df_slope) +
  # Crear áreas coloreadas para los períodos definidos por las fechas de vertical_lines
  geom_rect(data = data.frame(
    xmin = as.Date(vertical_lines[seq(1, length(vertical_lines) - 1, by = 2)]),
    xmax = as.Date(vertical_lines[seq(2, length(vertical_lines), by = 2)]),
    ymin = -Inf,
    ymax = Inf
  ), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
  fill = "#D3E8E4", alpha = 0.5) +  # Color pastel con transparencia
  
  # Eje para el volumen con color pastel
  geom_line(aes(x = dates, y = Volume100*100), color = "#A3C4BC", size = 1) +  # Pastel verde-azulado
  geom_point(aes(x = dates, y = Volume100*100), color = "#A3C4BC", size = 3) +
  
  # Eje para el slope con color pastel
  geom_line(aes(x = dates, y = Cone_Slope * 100), color = "#F7B7A3", size = 1) +  # Pastel rosado
  geom_point(aes(x = dates, y = Cone_Slope * 100), color = "#F7B7A3", size = 3) +
  
  scale_x_date(date_breaks = "3 months", date_labels = "%m-%Y",
               limits = as.Date(c('2020-03-01', '2024-04-01'))) +
  
  scale_y_continuous(
    name = "Volumen (%)",
    sec.axis = sec_axis(~ . / 100, name = "Slope")
  ) +
  # Etiquetas y título
  labs(
    x = "",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),  # Cambiar el ángulo a 60 grados y ajustar la posición
    plot.title = element_text(hjust = 0.5)
  ) + 
  theme_bw() + theme_classic()

Fig3a
